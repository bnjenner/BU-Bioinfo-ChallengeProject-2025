import os                                                                                                                        
import sys      
import argparse
import pysam
import pandas as pd


def load_bed(bed_path):
    """
    Load splice site centers from a BED file.
    Column 1 (chromStart, 0-based) is used as the site center.
    Column 3 (name) is used as the label if present, otherwise 'chrom:start'.

    Note: BED chromStart is 0-based; BAM positions are converted to 1-based
    (+1). Add 1 here to keep coordinates consistent.
    """
    bed = pd.read_csv(bed_path, sep="\t", header=None)
    has_name = bed.shape[1] >= 4

    splices = pd.DataFrame({
        "splice_name":   bed[3] if has_name else bed[0].astype(str) + ":" + bed[1].astype(str),
        "splice_center": bed[1] + 1,   # BED 0-based -> 1-based to match BAM pos
    })
    return splices


def bam_to_start_counts(bam_path):
    """Count positive-strand read starts per reference position."""
    counts = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if read.is_unmapped or read.is_reverse:
                continue
            pos = read.reference_start + 1  # convert to 1-based
            counts[pos] = counts.get(pos, 0) + 1

    return pd.DataFrame(sorted(counts.items()), columns=["pos", "N"])


def label_splice_sites(start_pos, splices, window):
    splices = splices.copy()
    splices["windowLeft"]  = splices["splice_center"] - window
    splices["windowRight"] = splices["splice_center"] + window

    # Cross join then filter to positions within any splice window
    matched = (
        start_pos.assign(_key=1)
        .merge(splices.assign(_key=1), on="_key")
        .drop(columns="_key")
        .query("windowLeft < pos < windowRight")
        .drop(columns=["windowLeft", "windowRight"])
        .copy()
    )

    # Left join back to preserve all positions (unmatched get NaN)
    labeled = start_pos.merge(
        matched[["pos", "splice_name", "splice_center"]],
        on="pos",
        how="left",
    )
    return labeled


def main():
    parser = argparse.ArgumentParser(
        description="Label BAM read start positions against ORF splice sites from a BED file."
    )
    parser.add_argument("bam",
        help="Path to input BAM file")
    parser.add_argument("--bed", "-b", required=True,
        help="BED file of splice sites (col 1 = 0-based start used as site center, col 3 = name)")
    parser.add_argument("--window", "-w", type=int, default=5,
        help="Window size in bp around each splice site center (default: 5, i.e. ±5 bp)")
    parser.add_argument("--output", "-o", default=None,
        help="Output CSV path (default: <bam_basename>.splice_labeled.csv)")
    args = parser.parse_args()

    if args.output is None:
        base = os.path.basename(args.bam).split(".")[0]
        args.output = f"{base}.splice_labeled.csv"

    splices = load_bed(args.bed)
    print(f"Loaded {len(splices)} splice sites from {args.bed}", file=sys.stderr)

    print(f"Reading {args.bam}...", file=sys.stderr)
    start_pos = bam_to_start_counts(args.bam)
    print(f"  {len(start_pos)} positions with positive-strand reads", file=sys.stderr)

    labeled = label_splice_sites(start_pos, splices, args.window)

    n_matched = labeled["splice_name"].notna().sum()
    print(f"  {n_matched} positions matched a splice site (window ±{args.window} bp)",
          file=sys.stderr)

    labeled.to_csv(args.output, index=False)
    print(f"Written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
