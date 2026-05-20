#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path
from collections import defaultdict


def parse_ground_truth(ground_truth_file):
    """
    Parse a tab-delimited ground truth file with rows like:

        subgenome_name    number_expected    S104_2538    S104_3381    S104_2677

    Returns:
        ground_truth: dict
            {
                "ORF6": {"S104_2538", "S104_3381", ...},
                "ORF7": {...}
            }

        subgenome_order: list
            The order of sgRNAs/subgenomes as they appear in the ground truth file.
    """

    ground_truth = {}
    subgenome_order = []

    with open(ground_truth_file, "r", newline="") as f:
        reader = csv.reader(f, delimiter="\t")

        for row in reader:
            if not row:
                continue

            # Skip comment lines
            if row[0].startswith("#"):
                continue

            subgenome_name = row[0].strip()

            # Optional: skip a header row if present
            if subgenome_name.lower() in {"subgenome_name", "subgenome", "sgrna"}:
                continue

            if len(row) < 2:
                raise ValueError(f"Malformed ground truth row: {row}")

            # row[1] is number_expected
            # row[2:] are expected read IDs
            expected_reads = set(x.strip() for x in row[2:] if x.strip())

            ground_truth[subgenome_name] = expected_reads
            subgenome_order.append(subgenome_name)

    return ground_truth, subgenome_order


def extract_read_id_from_header(header_line):
    """
    Extract read ID from a fasta header such as:

        >0__S101_31 left:30 start:69 end:26237 right:26516

    Returns:
        S101_31
    """

    header_line = header_line.strip()

    if header_line.startswith(">"):
        header_line = header_line[1:]

    # Keep only first whitespace-delimited field:
    # 0__S101_31
    first_field = header_line.split()[0]

    # Split on double underscore:
    # ["0", "S101_31"]
    if "__" in first_field:
        return first_field.split("__", 1)[1]

    # Fallback if the expected pattern is missing
    return first_field


def parse_fasta_read_ids(fasta_file):
    """
    Parse a fasta file and return the set of read IDs found in the headers.
    """

    read_ids = set()

    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                read_id = extract_read_id_from_header(line)
                read_ids.add(read_id)

    return read_ids


def parse_run_results(run_dir, fasta_extension=".fasta"):
    """
    Parse all sgRNA fasta files for one run.

    Expected structure:

        run_id/
          results/
            ORF6.fasta
            ORF7.fasta

    Returns:
        results: dict
            {
                "ORF6": {"S101_31", "S101_32", ...},
                "ORF7": {...}
            }
    """

    run_dir = Path(run_dir)
    results_dir = run_dir / "results/fasta"

    if not results_dir.exists():
        raise FileNotFoundError(f"Could not find results directory: {results_dir}")

    run_results = {}

    for fasta_file in sorted(results_dir.glob(f"*{fasta_extension}")):
        subgenome_name = fasta_file.split("-")[-1].replace(fasta_extension, "")
        read_ids = parse_fasta_read_ids(fasta_file)
        run_results[subgenome_name] = read_ids

    return run_results


def write_run_output(run_id, run_results, subgenome_order, output_dir):
    """
    Write one tab-delimited result file for a run in the format:

        subgenome_name    number_identified    read_id_1    read_id_2    ...

    The second column is named number_expected to match the ground truth format,
    but it contains the number of reads identified by the tool.
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / f"{run_id}_results.tsv"

    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")

        #writer.writerow(["subgenome_name", "number_expected"])

        for subgenome_name in subgenome_order:
            read_ids = sorted(run_results.get(subgenome_name, set()))
            writer.writerow([subgenome_name, len(read_ids), *read_ids])

        # Also include sgRNAs found by the tool but not present in ground truth
        extra_subgenomes = sorted(set(run_results) - set(subgenome_order))
        for subgenome_name in extra_subgenomes:
            read_ids = sorted(run_results[subgenome_name])
            writer.writerow([subgenome_name, len(read_ids), *read_ids])

    return output_file


def calculate_confusion_matrix(ground_truth, run_results):
    """
    Calculate TP, FP, FN, TN.

    Positive means:
        A particular read belongs to a particular sgRNA/subgenome.

    For each subgenome, compare:

        expected reads from ground truth
        identified reads from tool output

    TP:
        read expected for that subgenome and identified for that subgenome

    FP:
        read not expected for that subgenome but identified for that subgenome

    FN:
        read expected for that subgenome but not identified for that subgenome

    TN:
        read not expected for that subgenome and not identified for that subgenome

    Important:
        TN requires defining the full universe of possible read IDs.
        Here, the universe is defined as all reads appearing anywhere in either
        the ground truth or the run results.
    """

    all_subgenomes = set(ground_truth) | set(run_results)

    all_reads = set()
    for reads in ground_truth.values():
        all_reads.update(reads)
    for reads in run_results.values():
        all_reads.update(reads)

    tp = 0
    fp = 0
    fn = 0
    tn = 0

    for subgenome_name in all_subgenomes:
        expected_reads = ground_truth.get(subgenome_name, set())
        identified_reads = run_results.get(subgenome_name, set())

        for read_id in all_reads:
            expected = read_id in expected_reads
            identified = read_id in identified_reads

            if expected and identified:
                tp += 1
            elif not expected and identified:
                fp += 1
            elif expected and not identified:
                fn += 1
            else:
                tn += 1

    return tp, fp, fn, tn


def find_run_dirs(parent_dir):
    """
    Find run directories under a parent directory.

    A run directory is defined as any directory that contains:

        results/

    Example:

        parent_dir/
          run_001/
            results/
          run_002/
            results/
    """

    parent_dir = Path(parent_dir)

    run_dirs = []

    for child in sorted(parent_dir.iterdir()):
        if child.is_dir() and (child / "results").is_dir():
            run_dirs.append(child)

    return run_dirs


def main():
    parser = argparse.ArgumentParser(
        description="Analyze sgRNA identification results across multiple runs."
    )

    parser.add_argument(
        "-r",
        "--runs-dir",
        required=True,
        help="Directory containing run_id folders.",
    )

    parser.add_argument(
        "-g",
        "--ground-truth-dir",
        required=True,
        help="Directory containing one ground truth file per run.",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="Directory where per-run result files and summary CSV will be written.",
    )

    parser.add_argument(
        "--ground-truth-extension",
        default=".tsv",
        help="Extension for ground truth files. Default: .tsv",
    )

    parser.add_argument(
        "--fasta-extension",
        default=".fasta",
        help="Fasta file extension. Default: .fasta",
    )

    parser.add_argument(
        "--summary-name",
        default="confusion_matrix_summary.csv",
        help="Name of the output confusion matrix CSV.",
    )

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ground_truth_dir = Path(args.ground_truth_dir)

    run_dirs = find_run_dirs(args.runs_dir)

    if not run_dirs:
        raise RuntimeError(
            f"No run directories containing results/ found in {args.runs_dir}"
        )

    summary_file = output_dir / args.summary_name

    with open(summary_file, "w", newline="") as f:
        writer = csv.writer(f)

        writer.writerow(
            [
                "run_id",
                "ground_truth_file",
                "true_positives",
                "false_positives",
                "false_negatives",
                "true_negatives",
            ]
        )

        for run_dir in run_dirs:
            run_id = run_dir.name

            print(f"Processing run: {run_id}")

            ground_truth_file = ground_truth_dir / f"{run_id[:-9]}{args.ground_truth_extension}"

            if not ground_truth_file.exists():
                raise FileNotFoundError(
                    f"Could not find ground truth file for run {run_id}: "
                    f"{ground_truth_file}"
                )

            ground_truth, subgenome_order = parse_ground_truth(ground_truth_file)

            run_results = parse_run_results(
                run_dir,
                fasta_extension=args.fasta_extension,
            )

            run_output_file = write_run_output(
                run_id=run_id,
                run_results=run_results,
                subgenome_order=subgenome_order,
                output_dir=output_dir,
            )

            tp, fp, fn, tn = calculate_confusion_matrix(
                ground_truth=ground_truth,
                run_results=run_results,
            )

            writer.writerow(
                [
                    run_id,
                    str(ground_truth_file),
                    tp,
                    fp,
                    fn,
                    tn,
                ]
            )

            print(f"  Ground truth: {ground_truth_file}")
            print(f"  Wrote: {run_output_file}")
            print(f"  TP={tp}, FP={fp}, FN={fn}, TN={tn}")

    print()
    print(f"Done. Summary written to: {summary_file}")


if __name__ == "__main__":
    main()