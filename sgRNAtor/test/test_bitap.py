import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from sgRNAtor.search import sgRNAsearch

LEADER = "AACGCTGACC"  # 10 bases
MIN_MATCH = 6


def make_searcher(min_match=MIN_MATCH, max_edit=0):
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write(f">leader\n{LEADER}\n")
        fasta_path = f.name
    s = sgRNAsearch(fastq_files=[], leader=fasta_path)
    os.unlink(fasta_path)
    s._sgRNAsearch__build_bitmask(min_match, max_edit)
    return s


def search(s, read, min_match=MIN_MATCH, max_edit=0):
    return s._sgRNAsearch__bitap_partial_overlap(read, LEADER, min_match, max_edit)


def run(name, fn):
    try:
        fn()
        print(f"  PASS  {name}")
    except AssertionError as e:
        print(f"  FAIL  {name}: {e}")
    except Exception as e:
        print(f"  ERROR {name}: {type(e).__name__}: {e}")


# ── tests ──────────────────────────────────────────────────────────────────

def test_full_exact_match():
    s = make_searcher()
    read = LEADER + "TTTTTTTTTT"
    r = search(s, read)
    print("\t", read, "\n\t", LEADER, "\n\t", r)
    assert r["Match"], "expected match"
    assert r["overlap_length"] == len(LEADER), f"expected {len(LEADER)}, got {r['overlap_length']}"
    assert r["mismatches"] == 0
    assert r["anchor_start"] == 0
    assert r["read_position"] == len(LEADER) - 1


def test_partial_overlap_at_min_match():
    s = make_searcher()
    # Last MIN_MATCH chars of LEADER at 5' end of read
    read = LEADER[-MIN_MATCH:] + "TTTTTTTTTT"
    r = search(s, read)
    print("\t", read, "\n\t", LEADER, "\n\t", r)
    assert r["Match"], "expected match"
    assert r["overlap_length"] == MIN_MATCH, f"expected {MIN_MATCH}, got {r['overlap_length']}"
    assert r["anchor_start"] == len(LEADER) - MIN_MATCH


def test_wrong_partial_overlap_at_min_match():
    s = make_searcher()
    read = "TTTTTTTTTT" + LEADER[-MIN_MATCH:]
    r = search(s, read)
    print("\t", read, "\n\t", LEADER, "\n\t", r)
    assert not r["Match"], "expected no match"


def test_wrong_internal_overlap_at_min_match():
    s = make_searcher()
    read = "TTTTT" + LEADER[1:1+MIN_MATCH] + "TTTTT"
    r = search(s, read)
    print("\t", read, "\n\t", LEADER, "\n\t", r)
    assert not r["Match"], "expected no match"


def test_partial_overlap_below_min_match():
    s = make_searcher()
    short = MIN_MATCH - 1
    read = LEADER[-short:] + "TTTTTTTTTT"
    r = search(s, read)
    print("\t", read, "\n\t", LEADER, "\n\t", r)
    assert not r["Match"], "expected no match"


def test_no_match():
    s = make_searcher()
    read =  "TTTTTTTTTTTTTTTTTTTT"
    r = search(s, read)
    print("\t", read, "\n\t", LEADER, "\n\t", r)
    assert not r["Match"], "expected no match"


def test_one_mismatch_within_budget():
    s = make_searcher(max_edit=1)
    overlap = list(LEADER[-MIN_MATCH:])
    overlap[2] = "T" if overlap[2] != "T" else "A"  # introduce one substitution
    read = "".join(overlap) + "GGGGGGGGGG"
    r = search(s, read, max_edit=1)
    print("\t", read, "\n\t", LEADER, "\n\t", r)
    assert r["Match"], "expected match with 1 mismatch allowed"
    assert r["mismatches"] == 1, f"expected 1 mismatch, got {r['mismatches']}"


def test_one_mismatch_exceeds_budget():
    s = make_searcher(max_edit=0)
    overlap = list(LEADER[-MIN_MATCH:])
    overlap[2] = "T" if overlap[2] != "T" else "A"
    read = "".join(overlap) + "GGGGGGGGGG"
    r = search(s, read, max_edit=0)
    print("\t", read, "\n\t", LEADER, "\n\t", r)
    assert not r["Match"], "expected no match when mismatch exceeds budget"


if __name__ == "__main__":
    tests = [
        ("full exact match",              test_full_exact_match),
        ("partial overlap == min_match",  test_partial_overlap_at_min_match),
        ("wrong partial overlap",         test_wrong_partial_overlap_at_min_match),
        ("wrong internal overlap",        test_wrong_internal_overlap_at_min_match),
        ("partial overlap < min_match",   test_partial_overlap_below_min_match),
        ("no match",                      test_no_match),
        ("1 mismatch within budget",      test_one_mismatch_within_budget),
        ("1 mismatch exceeds budget",     test_one_mismatch_exceeds_budget),
    ]
    print(f"\nBitap search — {len(tests)} tests\n")
    for name, fn in tests:
        run(name, fn)
    print()
