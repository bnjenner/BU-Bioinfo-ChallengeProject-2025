import os
import sys
import random
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from sgRNAtor.search import sgRNAsearch, _myers_search, _suffix_prefix_overlap

LEADER = "AACGCTGACC"  # 10 bases (matches test_bitap.py for parity)
MIN_MATCH = 6


# ── helpers ──────────────────────────────────────────────────────────────────

def make_searcher():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write(f">leader\n{LEADER}\n")
        fasta_path = f.name
    s = sgRNAsearch(fastq_files=[], leader=fasta_path)
    os.unlink(fasta_path)
    return s


_S = make_searcher()


def leader_search(read, min_match=MIN_MATCH, max_edit=0, leader=LEADER):
    return _S._sgRNAsearch__myers_leader_search(read, leader, min_match, max_edit)


# ── brute-force oracles (obviously correct; define the spec) ──────────────────

def edit_dist(a, b):
    """Global Levenshtein distance."""
    prev = list(range(len(b) + 1))
    for i in range(1, len(a) + 1):
        cur = [i]
        for j in range(1, len(b) + 1):
            cost = 0 if a[i - 1] == b[j - 1] else 1
            cur.append(min(prev[j - 1] + cost, prev[j] + 1, cur[-1] + 1))
        prev = cur
    return prev[-1]


def ref_fit(pattern, text, k):
    """Oracle for _myers_search: 5'-most end where full pattern fits (free text ends)."""
    m, n = len(pattern), len(text)
    if m == 0:
        return None
    prev = [0] * (n + 1)                     # D[0][j] = 0 (free text start)
    for i in range(1, m + 1):
        cur = [i] + [0] * n                  # D[i][0] = i
        for j in range(1, n + 1):
            cost = 0 if pattern[i - 1] == text[j - 1] else 1
            cur[j] = min(prev[j - 1] + cost, prev[j] + 1, cur[j - 1] + 1)
        prev = cur
    for j in range(1, n + 1):
        if prev[j] <= k:
            return (j - 1, prev[j])
    return None


def ref_overlap(leader, text, min_match, k):
    """Oracle for _suffix_prefix_overlap: leader suffix (>=min_match) vs text prefix."""
    m = len(leader)
    if m == 0 or min_match > m:
        return None
    n = min(len(text), m + k)
    for e in range(n):
        best_d = None
        for a in range(0, m - min_match + 1):
            d = edit_dist(leader[a:], text[:e + 1])
            if d <= k and (best_d is None or d < best_d):
                best_d = d
        if best_d is not None:
            return (e, best_d)
    return None


def ref_leader_search(read, lead, min_match, k):
    m = len(lead)
    full = ref_fit(lead, read, k)
    part = ref_overlap(lead, read, min_match, k)
    cands = []
    if full:
        cands.append((full[0], 0, full[1]))
    if part:
        cands.append((part[0], 1, part[1]))
    if not cands:
        return {"Match": False}
    cands.sort(key=lambda x: (x[0], x[1]))
    end, _, dist = cands[0]
    ol = end + 1
    return {"Match": True, "overlap_length": ol, "mismatches": dist,
            "anchor_start": m - ol, "read_position": end}


def run(name, fn):
    try:
        fn()
        print(f"  PASS  {name}")
        return True
    except AssertionError as e:
        print(f"  FAIL  {name}: {e}")
    except Exception as e:
        import traceback
        print(f"  ERROR {name}: {type(e).__name__}: {e}")
        traceback.print_exc()
    return False


# ── semantic tests (ported from test_bitap.py) ───────────────────────────────

def test_full_exact_match():
    read = LEADER + "TTTTTTTTTT"
    r = leader_search(read)
    assert r["Match"]
    assert r["overlap_length"] == len(LEADER), r
    assert r["mismatches"] == 0
    assert r["read_position"] == len(LEADER) - 1


def test_partial_overlap_at_min_match():
    read = LEADER[-MIN_MATCH:] + "TTTTTTTTTT"
    r = leader_search(read)
    assert r["Match"]
    assert r["overlap_length"] == MIN_MATCH, r
    assert r["anchor_start"] == len(LEADER) - MIN_MATCH, r


def test_wrong_partial_overlap():
    read = "TTTTTTTTTT" + LEADER[-MIN_MATCH:]
    r = leader_search(read)
    assert not r["Match"], r


def test_wrong_internal_overlap():
    read = "TTTTT" + LEADER[1:1 + MIN_MATCH] + "TTTTT"
    r = leader_search(read)
    assert not r["Match"], r


def test_partial_overlap_below_min_match():
    short = MIN_MATCH - 1
    read = LEADER[-short:] + "TTTTTTTTTT"
    r = leader_search(read)
    assert not r["Match"], r


def test_no_match():
    r = leader_search("TTTTTTTTTTTTTTTTTTTT")
    assert not r["Match"], r


def test_one_mismatch_within_budget():
    overlap = list(LEADER[-MIN_MATCH:])
    overlap[2] = "T" if overlap[2] != "T" else "A"
    read = "".join(overlap) + "GGGGGGGGGG"
    r = leader_search(read, max_edit=1)
    assert r["Match"] and r["mismatches"] == 1, r


def test_one_mismatch_exceeds_budget():
    overlap = list(LEADER[-MIN_MATCH:])
    overlap[2] = "T" if overlap[2] != "T" else "A"
    read = "".join(overlap) + "GGGGGGGGGG"
    r = leader_search(read, max_edit=0)
    assert not r["Match"], r


# ── new tests: indels (the whole point of Myers) ─────────────────────────────

def test_full_leader_one_deletion():
    # drop one base inside the leader -> needs indel-aware match
    read = LEADER[:4] + LEADER[5:] + "TTTTTTTTTT"   # deletion at pos 4
    assert not leader_search(read, max_edit=0)["Match"]
    r = leader_search(read, max_edit=1)
    assert r["Match"] and r["mismatches"] == 1, r


def test_full_leader_one_insertion():
    read = LEADER[:4] + "G" + LEADER[4:] + "TTTTTTTTTT"  # insertion at pos 4
    assert not leader_search(read, max_edit=0)["Match"]
    r = leader_search(read, max_edit=1)
    assert r["Match"] and r["mismatches"] == 1, r


def test_partial_suffix_with_deletion():
    suf = LEADER[-MIN_MATCH:]
    read = suf[:2] + suf[3:] + "GGGGGGGGGG"          # deletion inside the suffix
    r = leader_search(read, min_match=MIN_MATCH - 1, max_edit=1)
    assert r["Match"] and r["mismatches"] == 1, r


def test_genomic_prefix_full_leader_body():
    # real-data shape: [genomic 5' context][full leader][body]
    read = "ACCAACCAAC" + LEADER + "GATTACAGATTACA"
    # exact match: full leader found after the genomic prefix, trims through its end
    r = leader_search(read, max_edit=0)
    assert r["Match"] and r["mismatches"] == 0, r
    assert r["read_position"] == len("ACCAACCAAC") + len(LEADER) - 1, r
    assert read[r["overlap_length"]:].startswith("GATTACA"), r  # body is left intact


def test_two_indels_budget_boundary():
    read = LEADER[:2] + LEADER[3:6] + "A" + LEADER[6:] + "TTTTTTTTTT"  # 1 del + 1 ins
    assert not leader_search(read, max_edit=1)["Match"]
    assert leader_search(read, max_edit=2)["Match"]


# ── unit fuzz: each algorithm vs its oracle ──────────────────────────────────

def _rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _mutate(rng, s, n_edits):
    s = list(s)
    for _ in range(n_edits):
        if not s:
            break
        op = rng.choice(("sub", "ins", "del"))
        i = rng.randrange(len(s))
        if op == "sub":
            s[i] = rng.choice("ACGT")
        elif op == "ins":
            s.insert(i, rng.choice("ACGT"))
        else:
            del s[i]
    return "".join(s)


def test_fuzz_myers_search_vs_oracle():
    rng = random.Random(1234)
    fails = 0
    for _ in range(4000):
        m = rng.randint(4, 14)
        pat = _rand_seq(rng, m)
        k = rng.randint(0, 3)
        # build a text that sometimes embeds a mutated copy of the pattern
        bg = _rand_seq(rng, rng.randint(0, 12))
        if rng.random() < 0.7:
            ins = _mutate(rng, pat, rng.randint(0, 4))
            text = bg + ins + _rand_seq(rng, rng.randint(0, 12))
        else:
            text = _rand_seq(rng, rng.randint(0, 30))
        got = _myers_search(pat, text, k)
        exp = ref_fit(pat, text, k)
        if got != exp:
            fails += 1
            if fails <= 5:
                print(f"    MISMATCH pat={pat} text={text} k={k} got={got} exp={exp}")
    assert fails == 0, f"{fails} myers_search fuzz mismatches"


def test_fuzz_suffix_prefix_overlap_vs_oracle():
    rng = random.Random(555)
    fails = 0
    for _ in range(4000):
        m = rng.randint(4, 14)
        lead = _rand_seq(rng, m)
        mm = rng.randint(2, m)
        k = rng.randint(0, 3)
        if rng.random() < 0.7:
            a = rng.randint(0, m - mm) if m - mm >= 0 else 0
            suf = _mutate(rng, lead[a:], rng.randint(0, 3))
            text = suf + _rand_seq(rng, rng.randint(0, 12))
        else:
            text = _rand_seq(rng, rng.randint(0, 20))
        got = _suffix_prefix_overlap(lead, text, mm, k)
        exp = ref_overlap(lead, text, mm, k)
        if got != exp:
            fails += 1
            if fails <= 5:
                print(f"    MISMATCH lead={lead} text={text} mm={mm} k={k} got={got} exp={exp}")
    assert fails == 0, f"{fails} suffix_prefix_overlap fuzz mismatches"


def test_fuzz_leader_search_vs_oracle():
    rng = random.Random(99)
    fails = 0
    for _ in range(4000):
        m = rng.randint(6, 14)
        lead = _rand_seq(rng, m)
        mm = rng.randint(3, m)
        k = rng.randint(0, 3)
        if rng.random() < 0.6:
            a = rng.randint(0, m - mm)
            frag = _mutate(rng, lead[a:], rng.randint(0, 3))
            if a == 0 and rng.random() < 0.5:
                frag = _rand_seq(rng, rng.randint(0, 8)) + frag   # genomic 5' context
            read = frag + _rand_seq(rng, rng.randint(0, 15))
        else:
            read = _rand_seq(rng, rng.randint(0, 30))
        got = leader_search(read, mm, k, leader=lead)
        exp = ref_leader_search(read, lead, mm, k)
        # compare on the fields that matter
        gm, em = got.get("Match"), exp.get("Match")
        ok = (gm == em) and (not gm or (got["overlap_length"] == exp["overlap_length"]
                                        and got["mismatches"] == exp["mismatches"]))
        if not ok:
            fails += 1
            if fails <= 8:
                print(f"    MISMATCH lead={lead} read={read} mm={mm} k={k}\n        got={got}\n        exp={exp}")
    assert fails == 0, f"{fails} leader_search fuzz mismatches"


if __name__ == "__main__":
    tests = [
        ("full exact match",              test_full_exact_match),
        ("partial overlap == min_match",  test_partial_overlap_at_min_match),
        ("wrong partial overlap",         test_wrong_partial_overlap),
        ("wrong internal overlap",        test_wrong_internal_overlap),
        ("partial overlap < min_match",   test_partial_overlap_below_min_match),
        ("no match",                      test_no_match),
        ("1 mismatch within budget",      test_one_mismatch_within_budget),
        ("1 mismatch exceeds budget",     test_one_mismatch_exceeds_budget),
        ("full leader + 1 deletion",      test_full_leader_one_deletion),
        ("full leader + 1 insertion",     test_full_leader_one_insertion),
        ("partial suffix + deletion",     test_partial_suffix_with_deletion),
        ("genomic prefix + leader + body", test_genomic_prefix_full_leader_body),
        ("two indels budget boundary",    test_two_indels_budget_boundary),
        ("FUZZ myers_search vs oracle",   test_fuzz_myers_search_vs_oracle),
        ("FUZZ suffix_overlap vs oracle", test_fuzz_suffix_prefix_overlap_vs_oracle),
        ("FUZZ leader_search vs oracle",  test_fuzz_leader_search_vs_oracle),
    ]
    print(f"\nMyers leader search — {len(tests)} tests\n")
    passed = sum(run(name, fn) for name, fn in tests)
    print(f"\n{passed}/{len(tests)} passed\n")
    sys.exit(0 if passed == len(tests) else 1)
