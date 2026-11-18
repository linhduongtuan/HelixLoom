"""A plain, separate implementation of the published conservation scores, used only to produce the golden numbers of the
Rust tests in `tests.rs` (run `python3 -I reference.py`; it prints the columns that the tests hold, and nothing in the
program runs it). It is checked against hand-worked cases at the top of its output, and what it is NOT: the authors' own
program, which was not run. Written from the papers (not from any program's code):

  * Henikoff & Henikoff 1994, "Position-based sequence weights", J. Mol. Biol. 243:574.
  * Capra & Singh 2007, "Predicting functionally important residues from sequence conservation", Bioinformatics 23:1875:
    the Jensen-Shannon divergence of a column from the background, times (1 - gap fraction), then a window of 3 columns on
    each side mixed in with lambda = 0.5.
  * Shannon entropy as the usual normalised conservation, 1 - H / log2(K).
Standard library only.
"""
import math
import sys

PROTEIN = "ARNDCQEGHILKMFPSTWYV"
# the background frequencies of the BLOSUM62 matrix, in the order above (BLOSUM's own: ... T W Y V), as Capra & Singh use them
# (normalised below). NOTE: the alphabet and this table must be in one order: the first version of this file had V W Y here
# and the numbers in the other order, and the golden numbers it made pinned that mistake (found in review).
BLOSUM62_BG = [0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062,
               0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072]
DNA = "ACGT"
TINY = 1e-7


def background(alphabet):
    if alphabet == PROTEIN:
        s = sum(BLOSUM62_BG)
        return [b / s for b in BLOSUM62_BG]
    return [1.0 / len(alphabet)] * len(alphabet)


def symbol(c, alphabet):
    c = c.upper()
    if alphabet == DNA and c == "U":
        c = "T"
    return alphabet.index(c) if c in alphabet else None


def weights(rows, alphabet):
    n, width = len(rows), len(rows[0])
    w = [0.0] * n
    for k in range(width):
        syms = [symbol(r[k], alphabet) for r in rows]
        present = [s for s in syms if s is not None]
        distinct = len(set(present))
        for i, s in enumerate(syms):
            if s is not None:
                w[i] += 1.0 / (distinct * present.count(s))
    total = sum(w)
    if total == 0:
        return [1.0 / n] * n
    return [x / total for x in w]


def column_score(col, w, alphabet, method, gap_penalty, max_gaps):
    syms = [symbol(c, alphabet) for c in col]
    gap = sum(wi for wi, s in zip(w, syms) if s is None)
    # the cutoff is the share of the rows with a gap; only the penalty is by weight
    if sum(1 for s in syms if s is None) / len(syms) > max_gaps + 1e-12:
        return None
    k = len(alphabet)
    counts = [0.0] * k
    for wi, s in zip(w, syms):
        if s is not None:
            counts[s] += wi
    total = sum(counts)
    if total == 0:
        return None
    p = [(c + TINY) / (total + k * TINY) for c in counts]
    if method == "entropy":
        h = -sum(x * math.log2(x) for x in p if x > 0)
        s = 1.0 - h / math.log2(k)
    else:
        q = background(alphabet)
        r = [(a + b) / 2 for a, b in zip(p, q)]
        s = 0.5 * sum(a * math.log2(a / m) for a, m in zip(p, r) if a > 0) \
            + 0.5 * sum(b * math.log2(b / m) for b, m in zip(q, r) if b > 0)
    if gap_penalty:
        s *= 1.0 - gap
    return s


def smooth(scores, window, lam):
    out = []
    for i, s in enumerate(scores):
        if s is None:
            out.append(None)
            continue
        near = [scores[j] for j in range(max(0, i - window), min(len(scores), i + window + 1))
                if j != i and scores[j] is not None]
        out.append(s if not near else (1 - lam) * s + lam * sum(near) / len(near))
    return out


def scores(rows, protein=True, method="jsd", window=3, lam=0.5, gap_penalty=True, use_weights=True, max_gaps=0.3):
    alphabet = PROTEIN if protein else DNA
    n = len(rows)
    w = weights(rows, alphabet) if use_weights else [1.0 / n] * n
    raw = [column_score([r[k] for r in rows], w, alphabet, method, gap_penalty, max_gaps) for k in range(len(rows[0]))]
    return smooth(raw, window, lam) if window > 0 else raw


def show(label, values):
    print(label, [None if v is None else round(v, 6) for v in values])


if __name__ == "__main__":
    # --- hand-checkable cases
    print("weights AD/AD/AE (hand: 7/24 7/24 10/24):", [round(x, 6) for x in weights(["AD", "AD", "AE"], PROTEIN)])
    show("DNA entropy AAAA/ACGT/AACC (hand: 1, 0, 0.5):",
         [column_score(c, [0.25] * 4, DNA, "entropy", True, 1.0) for c in ("AAAA", "ACGT", "AACC")])
    show("smooth [1,0,1] w=1 lam=.5 (hand: .5 .5 .5):", smooth([1.0, 0.0, 1.0], 1, 0.5))

    # --- golden alignments for the Rust tests
    protein = [
        "MKTAYIAKQRQISFVKSHFSRQ",
        "MKTAYIAKQRQISFVKSHFSRQ",
        "MKSAYIAKQRQVSFVKSHFSRQ",
        "MKTAYLAKQRQISFVRSHFARQ",
        "MRTAYIAKQK-ISFVKSHFSRQ",
        "MKTGYIAKQRQISFVKTHF-RQ",
        "MKTAYIVKQRQISFAKSHFSRE",
        "MQTAYIAKQRQ-SFVKSHFSRQ",
    ]
    dna = ["ACGTACGTAC", "ACGTACGAAC", "ACGAACGTAC", "ACGTTCGTAC", "ACG-ACGTAC", "ACGTACGTGC"]
    for name, rows, prot in (("protein", protein, True), ("dna", dna, False)):
        for method in ("jsd", "entropy"):
            show(f"{name} {method} raw (no window):", scores(rows, prot, method, window=0))
            show(f"{name} {method} window=3 lam=.5:", scores(rows, prot, method))
    show("protein jsd, no weights, no gap penalty, no window, max_gaps=1:",
         scores(protein, True, "jsd", window=0, gap_penalty=False, use_weights=False, max_gaps=1.0))
    print("weights protein:", [round(x, 6) for x in weights(protein, PROTEIN)])
