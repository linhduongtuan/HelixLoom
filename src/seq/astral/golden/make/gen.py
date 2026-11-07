"""Makes test data for the species-tree engine: random species trees, gene trees from the multispecies coalescent,
with polytomies and missing taxa on request. usage: gen.py OUTDIR SEED NTAX NGENES POLY MISSING [SCALE]"""
import random, sys, os, math

def random_species_tree(names, rng, scale):
    # a random rooted binary tree by joining random pairs; branch lengths in coalescent units
    nodes = [(n, 0.0) for n in names]  # (subtree, age)
    items = [(n, 0.0) for n in names]
    age = 0.0
    while len(items) > 1:
        i, j = rng.sample(range(len(items)), 2)
        a, b = items[i], items[j]
        age = max(a[1], b[1]) + rng.uniform(0.1, 1.0) * scale
        node = ("N", a[0], b[0], age)
        items = [x for k, x in enumerate(items) if k not in (i, j)] + [(node, age)]
    return items[0][0]

def newick_species(t, parent_age=None):
    if isinstance(t, str):
        return t
    _, a, b, age = t
    def part(c):
        if isinstance(c, str):
            return c
        return newick_species(c) + ":" + format(c[3], ".4f")
    def leaf_or(c):
        if isinstance(c, str):
            return c + ":" + format(age, ".4f")
        return newick_species(c) + ":" + format(age - c[3], ".4f")
    return "(" + leaf_or(a) + "," + leaf_or(b) + ")"

def species_plain(t):
    if isinstance(t, str):
        return t
    return "(" + species_plain(t[1]) + "," + species_plain(t[2]) + ")"

def simulate_gene(t, rng):
    # returns (list of lineages (gene subtree, time at node bottom)), run coalescence up the tree
    def run(node, parent_age):
        if isinstance(node, str):
            lineages = [(node, 0.0)]
            start = 0.0
        else:
            _, a, b, age = node
            lineages = run(a, age) + run(b, age)
            start = age
        t = start
        while len(lineages) > 1:
            k = len(lineages)
            wait = rng.expovariate(k * (k - 1) / 2.0)
            if parent_age is not None and t + wait > parent_age:
                break
            t += wait
            i, j = rng.sample(range(k), 2)
            a, b = lineages[i], lineages[j]
            lineages = [x for m, x in enumerate(lineages) if m not in (i, j)] + [("(%s:%.5f,%s:%.5f)" % (a[0], t - a[1], b[0], t - b[1]), t)]
        return lineages
    res = run(t, None)
    return res[0][0] + ";"

def collapse_and_drop(newick, rng, poly, missing, names):
    # parse tiny newick into nested lists to collapse edges and drop leaves
    pos = 0
    def parse():
        nonlocal pos
        if newick[pos] == "(":
            pos += 1
            kids = []
            while True:
                kids.append(parse())
                if newick[pos] == ",":
                    pos += 1
                    continue
                break
            assert newick[pos] == ")"
            pos += 1
            length = 0.0
            if newick[pos] == ":":
                pos += 1
                s = pos
                while newick[pos] not in ",);":
                    pos += 1
                length = float(newick[s:pos])
            return ["N", kids, length]
        s = pos
        while newick[pos] not in ":,);":
            pos += 1
        name = newick[s:pos]
        length = 0.0
        if newick[pos] == ":":
            pos += 1
            s = pos
            while newick[pos] not in ",);":
                pos += 1
            length = float(newick[s:pos])
        return ["L", name, length]
    tree = parse()
    drop = set(n for n in names if rng.random() < missing)
    if len(drop) > len(names) - 4:
        drop = set(list(drop)[: max(0, len(names) - 4)])
    def prune(n):
        if n[0] == "L":
            return None if n[1] in drop else n
        kids = [k for k in (prune(c) for c in n[1]) if k is not None]
        if not kids:
            return None
        if len(kids) == 1:
            kid = kids[0]
            kid[2] += n[2]
            return kid
        return ["N", kids, n[2]]
    tree = prune(tree)
    def collapse(n, is_root=False):
        if n[0] == "L":
            return n
        kids = [collapse(c) for c in n[1]]
        out = []
        for k in kids:
            if k[0] == "N" and rng.random() < poly:
                out.extend(k[1])
            else:
                out.append(k)
        return ["N", out, n[2]]
    tree = collapse(tree, True)
    def write(n):
        if n[0] == "L":
            return n[1] + ":%.5f" % n[2]
        return "(" + ",".join(write(c) for c in n[1]) + "):%.5f" % n[2]
    s = write(tree)
    s = s[: s.rfind(":")]
    return s + ";"

if __name__ == "__main__":
    out, seed, ntax, ngenes, poly, missing = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])
    scale = float(sys.argv[7]) if len(sys.argv) > 7 else 1.0
    rng = random.Random(seed)
    names = ["t%d" % i for i in range(ntax)]
    sp = random_species_tree(names, rng, scale)
    os.makedirs(out, exist_ok=True)
    open(os.path.join(out, "true.tre"), "w").write(species_plain(sp) + ";\n")
    lines = []
    for g in range(ngenes):
        raw = simulate_gene(sp, rng)
        lines.append(collapse_and_drop(raw, rng, poly, missing, names))
    open(os.path.join(out, "genes.tre"), "w").write("\n".join(lines) + "\n")
