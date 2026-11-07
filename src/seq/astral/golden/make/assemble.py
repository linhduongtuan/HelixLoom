"""Puts the data sets, the species trees that were scored and what the real ASTRAL printed into one text file each.
usage: assemble.py DATADIR REFDIR OUTDIR NAME..."""
import os, re, sys

data, ref, out = sys.argv[1], sys.argv[2], sys.argv[3]
names = sys.argv[4:]
os.makedirs(out, exist_ok=True)

def strip_lengths(line):
    return re.sub(r":[0-9.eE+-]+", "", line.strip())

for name in names:
    d, r = os.path.join(data, name), os.path.join(ref, name)
    lines = []
    lines.append("GENES")
    for g in open(os.path.join(d, "genes.tre")).read().split("\n"):
        if g.strip():
            lines.append(strip_lengths(g))
    labels = [("true", "true.tre", "q_true.tre")]
    k = 0
    while os.path.exists(os.path.join(r, "sp%d.tre" % (k + 1))):
        k += 1
        labels.append(("sp%d" % k, "sp%d.tre" % k, "q_sp%d.tre" % k))
    for label, sp, q in labels:
        src = os.path.join(d, sp) if label == "true" else os.path.join(r, sp)
        lines.append("SPECIES " + label)
        lines.append(open(src).read().strip())
        lines.append("ASTRAL " + label)
        # the annotated tree is the last line of the file ASTRAL wrote
        lines.append(open(os.path.join(r, q)).read().strip().split("\n")[-1])
        log = open(os.path.join(r, q.replace(".tre", ".log"))).read()
        score = re.search(r"Final quartet score is: (\d+)", log).group(1)
        norm = re.search(r"Final normalized quartet score is: ([0-9.eE+-]+)", log).group(1)
        lines.append("SCORE %s %s %s" % (label, score, norm))
    for kind in ("default", "p0", "x"):
        p = os.path.join(r, kind + ".log")
        if os.path.exists(p):
            log = open(p).read()
            score = re.search(r"Final quartet score is: (\d+)", log).group(1)
            norm = re.search(r"Final normalized quartet score is: ([0-9.eE+-]+)", log).group(1)
            lines.append("BEST %s %s %s" % (kind, score, norm))
            lines.append("TREE %s" % kind)
            lines.append(open(os.path.join(r, kind + ".tre")).read().strip().split("\n")[-1])
    open(os.path.join(out, name + ".txt"), "w").write("\n".join(lines) + "\n")
    print(name, sum(len(x) for x in lines))
