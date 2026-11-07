import random, sys
out, seed, ntax, count = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
rng = random.Random(seed)
lines = []
for _ in range(count):
    items = ["t%d" % i for i in range(ntax)]
    while len(items) > 1:
        a = items.pop(rng.randrange(len(items))); b = items.pop(rng.randrange(len(items)))
        items.append("(%s,%s)" % (a, b))
    lines.append(items[0] + ";")
open(out, "w").write("\n".join(lines) + "\n")
