#!/usr/bin/env python3
import argparse
from collections import defaultdict

def main():
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("output")
    p.add_argument("--max-hits", type=int, default=5)
    p.add_argument("--evalue", type=float, default=1e-5)
    a = p.parse_args()
    hits = defaultdict(list)
    with open(a.input, encoding="utf-8") as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 12 and f[0] != f[1] and float(f[10]) <= a.evalue:
                hits[f[0]].append(f)
    with open(a.output, "w", encoding="utf-8", newline="\n") as out:
        for query in sorted(hits):
            for f in sorted(hits[query], key=lambda x: (-float(x[11]), x[1]))[:a.max_hits]:
                out.write("\t".join(f) + "\n")

if __name__ == "__main__":
    main()
