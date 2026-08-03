#!/usr/bin/env python

# Replotting script for Part 1.
# Reads the _means.tsv (Position, Mean_Quality_Score) that qual_dist.py already
# wrote, and makes a nicer plot from it. not redoing the whole fastq

import argparse
import matplotlib.pyplot as plt

def get_args():
    parser = argparse.ArgumentParser(description="Plot per-position mean quality from a means.tsv")
    parser.add_argument("-f", "--file", help="input _means.tsv filename", required=True)
    parser.add_argument("-o", "--out", help="output NAME only, e.g. R1", required=True)
    return parser.parse_args()

args = get_args()
f = args.file        # the _means.tsv to read
o = args.out         # clean name, e.g. R1

outdir = "part1_output"

# read the two columns back out of the tsv (skip the header line)
positions = []
means = []
with open(f, "r") as tsv:
    next(tsv)                        # skip header row
    for line in tsv:
        pos, mean = line.strip("\n").split("\t")
        positions.append(int(pos))
        means.append(float(mean))

# plot:
plt.bar(positions, means)
plt.xlabel("Base Position")
plt.ylabel("Mean Quality Score")
plt.title(f"Mean Quality Score per Position for {o}")
plt.savefig(f"{outdir}/{o}_qual.png")
 
print(f"wrote {outdir}/{o}_qual.png")