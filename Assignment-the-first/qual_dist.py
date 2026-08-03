#!/usr/bin/env python

# Assignment the First, Part 1
# For ONE fastq file, find the mean quality score at each base position across
# all reads. Uses the running-sum accumulator from PS4 (one zero per position,
# add each read's scores in, divide by number of records at the end). This is
# memory-safe: never hold all the reads at once, just the running totals.
# Writes the per-position means to a .tsv, then plots them.
 
import argparse
import gzip                       # needed to read the .gz files directly in python 
import matplotlib.pyplot as plt
import bioinfo                    # for bioinfo.convert_phred()
 
 
def get_args():
    parser = argparse.ArgumentParser(description="Mean quality score per base position for one fastq file")
    parser.add_argument("-f", "--file", help="input fastq.gz filename", required=True)
    parser.add_argument("-o", "--out", help="output name prefix for tsv and png", required=True)
    parser.add_argument("-l", "--length", help="read length (101 for bio reads, 8 for indexes)", type=int, required=True)
    return parser.parse_args()
 
args = get_args()
f = args.file        # input fastq.gz
o = args.out         # output name to prefix for tsv and png file creation
read_length = args.length    # how many positions per read
 
outdir = "part1_output"          # the script owns where files go
 
# make the accumulator: one bucket per position, all zeros to start
# [0.0] * read_length gives a row of read_length zeros (101 for reads, 8 for indexes)
my_list = [0.0] * read_length
 

 
# through the file one line at a time, summing quality scores into the buckets
num_records = 0
with gzip.open(f, "rt") as fq:        # "rt" = read text, get strings not raw bytes
    i = 0
    for line in fq:
        line = line.strip("\n")
        if i % 4 == 3:               # quality line
            for pos, letter in enumerate(line):
                my_list[pos] += bioinfo.convert_phred(letter)
            num_records += 1
        i += 1

# turn the running sums into means: divide each bucket by the number of records
for pos in range(len(my_list)):
    my_list[pos] = my_list[pos] / num_records
 
 
# write the means out to a tsv (small file, so we can re-plot later without rerunning)
with open(f"{outdir}/{o}_means.tsv", "w") as out:
    out.write("Position\tMean_Quality_Score\n")
    for pos in range(len(my_list)):
        out.write(f"{pos}\t{my_list[pos]}\n")
 
 
# plot the means
x = range(len(my_list))
y = my_list
plt.bar(x, y)
plt.xlabel("Base Position")
plt.ylabel("Mean Quality Score")
plt.title(f"Mean Quality Score per Position for {o}")
plt.savefig(f"{outdir}/{o}_qual.png")
 
print(f"{num_records} records processed")
print(f"wrote {o}_means.tsv and {o}_qual.png")