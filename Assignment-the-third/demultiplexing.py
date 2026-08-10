import argparse
import bioinfo   # for reverse_complement; no quality cutoff needed (min Hamming distance = 3)
import gzip
import itertools

def get_record(fh) -> list[str]:
    '''Read 4 lines of a file (one read in a FASTQ record) from an open file handle.
    Returns each line as a different entry in a list: [header, sequence, plus, quality].
    so each can be used/referenced individually'''
    header = fh.readline().strip("\n")     # line 1: header
    seq = fh.readline().strip("\n")        # line 2: sequence
    plus = fh.readline().strip("\n")       # line 3: plus sign
    qual = fh.readline().strip("\n")       # line 4: quality
    return [header, seq, plus, qual]

def get_args():
    parser = argparse.ArgumentParser(description="Demultiplex dual-indexed FASTQ reads by index-pair and report on index hopping")
    parser.add_argument("-r1", "--read1", help="input R1 FASTQ (biological read 1, forward strand)", required=True)
    parser.add_argument("-r2", "--read2", help="input R2 FASTQ (index 1, pairs forward read or R1 file)", required=True)
    parser.add_argument("-r3", "--read3", help="input R3 FASTQ (index 2, pairs with reverse read or R4 file and is reverse complement)", required=True)
    parser.add_argument("-r4", "--read4", help="input R4 FASTQ (biological read 2, reverse strand)", required=True)
    parser.add_argument("-i", "--indexes", help="known indexes file for valid barcodes/indexes used (indexes.txt)", required=True)
    parser.add_argument("-o", "--outdir", help="output directory for demultiplexed files", required=True)
    return parser.parse_args()

args = get_args()
r1 = args.read1      # biological read 1 (forward)
r2 = args.read2      # index 1
r3 = args.read3      # index 2
r4 = args.read4      # biological read 2 (reverse)
index_file = args.indexes   # the known indexes
outdir = args.outdir     # where all output files get written

# Read the 24 known indexes into a set, set is faster lookup compared to a list, no key/value so no dict
# Store them in forward orientation (as they appear in the txt file / as R2 gives them).
barcodes = set()
with open(index_file, "r") as fh:
    next(fh)                # skip the header row
    for line in fh:
        columns = line.strip("\n").split("\t")   # split the row on tabs into columns
        sequence = columns[4]         # 5th column is the index sequence (columns start at 0)
        barcodes.add(sequence)        # add it to the set

# print(len(barcodes))
# print(barcodes)

# make counting dictonary to track things as it goes
# every ordered pair of known indexes (both halves drawn from barcodes),
# including an index paired with itself = 576 combos total have to use itertools product
counting_dict = {}
for first_index, second_index in itertools.product(barcodes, repeat=2):
    index_pair = first_index + "-" + second_index
    counting_dict[index_pair] = 0
counting_dict["unknown"] = 0 #set manually as one category key not coming from index sequences.

# make dictonaries for matched files to write to per index per read (48 total)
# key = index, value = the open output file for that index
r1_FHdict = {}
r2_FHdict = {}

for index in barcodes:              #go through each index sequence (now in barcodes set from above)
    r1_FHdict[index] = open(f"{outdir}/{index}_R1.fastq", "w") #open a file for writing in the output files location, name it the <index sequence>_R1.fastq and make this the value at the specific <index sequence> key location
    r2_FHdict[index] = open(f"{outdir}/{index}_R2.fastq", "w") #at each spot/key in the dictonary, so for each <index seq> key in the dict, open a file named as outlined and set it as the value (associated with that <index seq> as the key)

# create unknown and hopped files
hopped_R1 = open(f"{outdir}/hopped_R1.fastq", "w")
hopped_R2 = open(f"{outdir}/hopped_R2.fastq", "w")
unknown_R1 = open(f"{outdir}/unknown_R1.fastq", "w")
unknown_R2 = open(f"{outdir}/unknown_R2.fastq", "w")


# print(len(r1_FHdict))
# print(len(r2_FHdict))
# print(r1_FHdict) #these look insane?
# print(r2_FHdict) #these look insane?

# WHILE USING get_record- each record becomes a list of its 4 FASTQ lines: [0]=header, [1]=sequence, [2]=plus, [3]=quality score/phred
r1_record = [...]   # current record from R1 (biological forward)
r2_record = [...]   # current record from R2 (index 1)
r3_record = [...]   # current record from R3 (index 2)
r4_record = [...]   # current record from R4 (biological reverse)

#Step 1: open the four input files + loop skeleton that reads records and stops

# open the four input files for reading (rt = read text, since they are gzipped)
r1_fh = gzip.open(r1, "rt")         #files can be referenced with arparse argument variable assignments from the beginning
r2_fh = gzip.open(r2, "rt")
r3_fh = gzip.open(r3, "rt")
r4_fh = gzip.open(r4, "rt")

# temporary counter to test that reading works, delete after Step 1
# record_count = 0

while True:
    r1_record = get_record(r1_fh)   # each record is [header, seq, plus, qual]
    r2_record = get_record(r2_fh)
    r3_record = get_record(r3_fh)
    r4_record = get_record(r4_fh)

    if r1_record[0] == "":     # header came back empty -> files are out, stop
        break

    # record_count += 1          # temporary: count each record read
    
# temporary: confirm the loop read the right number of records
# print(record_count)
    index1 = r2_record[1]                               # R2 sequence is already index/barcode 1, position 1 is sequence so grab and set as index1
    index2 = bioinfo.reverse_complement(r3_record[1])   # R3 sequence is index2, but needs reverse complimented to be correct orientation
    pair = index1 + "-" + index2                        # create header format: index1-index2 with literal - between
    #temporary
    # print(pair) # test for expected 8 from test files

    # Sort each read into one of three buckets: unknown, hopped, or matched.
    # Order matters: rule out unknown first, then the survivors are all valid, then it just becomes either they match or don't
    # indexes and split cleanly into matched (same) vs hopped (different).

    if index1 not in barcodes or index2 not in barcodes:    # unknown: at least one index is not one of the 24.
        # catches any index with an N, since no real index contains an N, an N-containing index will never be in the barcodes set. (no separate N check needed.)
        print("unknown", pair)          # temporary: check sorting
        new_r1_header = r1_record[0] + " " + pair
        new_r4_header = r4_record[0] + " " + pair
        unknown_R1.write(new_r1_header + "\n" + r1_record[1] + "\n" + r1_record[2] + "\n" + r1_record[3] + "\n")
        unknown_R2.write(new_r4_header + "\n" + r4_record[1] + "\n" + r4_record[2] + "\n" + r4_record[3] + "\n")
        counting_dict["unknown"] += 1

    elif index1 == index2:                                  # matched: both indexes valid AND identical -> real dual-matched pair
        print("matched", pair)          # temporary
        new_r1_header = r1_record[0] + " " + pair
        new_r4_header = r4_record[0] + " " + pair
        r1_FHdict[index1].write(new_r1_header + "\n" + r1_record[1] + "\n" + r1_record[2] + "\n" + r1_record[3] + "\n")
        r2_FHdict[index1].write(new_r4_header + "\n" + r4_record[1] + "\n" + r4_record[2] + "\n" + r4_record[3] + "\n")
        counting_dict[pair] += 1

    else:                                                   # hopped: both indexes valid but different from each other -> index hop
        print("hopped", pair)           # temporary
        new_r1_header = r1_record[0] + " " + pair
        new_r4_header = r4_record[0] + " " + pair
        hopped_R1.write(new_r1_header + "\n" + r1_record[1] + "\n" + r1_record[2] + "\n" + r1_record[3] + "\n")
        hopped_R2.write(new_r4_header + "\n" + r4_record[1] + "\n" + r4_record[2] + "\n" + r4_record[3] + "\n")
        counting_dict[pair] += 1




# close all output files 
for fh in r1_FHdict.values():
    fh.close()
for fh in r2_FHdict.values():
    fh.close()
hopped_R1.close()
hopped_R2.close()
unknown_R1.close()
unknown_R2.close()

# close input files
r1_fh.close()
r2_fh.close()
r3_fh.close()
r4_fh.close()

# write raw per-pair counts to a tsv for data export/manipulation after the fact
# any later formatting or figures derive from this saved data
with open(f"{outdir}/counts.tsv", "w") as tsv:
    tsv.write("index_pair\tcount\n")          # header row
    for index_pair in counting_dict:          # every key in the dict
        tsv.write(index_pair + "\t" + str(counting_dict[index_pair]) + "\n")


# tally totals from the counting dict for the summary
total_reads = 0
matched_total = 0
hopped_total = 0

for index_pair in counting_dict:
    count = counting_dict[index_pair]
    total_reads += count                     # every count adds to the grand total

    if index_pair == "unknown":
        continue                             # unknown has no two halves to compare, skip it here otherwise stalls
    first, second = index_pair.split("-")    # split the pair back into its two halves
    if first == second:
        matched_total += count               # halves equal -> matched so plus 1 count
    else:
        hopped_total += count                # halves differ -> hopped so plust 1 count

unknown_total = counting_dict["unknown"]     # grab unknown directly no need to count, already total/value in dict

# temporary: print(total_reads, matched_total, hopped_total, unknown_total)

# write the human-readable summary to a markdown file
with open(f"{outdir}/summary.md", "w") as summary:
    summary.write("# Demultiplexing Summary\n\n")
    summary.write(f"Total read-pairs: {total_reads}\n\n")
    # each total as a count and as a percent of all reads, rounded to 2 decimals
    summary.write(f"Matched: {matched_total} ({matched_total/total_reads*100:.2f}%)\n\n")
    summary.write(f"Index-hopped: {hopped_total} ({hopped_total/total_reads*100:.2f}%)\n\n")
    summary.write(f"Unknown: {unknown_total} ({unknown_total/total_reads*100:.2f}%)\n\n")

    # per-sample breakdown: only matched pairs (halves equal) that actually got reads
    summary.write("## Percentage of reads per sample (matched pairs)\n\n")
    for index_pair in counting_dict:
        if index_pair == "unknown":
            continue                        # unknown has no dash to split, skip it otherwise stalls
        first, second = index_pair.split("-")
        if first == second and counting_dict[index_pair] > 0:   # matched and non-zero
            count = counting_dict[index_pair]
            summary.write(f"{index_pair}: {count} ({count/total_reads*100:.2f}%)\n\n")