#!/usr/bin/env python

# to figure out the hamming distance (minum) for all the actual indexes in this set
# is a quality cutoff score even needed?


from itertools import combinations

# read the 24 index sequences out of indexes.txt
indexes = []                                  # empty list to collect the index sequences
with open("/projects/bgmp/hodapp/bioinfo/Bi622/Demultiplex/indexes.txt", "r") as fh:
    next(fh)                                  # skip the header line
    for line in fh:
        columns = line.strip("\n").split("\t")   # split the row on tabs into columns
        index_seq = columns[4]                # index sequence is the 5th column (index 4)
        indexes.append(index_seq)             # add it to the list

# compare every unique pair and track the smallest mismatch count (in all the permutations as the minimum)
min_dist = 9        # bigger than any possible 8bp distance, so the first real pair replaces it
for ind1, ind2 in combinations(indexes, 2):        # combinations generates each unique pair (combo 2) one at a time and keeps its place; the for just asks for the next thing as per usual each loop and unpacks it into ind1, ind2 when it is handed a pair/tuple here
    dist = 0                          # counter for mismatches, starts at 0
    for x, y in zip(ind1, ind2):            # trace through both indexes position by position
        if x != y:                    # if the two characters differ
            dist += 1           # mismatch count, the positions where index1 and index2 differ, true from above so add 1, then next position... 
    # zip(ind1, ind2) compares the two 8-base strings position by position, x != y is True (1) wherever they differ, count the differs. So for one pair it returns the mismatch count.
    if dist < min_dist:   # closer than anything seen so far?
        min_dist = dist                       # remember it

print(f"number of indexes read: {len(indexes)}")
print(f"minimum Hamming distance among the known indexes: {min_dist}")