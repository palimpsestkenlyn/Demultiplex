## Goal: To sort individual reads into unique FASTQ files based on the index matched (and index status unknown or hopped) and calculate frequency of index hopping

import bioinfo   #to use qual_score (which uses convert_phred) to turn each quality character into a score. Needed to compare to qual_cutoff later

def reverse_complement(seq: str) -> str:
    '''take a sequence of DNA or RNA and return the reverse complement of the strand'''
    returns the reverse complement of a string
Input: NCTTCGAC
Expected output: GTCGAAGN

Set a quality score cutoff for the index reads = qual_cutoff # based on Part 1 per base quality distribution

Loop through the actual known 24 index file to create a set: barcodes_set

Create a dictonary to count all categorization configurations (i1 - i1, i1 - i2, i1-i3.....576 index combos, and an unknown which captures everything NOT a valid index from barcodes_set, so 577 total keys)
counting_dict = {} make emptry dictionary and initialize all possible combinations to zero 
    #keys: unknown, then every index combo (don't need a hopped because everything will be in an index combo...)
    #value: number of each
    Loop through barcodes_set to initialize all possible combos as keys with value of 0

index1 = sequence line from R2 record, index or barcode sequence on forward strand
index2 = reverse_complement(sequence line from R3 record), index or barcode sequence on reverse strand
pair_key = index1 + "-" + index2 #keep the format, two index sequences seperated by literal - to match the required header output for later

Loop through all 4 files (R1, R2, R3, R4) at the same time:
    if either index1(R2) or index2(R3) for the read contain N's
        write out the full record (all 4 lines) to the correct UNKNOWN file (R1 or R4 file based on input file) and append to the header line index1-index2
        and add 1 to counting_dict["unknown"]
    elif index1(R2) or index2(R3) not in barcodes_set:
         write out the full record (all 4 lines) to the correct UNKNOWN file (R1 or R4 file based on input file header) and append to the header line index1-index2
         and add 1 to counting_dict["unknown"]
    elif the quality score of index1(R2) or index2(R3) is below qual_cutoff:
        write out the full record (all 4 lines) to the correct UNKNOWN file (R1 or R4 file based on input file header) and append to the header line index1-index2
        add 1 to counting_dict["unknown"]
    else:
        if index1(R2) == index2(R3) # if index1 and index2 match:
            write out the full record (all 4 lines) to the appropriate file for that index (and R1 or R4 based on the input file) and append to the header line index1-index2
            and add 1 to counting_dict[pair_key]
        else they do not match:
            so write out the full record (all 4 lines) to the appropriate hopped file (R1 or R4 based on input file) and append to the header line index1-index2
            and add 1 to counting_dict[pair_key]

Return/Report total matched index reads (sum of all keys where the two indexs' are equal)
Return/Report total hopped (sum of all keys where the two indexs' differ)
Return/Report the unknown count 
Return/Report each of the three as a percentage of total read pairs (after counting the total of read pairs)


output files: 52 in total
    names:
        index_read1 or index_read2 BUT ACTUALLY read 1 or read 2 (so from R1 or R4 source file) = 48 files
        hopped_R1 or hopped_R2 = 2 files
        unknown_R1 or unknown_R2 = 2 files



