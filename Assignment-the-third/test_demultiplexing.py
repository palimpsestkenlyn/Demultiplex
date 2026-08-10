import gzip

def get_record(fh) -> list[str]:
    '''Read 4 lines of a file (one read in a FASTQ record) from an open file handle.
    Returns each line as a different entry in a list: [header, sequence, plus, quality].
    so each can be used/referenced individually'''
    header = fh.readline().strip("\n")     # line 1: header
    seq = fh.readline().strip("\n")        # line 2: sequence
    plus = fh.readline().strip("\n")       # line 3: plus sign
    qual = fh.readline().strip("\n")       # line 4: quality
    return [header, seq, plus, qual]

fh = gzip.open("/projects/bgmp/hodapp/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test_R2.fastq.gz", "rt")
record = get_record(fh)
fh.close()
assert record[0].startswith("@"), "header should start with @"
assert record[2] == "+", "third line should be +"
assert len(record[1]) == len(record[3]), "seq and qual must be same length"
print("get_record passed")