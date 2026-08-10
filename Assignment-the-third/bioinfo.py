#!/usr/bin/env python

# Author: <Kenlyn Hodapp> <optional@email.address>

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.5"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = {"A","T","G","C","N","a","t","g","c","n"}
RNA_bases = {"A","U","G","C","N","a","u","g","c","n"}

DNAcomp_dict = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"} # make a replacement dictionary to list the changes of bases, A for T, 
RNAcomp_dict = {"A": "U", "U": "A", "C": "G", "G": "C", "N": "N"}

def reverse_complement(seq: str) -> str:
    '''take a sequence of DNA and return the reverse complement of the strand'''
    complemented_bases = []              # empty list to collect each complemented base (faster than "" bc strings are immutable so new strings being constantly created. list mutable)
    for base in reversed(seq):      # walk the sequence backwards
        complement = DNAcomp_dict[base]  # look up this base's partner
        complemented_bases.append(complement)  # add it to the list
    revcomp = "".join(complemented_bases)  # glue the list into one string at the end
    return revcomp



def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score for phred+33 encoding.'''
    if not 33 <= ord(letter) <= 75:
        raise ValueError(f"'{letter}' is not a valid phred+33 character")
    return ord(letter) - 33

def qual_score(phred_score: str) -> float:
    """calculates the average phred quality score for a phred score line/string"""
    total = 0
    for letter in phred_score:
        score = convert_phred(letter)
        total += score
    average = total / len(phred_score)
    return average

def validate_base_seq(seq: str, RNAflag: bool=False) -> bool:
    '''This function takes a string. Returns True if string is composed of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    valid_bases = RNA_bases if RNAflag else DNA_bases
    for base in seq:
        if base not in valid_bases:
            return False
    return True

def gc_content(seq: str) -> float:
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    if not (validate_base_seq(seq) or validate_base_seq(seq, True)):
        raise AssertionError("String is not DNA or RNA")
    seq = seq.upper()
    Gs = seq.count("G")
    Cs = seq.count("C")
    return (Gs+Cs)/len(seq)

def calc_median(lst: list) -> float:
    '''Given a list (lst) of numbers, returns the median value of the list'''
    # median = whatever number is sitting in the middle ONCE sorted
    # list is already sorted coming in, so don't need to sort here, BUT if it wasn't sorted it would return the wrong value
    
    # check 1: checking EVERY item in lst is a number (int or float)
    # isinstance(x, (int, float)) asks "is this ONE value a number?" for each x
    # all(...) wraps that and asks "is that true for EVERY item in the list?" but if NOT all are numbers then raise
    # if even ONE item fails, all() = False, so "not False" = True, condition fires
    if not all(isinstance(x, (int, float)) for x in lst):
        # raise STOPS the function immediately, right here, nothing below runs which is different from printing the error which just displays a message and the function keeps going anyway, and possibly waste time or lead to values used later which are bad
        raise TypeError("calc_median requires numeric values")

    ordered = sorted(lst)

    n = len(ordered)        # how many numbers total
    mid = n // 2        # the position/index of the middle spot. // = floor division, drops decimals as I need an actual spot or index
    
    if n % 2 == 1:       # odd count means there IS one true middle number, no averaging needed, just return the mid
                            # ex: [1,2,100] so n=3, mid=1, lst[1] = 2 , thats just the middle one, done
        return ordered[mid]
    else:             # if its not odd, its even...even count means no single middle exists, so average the two closest to center
                            # ex: [1,2] so n=2, mid=1, so average lst[0] and lst[1] which would be (1+2)/2 = 1.5
        return (ordered[mid - 1] + ordered[mid]) / 2

def oneline_fasta(input_fasta, output_fasta):
    '''takes a fasta with sequences split over multiple lines and writes a new fasta where each entry is one header line followed by one single sequence line'''

    with open(input_fasta) as fh, open(output_fasta, "w") as outfh:
        sequence = ""                        # collects the sequence lines for the current protein

        for line in fh:
            line = line.strip()

            if line.startswith(">"):             # hit a new header
                if sequence != "":               # if were building a protein (previous), write it out first
                    outfh.write(sequence + "\n")
                    sequence = ""                # reset for the next protein
                outfh.write(line + "\n")         # write this header line
            else:
                sequence = sequence + line       # glue sequence lines together, and hold them in sequence... until next header line/loop write them out there and reset

        outfh.write(sequence + "\n")             # write the very last protein's sequence

if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    assert convert_phred("!") == 0, "wrong phred score for '!'"
    assert convert_phred("J") == 41, "wrong phred score for 'J'"
    print("Your convert_phred function is working! Nice job")
    
    assert qual_score("K") == 42.0, "wrong average phred score for 'K'"
    assert qual_score("KKKK") == 42.0, "wrong average phred score for 'KKKK'"
    assert qual_score("FFFKKK") == 39.5, "wrong average phred score for 'FFFKKK'"
    assert qual_score("$!@JK") == 23.4, "wrong average phred score for '@$!@JK'"
    assert qual_score("!") == 0.0, "wrong average phred score for '!'"
    print("qual_score is working")

    assert validate_base_seq("GGCCTTAA") == True, "fails on valid DNA"
    assert validate_base_seq("ggccttaa") == True, "fails on lowercase"
    assert validate_base_seq("GGCCUUAA", True) == True, "fails on valid RNA"
    assert validate_base_seq("GGCCZTAA") == False, "does not catch invalid character"
    print("validate_base_seq is working")
 
    assert gc_content("AATATATATATATTTATA") == 0, "wrong GC content"
    assert gc_content("GGGGGGCCCCCCCGCGCGC") == 1, "wrong GC content"
    assert gc_content("CCGT") == 0.75, "wrong GC content"
    assert gc_content("CGATTAGCCGTA") == 0.5, "wrong GC content"
    print("GC_content is working")

    assert calc_median([10,20,30,40,50]) == 30, "wrong median for odd list"
    assert calc_median([2,4,6,8]) == 5.0, "wrong median for even list"
    assert calc_median([42]) == 42, "wrong median for single item"
    assert calc_median([0,0,0,0,99]) == 0, "wrong median for skewed list"
    print("calc_median is working")