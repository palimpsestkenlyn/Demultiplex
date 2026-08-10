## Overview


* Start date: 7.21.26
define the problem, input and output, functions skeleton(doc string, input/output), pseudocode...ALL due sat at 10 am. AFTER 10 am sat, review other students psuedocode due TUESDAY next week

NOTHING FOR PART 2 in answers MD file....
## Questions/Process:
* ==8.4.26- for my `get_record` function, how would a type annotation work here? type annotations are like PYTHON types right? not like my logic so I can't say the input is a fastq file or can i?==
* 7.23.26: Do the barcode appends need the barcode-reverse complement seq? or just barcode-barcode in the order that it would appear in the txt file/key?
	* Ok so regular index orders NOT reverse complement because that way I can see what the index hop pairs are.
* For cutoff scores (in part 2) are we saying a cutoff score for the BARCODES only? or are we also doing some comparison of the quality score of the reads as well?  Do part 1 cutoff and part 2 cutoff have overlap or is it?
	* **Part 1 averages down the columns.** Meaning one number PER POSITION, computed across every read in the file. So for barcodes, eight numbers total, = histogram bars. 
	* **Part 2 averages across a row.** One READ at a time, meaning for barcode its own 8 scores collapsed to one number like an average? Using qual_score and then compare to qual_cutoff which I set? 
		* If so it seems like the error issue in a barcode is VERY different than in a read. in the length of a read a few wrong base calls seems somewhat diluted by coverage or matches at other positions. In an index, how impactful is one base being called wrong? IT DEPENDS ON HAMMING DISTANCE:
	* Do we NEED a qual_cutoff, if we have looked at hamming distance, and the pool of INDEXES that are being used have sufficient hamming distance. then even 1 or 2 base changes (for a hamming distance of 3 or 4) would not end up resulting in the misidentification OF the swapped index to looking like a valid index. soooo the quality cutoff could not even need to be an average, but could just be like are there less than Ns or less.... depends on the specific hamming distances. SO THINK THIS THROUGH
		* SO DO QUALTIFY CUTOFF BEFORE EVERYTHING, throw out, and then the other unknowns and only THEN deal with the data.
## Initial Data Exploration: 
* Start date: 7.23.26
* examined files in /projects/bgmp/shared/2017_sequencing/ on talapas to get a sense of what is being worked with. 
* 4 fastq files and a indexes.txt file
````
-rw-r-xr--+ 1 coonrod  is.racs.pirg.bgmp 20927576171 Jul 30  2018 1294_S1_L008_R1_001.fastq.gz
-rw-r-xr--+ 1 coonrod  is.racs.pirg.bgmp  2702499002 Jul 30  2018 1294_S1_L008_R2_001.fastq.gz
-rw-r-xr--+ 1 coonrod  is.racs.pirg.bgmp  2992503564 Jul 30  2018 1294_S1_L008_R3_001.fastq.gz
-rw-r-xr--+ 1 coonrod  is.racs.pirg.bgmp 21814991385 Jul 30  2018 1294_S1_L008_R4_001.fastq.gz
````
* R1 sequence data from forward read
* R2 barcode reads for R1 reads
* R3 barcode reads from R4 reads (reverse compliment to the R2 barcodes)
* R4 sequence data from reverse reads
* indexes.txt the 24 unique barcodes as well as other uncessary data. all barcodes 8 nt in length
* Head files to confirm data layout: (confirmed R1 format matches R4, R2 matches R3)
```
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz | head -4
@K00337:83:HJKJNBBXX:8:1101:1265:1191 1:N:0:1
GNCTGGCATTCCCAGAGACATCAGTACCCAGTTGGTTCAGACAGTTCCTCTATTGGTTGACAAGGTCTTCATTTCTAGTGATATCAACACGGTGTCTACAA
+
A#A-<FJJJ<JJJJJJJJJJJJJJJJJFJJJJFFJJFJJJAJJJJ-AJJJJJJJFFJJJJJJFFA-7<AJJJFFAJJJJJF<F--JJJJJJF-A-F7JJJJ

zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | head -4
@K00337:83:HJKJNBBXX:8:1101:1265:1191 2:N:0:1
NCTTCGAC
+
#AA<FJJJ
```
  * files have 1,452,986,940 lines ÷ 4 = 363,246,735 reads (did not do this on command line as it takes a long time, was done in class)
  * Figure out read length in each file: wc -c counts newline character at end of line, so actual length is return minus 1:
  ````
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz | head -2 | tail -1 | wc -c
102 
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | head -2 | tail -1 | wc -c
9
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | head -2 | tail -1 | wc -c
9 
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz | head -2 | tail -1 | wc -c
102
  ````
  * Can tell immediately this is phred+33 encoding because there are symbols like # and < which are ASCII symbols that fall below the range of what would be seen for phred+64. Phred+33 uses ASCII ranges of 33 to 126, and seeing # or < which correspond to 35 and 60 which are within this range confirms this is phred+33. Looking at one file, you could assume that all would match from the same run, but just in case all 4 files show # within the first few reads, which would not be seen in phred+64. 
  * Create test files for further development: grab 8 reads (want to test 2 indexes and various combos of hopped, matched, unknown) and redirect to new fastq files. Can manually make adjustments to data there to ensure all categories met for testing. Repeat for all 4 files...
  ```
  zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz | head -32 > TEST-input_FASTQ/test_R1.fastq
  ```
7.31.26: 
Create test files that hit the 3 categories in various combos, grabbed real reads/exports from files, than manually manipulated ONLY the barcodes/index sequences in the R2/R3 files for testing. 
```
Reads 1 and 2, dual matched (R2 index and R3 = reverse complement of the same index):

Read 1, matched B1:  GTAGCGTA-GTAGCGTA
R2: `GTAGCGTA` (B1)  
R3: `TACGCTAC` (revcomp of B1)

Read 2, matched A5:  CGATCGAT-CGATCGAT
R2: `CGATCGAT` (A5)  
R3: `ATCGATCG` (revcomp of A5)

Reads 3 and 4, index-hopped (both real indexes, but different from each other):

Read 3, hop B1/C1:  GTAGCGTA-GATCAAGG
R2: `GTAGCGTA` (B1)  
R3: `CCTTGATC` (revcomp of C1, which is GATCAAGG)

Read 4, hop A5/B9:  CGATCGAT-AACAGCGA
R2: `CGATCGAT` (A5)  
R3: `TCGCTGTT` (revcomp of B9, which is AACAGCGA)

Reads 5 and 6, unknown via an N in the barcode, left as is:

Read 5, N in R2:  NACCGGAT-TACCGGAN
R2: `NACCGGAT`   
R3: `NTCCGGTA` 

Read 6, N in R3:  NGTTCCGT-TGTTCCGN
R2: `NGTTCCGT`   
R3: `NCGGAACA` 

Reads 7 and 8, unknown via a barcode that isn't on the 24-index list:

Read 7, junk in R2:  AAAAAAAA-GTAGCGTA
R2: `AAAAAAAA` (not a real index)  
R3: `TACGCTAC` (revcomp of B1, real, but doesn't matter, the pair is already unknown)

Read 8, junk in R3:  AAAAAAAA-AAAAAAAA
R2: `AAAAAAAA`   
R3: `TTTTTTTT` (revcomp of AAAAAAAA so will match, but not in index dict)
```
8.1.26
## 8.1.26 - Part 1: Quality Score Distribution

Wrote `qual_dist.py`. Outputs a `_means.tsv` of the
per-position means and a `_qual.png` plot. The tsv is saved so plots can be
regenerated without rerunning the full job.

Wrote `plot_qual.py` to replot from a saved `_means.tsv` for adjusting plot appearance.

Tested `qual_dist.py` on the small test files (gzipped test_R1) to confirm it
runs and produces the tsv and png outputs correctly.

Ran via sbatch. Confirmed working on R1 (bio read, -l 101), then submitted
R2 and R3 (indexes, -l 8) and R4 (bio read, -l 101) as separate parallel jobs. slurm scripts (qual_R1.sh, qual_R2.sh, qual_R3.sh, qual_R4.sh,) found in /projects/bgmp/hodapp/bioinfo/Bi622/Demultiplex/Assignment-the-first/

slurm out info:
R1
```
WARN cache for Repodata at /home/hodapp/.cache/rattler/cache/repodata is on a network/parallel filesystem (NFS/SMB/FUSE/BeeGFS/Lustre/GPFS/CephFS), redirected to /tmp/pixi-cache-hodapp/repodata for this run. Set [cache.repodata] in config.toml or PIXI_CACHE_DIR to override, or [cache.netfs-redirect] = "never" to keep the original path.
363246735 records processed
wrote part1_output/R1_means.tsv and part1_output/R1_qual.png
	Command being timed: "pixi run python qual_dist.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -o part1_output/R1 -l 101"
	User time (seconds): 4404.01
	System time (seconds): 0.92
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:13:42
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 71572
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 23516
	Voluntary context switches: 2673
	Involuntary context switches: 10727
	Swaps: 0
	File system inputs: 0
	File system outputs: 16
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```
R2
```
 WARN cache for Repodata at /home/hodapp/.cache/rattler/cache/repodata is on a network/parallel filesystem (NFS/SMB/FUSE/BeeGFS/Lustre/GPFS/CephFS), redirected to /tmp/pixi-cache-hodapp/repodata for this run. Set [cache.repodata] in config.toml or PIXI_CACHE_DIR to override, or [cache.netfs-redirect] = "never" to keep the original path.
363246735 records processed
wrote R2_means.tsv and R2_qual.png
	Command being timed: "pixi run python qual_dist.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -o R2 -l 8"
	User time (seconds): 658.81
	System time (seconds): 0.19
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 11:01.95
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 72632
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 23339
	Voluntary context switches: 1263
	Involuntary context switches: 1562
	Swaps: 0
	File system inputs: 0
	File system outputs: 16
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```
R3
```
 WARN cache for Repodata at /home/hodapp/.cache/rattler/cache/repodata is on a network/parallel filesystem (NFS/SMB/FUSE/BeeGFS/Lustre/GPFS/CephFS), redirected to /tmp/pixi-cache-hodapp/repodata for this run. Set [cache.repodata] in config.toml or PIXI_CACHE_DIR to override, or [cache.netfs-redirect] = "never" to keep the original path.
363246735 records processed
wrote R3_means.tsv and R3_qual.png
	Command being timed: "pixi run python qual_dist.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -o R3 -l 8"
	User time (seconds): 669.13
	System time (seconds): 0.18
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 11:11.20
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 72788
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 23382
	Voluntary context switches: 112
	Involuntary context switches: 1584
	Swaps: 0
	File system inputs: 0
	File system outputs: 16
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```
R4
```
WARN cache for Repodata at /home/hodapp/.cache/rattler/cache/repodata is on a network/parallel filesystem (NFS/SMB/FUSE/BeeGFS/Lustre/GPFS/CephFS), redirected to /tmp/pixi-cache-hodapp/repodata for this run. Set [cache.repodata] in config.toml or PIXI_CACHE_DIR to override, or [cache.netfs-redirect] = "never" to keep the original path.
363246735 records processed
wrote R4_means.tsv and R4_qual.png
	Command being timed: "pixi run python qual_dist.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -o R4 -l 101"
	User time (seconds): 4443.22
	System time (seconds): 0.93
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:14:15
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 71488
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 22983
	Voluntary context switches: 113
	Involuntary context switches: 10404
	Swaps: 0
	File system inputs: 0
	File system outputs: 16
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```

Answered questions in Answers.md

to find indexes with N's did a sed to grab just the sequence lines (since that works on talapas unlike my mac, oh fun) starting line 2 every 4th line. then count with grep which ones have N. only need the two index files
```
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | sed -n '2~4p' | grep -c "N"
3976613

zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | sed 
-n '2~4p' | grep -c "N"
3328051
```
Ok after answering the Answers.md questions about quality cutoff I just cannot. This makes no sense Part 1 is so weird. Hamming distance is the right measure for if I even need a qual cutoff so I might as well just do it. Unless there is some thing in this actual set of indexes part 1 of this whole thing seems odd and not helpful to me. 
## Hamming distance of known indexes

Wrote [min_hamming.py](/projects/bgmp/hodapp/bioinfo/Bi622/Demultiplex/hamming_dist.py) to compute the minimum Hamming distance among the 24 known indexes. Reads the index sequences from indexes.txt (5th column), compares every unique pair with itertools.combinations, and then counts mismatched positions per pair, (wow zip is sure handy) and tracks the smallest minimum among the whole set by adding up each place where there is a mismatch (if x != y: then increase count 1) 
**Result:** Minimum Hamming distance = 3. 
This means the closest two indexes differ in 3 positions, so no single base error (or even 2!) can turn one valid index into another. Supports the argument that a per-base quality cutoff on indexes is unnecessary for correct demultiplexing and was maybe just for fun? Or confusion. There is no need for an quality cutoff score, anything with an N is getting put in unknown, and with a minimum hamming distance of 3, this is not going to be a reality enough to impact. 
## Pipeline / approach
- for demultiplexing.py script
	Step 1: open the four input files + the loop skeleton that reads records and stops.* Get the while loop going through all four files, grabbing records, stopping at the end. No sorting yet. Test: make it count how many records it read, confirm that matches test files (8 reads). This proves the reading record function works.
	
	Step 2: pull the indexes out and build the pair.** Inside the loop, get index1 from R2, index2 from revcomp'd R3, make the pair in the header format. Test: print them, compare to outputs of headers in test files.
	
	**Step 3: the sorting logic (the three buckets).** The if/elif/else that decides matched vs hopped vs unknown. Test: print which bucket each read lands in, check against test cases 
	
	**Step 4: build the modified header.** Stick the index-pair onto the R1 and R4 headers. Test: print, confirm the format looks right.
	
	**Step 5: actually write to files + tally.** Wire the sorted reads into the right file handles, add the counts. Test: look in test_output, confirm reads landed in the right files.
	
	**Step 6: close files + report.** Close all 52, print the summary? tbd on summary format and outputs
  

## Environment

* python 3.14.6
- matplotlib 3.11.1
- see pixi. files uploaded at repo https://github.com/palimpsestkenlyn/Demultiplex
	- 1bb097969d3f70c9872148b4320b03cfd6633b64 for pixi.toml

## Reference data

```
input files at: /projects/bgmp/shared/2017_sequencing/
Output files at: /scratch/bgmp/hodapp/ ** either test_output/ OR deumux_output/

-r1=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz \

-r2=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz \

-r3=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz \

-r4=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz \

-i=/projects/bgmp/shared/2017_sequencing/indexes.txt \

-o=/scratch/bgmp/hodapp/demux_output
```

  

---

  

## Log

### 8.2.26 — Starting actual code for demux, assignment the 3rd
- **Logic:** set up argparse: going to use 6 arguments or flags, to make as general as possible. lots of inputs, but more flexible and useful this way. so r1-r4 for input files, and then the index file as an input as well as the directory the outputs get written to.

### 8.4.26 — 
- **ran:** testing the barcode set, and file handle dictonaries (r1 and r2) are working and creating the files to write out to as expected. created directory in my scratch to test blank file creation to
- ```
  python demultiplexing.py -r1=/projects/bgmp/hodapp/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test_R1.fastq -r2=/projects/bgmp/hodapp/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test_R2.fastq -r3=/projects/bgmp/hodapp/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test_R3.fastq -r4=/projects/bgmp/hodapp/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test_R4.fastq -i=/projects/bgmp/hodapp/bioinfo/Bi622/Demultiplex/indexes.txt -o=/scratch/bgmp/hodapp/test_output
  ```


- **got:**  24 for length of barcodes set, and each r1 and r2 dict of file handles, as well as the expected 52 files in my scratch. Success! tried printing the file handle dictonaries and they look crazy. 
- ```
  ls /scratch/bgmp/hodapp/test_output
AACAGCGA_R1.fastq  AGAGTCCA_R2.fastq  ATCGTGGT_R1.fastq  CGATCGAT_R2.fastq  CTCTGGAT_R1.fastq  GATCTTGC_R2.fastq  GTCCTAAG_R1.fastq  TACCGGAT_R2.fastq  TCGACAAG_R1.fastq  TCGGATTC_R2.fastq  unknown_R1.fastq
AACAGCGA_R2.fastq  AGGATAGC_R1.fastq  ATCGTGGT_R2.fastq  CGGTAATC_R1.fastq  CTCTGGAT_R2.fastq  GCTACTCT_R1.fastq  GTCCTAAG_R2.fastq  TAGCCATG_R1.fastq  TCGACAAG_R2.fastq  TCTTCGAC_R1.fastq  unknown_R2.fastq
ACGATCAG_R1.fastq  AGGATAGC_R2.fastq  CACTTCAC_R1.fastq  CGGTAATC_R2.fastq  GATCAAGG_R1.fastq  GCTACTCT_R2.fastq  hopped_R1.fastq    TAGCCATG_R2.fastq  TCGAGAGT_R1.fastq  TCTTCGAC_R2.fastq
ACGATCAG_R2.fastq  ATCATGCG_R1.fastq  CACTTCAC_R2.fastq  CTAGCTCA_R1.fastq  GATCAAGG_R2.fastq  GTAGCGTA_R1.fastq  hopped_R2.fastq    TATGGCAC_R1.fastq  TCGAGAGT_R2.fastq  TGTTCCGT_R1.fastq
AGAGTCCA_R1.fastq  ATCATGCG_R2.fastq  CGATCGAT_R1.fastq  CTAGCTCA_R2.fastq  GATCTTGC_R1.fastq  GTAGCGTA_R2.fastq  TACCGGAT_R1.fastq  TATGGCAC_R2.fastq  TCGGATTC_R1.fastq  TGTTCCGT_R2.fastq
  ```

- **notes:** why does printing my file handle dictonaries look like this?
- ```
  print(r1_FHdict)
  {'TCTTCGAC': <_io.TextIOWrapper name='/scratch/bgmp/hodapp/test_output/TCTTCGAC_R1.fastq' mode='w' encoding='UTF-8'>, 'CACTTCAC': <_io.TextIOWrapper name='/scratch/bgmp/hodapp/test_output/CACTTCAC_R1.fastq' mode='w' encoding='UTF-8'>, 'TCGAGAGT': <_io.TextIOWrapper .....
  ```
- **Realizing I want another function, if I am going to create the 4 lines for each file, should be function**
	- Going to write a get_record function: no I am not, did, left it in draft. it confuses more than helps me. maybe if I was going to do this again but it feels more clear to me to actually write this section out, yes a bit repetitive but at least I can follow it.
- An assert for this function that I didn't put into bioinfo.py is weird don't want it to run every time. Did lots of print and testing on test files, but to confirm new function works made test_demultiplexing.py (/projects/bgmp/hodapp/bioinfo/Bi622/Demultiplex/Assignment-the-third/test_demultiplexing.py). I am going to call that good enough for a unit test, don't want to redo the current logic. 
- ```
  python test_demultiplexing.py 
get_record passed
  ```

### 8.4.26 part2 — 
#### Trying to approach the 4 files at once, different options
-  Overall logic problem: all 4 files information match ONLY by location (so read 1 is the 1st 4 lines in each file) so they have to be examined at the same time to keep things together...tricky!

 
	-  **for line in fh:** (the loop used before) Walks ONE file top to bottom, handing over one line at a time. Right tool when doing the SAME thing to every line in ONE file (counting quality, tallying k-mers, summing lengths). Problem here: a for-loop wants to own one file and run it to the end. Can't easily make four of them step forward together in sync. So it does not fit the four-files-at-once job. 
  
	-  **readline()**  Reads ONE line, then stops, and REMEMBERS where it stopped (same counter type thing seen in for blank in blank = iteration). Call it again -> next line. Call it 4 times in a row -> one full record, because records are 4 lines in fixed order. The order of the 4 calls is what puts the right line in each variable: 1st call = header, 2nd = sequence, 3rd = plus, 4th = quality. Never ask for "line 2" by number. Just keep asking for the NEXT one, the file tracks what next means. Here it fits. Pull one record from each of the four files, all four move forward one step together, repeat. allows each to be seperated into distinct objects which can then be used for things like appending header, etc...
  
	-  **zip(a, b, c, d)**  zip steps several things forward together so does offer like synchronization. Used it in the hamming code. But the catch for FASTQ: zip on files steps ONE LINE at a time, not one record. So each turn gives line 1 of each file, then line 2 of each, etc. Would have to track "am I on the header/seq/plus/qual line" and regroup lines back into records. Feels like a headache, and readline seems simpler. 

	-  **the takeaway:** 
		- one file, same thing to every line -> for line in fh: 
		- several files marched forward together -> readline() Picked readline for demux because it keeps four files in sync 
- 
#### Writing demultiplexing.py with some idea of logic now?
- Start writing script and testing: open 4 input files (gzip on all 4 test input files so I can from t he start be testing with gzip.open rt flag)
- Ran test with temporary print on how many records to see if get_record was working. yes got 8 and test files have 32 lines 8 reads. wrote while true loop with a break for empty r1_record at[0] ==WILL THIS WORK TO STOP ALL? I WOULD ASSUME BC THEY ARE HAPPENING Simultaneously?== 
- add to the loop, to build the index pairs in the format of the required header line, rev comp R3 seq.. test to print. Success! all 8 pairings look correct
- ```
  python demultiplexing.py -r1=/projects/bgmp/hodapp/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test_R1.fastq.gz -r2=/projects/bgmp/hodapp/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test_R2.fastq.gz -r3=/projects/bgmp/hodapp/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test_R3.fastq.gz -r4=/projects/bgmp/hodapp/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test_R4.fastq.gz -i=/projects/bgmp/hodapp/bioinfo/Bi622/Demultiplex/indexes.txt -o=/scratch/bgmp/hodapp/test_output
GTAGCGTA-GTAGCGTA
CGATCGAT-CGATCGAT
GTAGCGTA-GATCAAGG
CGATCGAT-AACAGCGA
NACCGGAT-TACCGGAN
NGTTCCGT-TGTTCCGN
AAAAAAAA-GTAGCGTA
AAAAAAAA-AAAAAAAA
  ```
- write an if/elif/else to sort into the 3 buckets membership, after thinking this through, I don't think I need a specific N check like in my psuedocode. its a bit like a riddle but if I just check for if index1 OR index2 are not in the barcode set, then it will automatically catch all the Ns anyways. bc no indexes in the barcode set have Ns. so skip the N check, and just 3 branches.
- temporary print statements to check, yay!
- ```
  python demultiplexing.py -r1=/projects/bgmp/hodapp/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test_R1.fastq.gz -r2=/projects/bgmp/hodapp/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test_R2.fastq.gz -r3=/projects/bgmp/hodapp/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test_R3.fastq.gz -r4=/projects/bgmp/hodapp/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test_R4.fastq.gz -i=/projects/bgmp/hodapp/bioinfo/Bi622/Demultiplex/indexes.txt -o=/scratch/bgmp/hodapp/test_output
matched GTAGCGTA-GTAGCGTA
matched CGATCGAT-CGATCGAT
hopped GTAGCGTA-GATCAAGG
hopped CGATCGAT-AACAGCGA
unknown NACCGGAT-TACCGGAN
unknown NGTTCCGT-TGTTCCGN
unknown AAAAAAAA-GTAGCGTA
unknown AAAAAAAA-AAAAAAAA
  ```


### 8.6.26 Statistics, write outputs and run script after testing

Close all files: use the FHdicts and values to close those, then manually list for hooped and unknown. Also close the input files. 

Write out counts/summary for full per-pair into to tsv so if there are display or reports later, don't need to rerun.

set counts, and add for the different categories to report on. do initial output markdown at least with summary of statistics. Use TSV file for more, graphs and such or another script if more time.

* tally totals from the counting dict for the summary
* have to use continue as the unknown has no two halves to compare, skip it  otherwise stalls when looking at separating out the two halves of indexes...
* just grab unknown directly no need to count, already total/value in dict
* write the summary to a markdown file (add two new lines so there are spaces between so its easier to read)
	* each total as a count and as a percent of all reads, rounded to 2 decimals\
	* per-sample breakdown: only matched pairs (halves equal) that actually got reads

- **ran:**
	demux.sh with:
	```
	/usr/bin/time -v python demultiplexing.py \

-r1=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz \

-r2=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz \

-r3=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz \

-r4=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz \

-i=/projects/bgmp/shared/2017_sequencing/indexes.txt \

-o=/scratch/bgmp/hodapp/demux_output
	```

* Usage summary:
	* time: 47:36.07
	* Exit status: 0
	* maximum resident set size: 250040
	* percent of CPU: 70%
 ```
 Command being timed: "python demultiplexing.py -r1=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -r2=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -r3=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -r4=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -i=/projects/bgmp/shared/2017_sequencing/indexes.txt -o=/scratch/bgmp/hodapp/demux_output"
        User time (seconds): 2000.49
        System time (seconds): 20.48
        Percent of CPU this job got: 70%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 47:36.07
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 250040
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 0
        Minor (reclaiming a frame) page faults: 33214
        Voluntary context switches: 50472
        Involuntary context switches: 12154
        Swaps: 0
        File system inputs: 0
        File system outputs: 0
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
        Page size (bytes): 4096
        Exit status: 0
 ```

NOTE FOR FUTURE: 
  multithread gzip:
	  pigz <files> NOT on a login node....
	  request our new nodes: so others are not using it if its full kick them off
		  #SBATCH --constraint=turin
	THIS will allow for the ability to create compressed output files so as not to take up so much space. Read up on pigz!

