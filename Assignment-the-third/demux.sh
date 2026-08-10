#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --job-name=demux
#SBATCH --output=demux_%j.out

/usr/bin/time -v python demultiplexing.py \
-r1=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz \
-r2=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz \
-r3=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz \
-r4=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz \
-i=/projects/bgmp/shared/2017_sequencing/indexes.txt \
-o=/scratch/bgmp/hodapp/demux_output