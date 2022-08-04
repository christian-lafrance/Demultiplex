#!/bin/bash
#SBATCH --partition=bgmp       ### Partition (like a queue in PBS)
#SBATCH --job-name=demultiplex      ### Job Name
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node
#SBATCH --account=bgmp      ### Account used for job submission
#SBATCH --error=demultiplex.err          ### File in which to store job error messages

conda activate bgmp_py310

myDir="/projects/bgmp/clafranc/bioinfo/Bi622/Demultiplex/Assignment-the-third"
dataDir="/projects/bgmp/shared/2017_sequencing/"


/usr/bin/time -v ./demultiplex.py -pe True -r1 $dataDir/1294_S1_L008_R1_001.fastq.gz \
-r2 $dataDir/1294_S1_L008_R4_001.fastq.gz -i1 $dataDir/1294_S1_L008_R2_001.fastq.gz \
-i2 $dataDir/1294_S1_L008_R3_001.fastq.gz -index $dataDir/indexes.txt \
-q 30 -n 3 -o output_files/
