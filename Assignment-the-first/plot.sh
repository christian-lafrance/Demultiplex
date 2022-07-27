#!/bin/bash
#SBATCH --partition=bgmp       ### Partition (like a queue in PBS)
#SBATCH --job-name=i2_plot      ### Job Name
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node
#SBATCH --account=bgmp      ### Account used for job submission
#SBATCH --error=i2_plots.err          ### File in which to store job error messages

conda activate bgmp_py310

myDir="/projects/bgmp/clafranc/bioinfo/Bi622/Demultiplex/Assignment-the-first"
dataDir="/projects/bgmp/shared/2017_sequencing/"

echo "index 2"

/usr/bin/time -v $myDir/mean_phred.py -f $dataDir/1294_S1_L008_R3_001.fastq.gz -l 8 -p Index_2

