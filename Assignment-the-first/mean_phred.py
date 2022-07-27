#!/usr/bin/env python

# srun --account=bgmp --partition=bgmp --nodes=1 --ntasks-per-node=1 --time=2:00:00 --cpus-per-task=1 --pty bash

import bioinfo 
import matplotlib.pyplot as plt
import argparse
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="mean_phred.py")
    parser.add_argument("-f", help="fastq file name", required=True, type=str)
    parser.add_argument("-l", help="read length", required=True)
    parser.add_argument("-p", help="plot name without .png", required=True, type=str)
    return parser.parse_args()
args = get_args()

read_len = int(args.l)


def mean_phred(file: str, read_len: int) -> list:
    '''
    Accepts as an argument the file name of a FASTQ file. 
    Parses the FASTQ file and extracts the quality scores. 
    Converts the quality scores to numerical values and
    returns a list of mean quality score at each read
    position. 
    '''

    record = 0 # stores the number of the current record in 0 base counting. 
    with gzip.open(file, "rt", encoding="utf-8") as fq:

        qual_totals = [0 for i in range(read_len)]

        count = 1 # stores the current line in the file, 1-4. 

        for line in fq:
            if count == 4:
                line = line.strip("\n")
                for pos, score in enumerate(line):
                    qual_totals[pos] += bioinfo.convert_phred(score)                            

                record += 1            
            count += 1

            if count == 5:
                count = 1
                
    qual_means = [i/record for i in qual_totals]

    return qual_means


def generate_plot(qual_means: list, read_len: int) -> None:
    '''
    Accepts as an argument a list of mean quality
    scores. The index of each quality score in the list
    should correspond to the read position. 
    Generates and saves a plot of the mean quality score
    at each read position as a png file. 
    '''

    pos = [i for i in range(read_len)]

    fig, ax = plt.subplots()
    ax.bar(pos, qual_means)
    ax.set_ylabel("Mean Quality Score")
    ax.set_xlabel("Read Position")
    ax.set_title(f"Mean Quality Scores for {args.p}")

    plt.savefig(args.p)


qual_means = mean_phred(args.f, read_len)
generate_plot(qual_means, read_len)

        


