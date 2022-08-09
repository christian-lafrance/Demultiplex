#!/usr/bin/env python

import argparse
import bioinfo
import gzip
import numpy as np
from itertools import zip_longest
import matplotlib.pyplot as plt


def get_args():
    parser = argparse.ArgumentParser(description="demultiplex.py")
    parser.add_argument("-pe", help="Bool indicating paired-end reads", required=False, type=bool)
    parser.add_argument("-r1", help="Read 1 file", required=True, type=str)
    parser.add_argument("-r2", help="Read 2 file", required=False, type=str)
    parser.add_argument("-i1", help="index 1 file", required=True, type=str)
    parser.add_argument("-i2", help="index 2 file", required=True, type=str)
    parser.add_argument("-index", help="TSV file with expected indexes", required=True, type=str)
    parser.add_argument("-q", help="quality score cut off for index reads, default is 30", required=False, type=int, default=30)
    parser.add_argument("-n", help="Number of bases allowed to have quality score below the threshold. Default is 3.", required=False, type=int, default=3)
    parser.add_argument("-o", help="directory name for output files", required=True, type=str)
    return parser.parse_args()
args = get_args()


def generate_barcode_list(index_file: str) -> set:
    '''
    Takes a TSV file name as an argument. Opens the file that contains 
    expected barcodes and generate a list of the expected barcodes. This 
    list will be searched by the generate_output_files function to check 
    if the current barcode is expected. 
    '''
    barcodes = set() # stores each valid barcode as a string. 

    with open(index_file, "r") as f:
        header = True
        for line in f:
            if header == False:
                current_line = line.split()
                barcode = current_line[-1].strip("\n")
                barcodes.add(barcode)
            header = False

    return barcodes


def generate_output_files(barcodes: set, is_pairedend: bool) -> dict:
    '''
    Takes a list of expected barcodes as an input and generates the
    expected files. One file per barcode and read. Includes unknown file
    and hopped index file per read. Total of 52 files. Generates a dictionary
    with a string of the filename as a key and the filehandle as a value. 
    Use this to refer to the file. 
    '''
    # generate output files:
    output_files_dict = {} # dictionary to store file name as key and handle as value. 
    for barcode in barcodes:
        key = f"{barcode}_R1.fq"
        output_files_dict[key] = open(f"{args.o}{barcode}_R1.fq", "w")

        if is_pairedend:
            key = f"{barcode}_R2.fq"
            output_files_dict[key] = open(f"{args.o}{barcode}_R2.fq", "w")
    
    output_files_dict["hopped_R1.fq"] = open(f"{args.o}hopped_R1.fq", "w")
    output_files_dict["unknown_R1.fq"] = open(f"{args.o}unknown_R1.fq", "w")

    if is_pairedend:
        output_files_dict["hopped_R2.fq"] = open(f"{args.o}hopped_R2.fq", "w")
        output_files_dict["unknown_R2.fq"] = open(f"{args.o}unknown_R2.fq", "w")

    return output_files_dict


def write_record(prefix: str, index1_seq: str, index2_seq: str, 
                current_records: dict, output_files_dict: dict) -> None:
    '''
    General function used to write files based on the current record and
    the index sequences. Returns None. 
    '''
    # update headers with indexes
    current_records["r1"][0] += f" {index1_seq}-{index2_seq}"
    current_records["r2"][0] += f" {index1_seq}-{index2_seq}"

    # write current records to unknown fq files. 
    output_files_dict[f"{prefix}_R1.fq"].write("\n".join(current_records["r1"]))
    output_files_dict[f"{prefix}_R1.fq"].write("\n")
    output_files_dict[f"{prefix}_R2.fq"].write("\n".join(current_records["r2"]))
    output_files_dict[f"{prefix}_R2.fq"].write("\n")

    return


def parse_and_write_files(barcodes, output_files_dict) -> tuple[dict:dict]:
    '''
    Need to update to make functional with is_pairedend bool.

    Opens all read, index, and output files, parses the input files,
    and checkes for hopped, unknown, and dual-matched indexes. Writes
    records to the appropriate output file. 

    Returns two dictionaries: 
    stats_dict which includes stats on hopped, unknown, and matched records.
    hopped_dict which includes hopped index pairs and numnber of occurences.
    '''
    with gzip.open(args.r1, "rt", encoding="utf-8") as r1_file, gzip.open(args.r2, "rt", encoding="utf-8") as r2_file, gzip.open(args.i1, "rt", encoding="utf-8") as i1_file, gzip.open(args.i2, "rt", encoding="utf-8") as i2_file:
        
        line_counter = 0 # counts 0-4 to assemble each record one at a time in memory. 
        hopped_dict = {} # stores hopped indexes and counts: index1-index2: count.
        stats_dict = {}  # stores counts for hopped, unknown, and mapped records. 
        matched_count = {} # stores counts of dual-matched records: index1 : count. 
        current_records = {"r1": ["", "", "", ""], "r2": ["", "", "", ""], 
                            "i1": ["", "", "", ""], "i2": ["", "", "", ""]} # stores the current record. 
        
        for r1, r2, i1, i2 in zip(r1_file, r2_file, i1_file, i2_file):
            
            # Assemble the current record using the line counter. 
            current_records["r1"][line_counter] = r1.strip("\n") 
            current_records["r2"][line_counter] = r2.strip("\n")
            current_records["i1"][line_counter] = i1.strip("\n")
            current_records["i2"][line_counter] = i2.strip("\n")

            line_counter += 1
        
            if line_counter == 4:
                '''
                The next record is fully assembled. Reset line_counter, filter
                for quality, check index pairs, and write to files. Reset current
                record dictionary at the end. 
                dict structure example:
                {"i1":[header, seq, +, qscores]}
                '''
                line_counter = 0

                # check if indexes are valid:
                index1_seq = current_records["i1"][1]
                index2_seq = bioinfo.rev_comp(current_records["i2"][1])

                if "N" in index1_seq or "N" in index2_seq:
                    '''Write to unknown file'''
                    write_record("unknown", index1_seq, index2_seq, current_records, output_files_dict)

                    # collect unknown stats. 
                    if "unknown_indexes" in stats_dict:
                        stats_dict["unknown_indexes"] += 1
                    else:
                        stats_dict["unknown_indexes"] = 1

                elif index1_seq in barcodes and index2_seq in barcodes:
                    if index1_seq != index2_seq:
                        '''They are hopped. Write to hopped file. '''
                        write_record("hopped", index1_seq, index2_seq, current_records, output_files_dict)

                        # collect hopped stats. 
                        if "hopped_indexes" in stats_dict:
                            stats_dict["hopped_indexes"] += 1 
                        else:
                            stats_dict["hopped_indexes"] = 1

                        key = f"{index1_seq}-{index2_seq}"
                        if key in hopped_dict:
                            hopped_dict[key] += 1
                        else:
                            hopped_dict[key] = 0

                    else:
                        '''
                        Index pairs are expected and not hopped. Filter for quality score. 
                        If below the specified quality cut off (or default of 30) write 
                        to unknown file. If above the cut off, then write to designated
                        file with index in the name. Will only tolerate a certain number
                        of bases below the specified threshold before writting to the
                        unknown file. This number of bases can be configured with -n. 
                        '''
                        # filter for index seq containing N and qual score:
                        # by iterating over base
                        index1_qual = current_records["i1"][3]
                        index2_qual = current_records["i2"][3]

                        for i, base in enumerate(index1_seq): 
                            low_qual_base_count = 0
                            i1_score = bioinfo.convert_phred(index1_qual[i])
                            i2_score = bioinfo.convert_phred(index2_qual[i])
                            if i1_score < args.q or i2_score < args.q:
                                low_qual_base_count += 1
                            if low_qual_base_count > args.n:
                                write_record("unknown", index1_seq, index2_seq, current_records, output_files_dict)                                
                                if "unknown_indexes" in stats_dict:
                                    stats_dict["unknown_indexes"] += 1
                                else:
                                    stats_dict["unknown_indexes"] = 1
                                break

                            elif i == len(index1_seq)-1:
                                '''All scores are above the threshold. Write reads to dual-mapped files.'''
                                write_record(index1_seq, index1_seq, index2_seq, current_records, output_files_dict)

                                if index1_seq in matched_count:
                                    matched_count[index1_seq] += 1
                                else:
                                    matched_count[index1_seq] = 1

                else:
                    '''
                    At least one index is not expected. Update headers and write
                    to unknown files. Use write_record function. 
                    '''
                    write_record("unknown", index1_seq, index2_seq, current_records, output_files_dict)

                    # collect unknown stats. 
                    if "unknown_indexes" in stats_dict:
                        stats_dict["unknown_indexes"] += 1
                    else:
                        stats_dict["unknown_indexes"] = 1

    # close files when done
    for file in output_files_dict:
        output_files_dict[file].close()

    return stats_dict, hopped_dict, matched_count


def generate_user_report(summary_stats: dict, hopped_dict: dict, matched_count: dict, barcode_set: list) -> None:
    '''
    Generates Demultiplex_summary.tsv which includes number of hopped indexes,
    number of unknown indexes, number of dual-matched indexes, and their
    percentages. Also generates a heatmap of hopped indexes and the number
    of occurences. 
    '''
    # plot heatmap
    barcodes = sorted(list(barcode_set))
    l = len(barcodes)
    values = [[0 for i in range(l)] for j in range(l)] 

    for k, val in hopped_dict.items():
        axis = k.split("-")
        
        i1 = axis[0]
        i2 = axis[1]
        
        v = barcodes.index(i1)
        h = barcodes.index(i2)
            
        values[v][h] = val

    fig, ax = plt.subplots(figsize= (15, 15))
    im = ax.imshow(values)

    ax.set_xticks(np.arange(l), labels=barcodes)
    ax.set_yticks(np.arange(l), labels=barcodes)

    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    for i in range(l):
        for j in range(l):
            text = ax.text(j, i, values[i][j], ha="center", va="center", color="w")

    ax.set_title("Heatmap of Index Hopping")
    plt.savefig("Index_hopping_heatmap.png")


    '''write summary file'''
    total = 0
    matched = 0
    for key, value in stats_dict.items():
        total += value
    for value in matched_count.values():
        total += value
        matched += value

    try:
        hopped = stats_dict["hopped_indexes"]
    except KeyError:
        hopped = 0

    hopped_perc = round(hopped/total, 3)
    unknown = stats_dict["unknown_indexes"]
    unknown_perc = round(unknown/total, 3)

    with open("Demultiplex_summary.txt", "w") as out:
        out.write("=== User parameters for barcode identification ===\n")
        out.write(f"quality_cutoff: {args.q}\n")
        out.write(f"# of bases allowed below cutoff: {args.n}\n\n")
        out.write("========== Record counts ==========\n")
        for name, count in stats_dict.items():
            out.write(f"{name}: {count}\n")
        out.write(f"dual_matched_count: {matched}\n")
        out.write(f"total_records: {total}\n\n")
        out.write("========== Dual-matched ==========\n")
        out.write("Index:\t\t%_of_total\t%_of_matched\n")
        for name, count in matched_count.items():
            out.write(f"{name}:\t\t{round(count/total*100, 3)}\t\t{round(count/matched*100, 3)}\n")

        out.write("\n\n========== Hopped indexes ==========\n")
        out.write("Index pair: count\n")

        '''
        Include hopped indexes and counts in the report. Need to sort by the
        values of hopped_dict. To do this, iterate over the dictionary, creating
        a list of tuples with keys and values swapped. Sorted will then sort
        by the counts. 
        '''
        hopped_list = []
        for k, v in hopped_dict.items():
            hopped_list.append((v, k))
        for item in sorted(hopped_list, reverse=True):
            indexes = item[1]
            count = item[0]
            out.write(f"{indexes}: {count}\n")

    return


# FUNCTION CALLS
is_pairedend = args.pe

barcodes = generate_barcode_list(args.index)
output_files_dict = generate_output_files(barcodes, is_pairedend)
stats_dict, hopped_dict, matched_count = parse_and_write_files(barcodes, output_files_dict)
generate_user_report(stats_dict, hopped_dict, matched_count, barcodes)

