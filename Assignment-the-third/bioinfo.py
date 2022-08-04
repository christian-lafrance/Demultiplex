#!/usr/bin/env python

# Author: Christian La Franmce clafranc@uoregon.edy
# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - 
#https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module
'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''
__version__ = "0.6"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning
DNA_bases = "ATGCN"
RNA_bases = "AUGCN"

def convert_phred(letter:str) -> int:
    '''
    Takes as an argument a single ASCII character and 
    returns the corresponding Phred 33 score. 
    '''
    score = ord(letter) - 33
    return score


def qual_score(phred_score: list) -> float:
    """
    Takes as an argument a list of ASCII characters, converts 
    each into the corresponding numerical score and returns 
    their average.
    """
    total = 0
    for score in phred_score:
        total += convert_phred(score)
        
    return total/len(phred_score)


def validate_base_seq(seq:str, RNA=False) -> bool:
    '''
    Takes as an argument a sequence and returns a bool, True if
    the sequence contains valid bases and False if not. An optional
    second argument specifies if the input sequence is RNA (bool).
    Default is False.
    '''

    seq = seq.upper()

    return len(seq) == seq.count("A") + seq.count("U" if RNA else "T") + seq.count("G") + seq.count("C") 



def gc_content(seq:str) -> float:
    '''
    Takes as an argument a string of DNA bases and calculates
    the GC content. 
    '''
    seq = seq.upper()

    gc_content = (seq.count("G") + seq.count("C"))/len(seq)

    return gc_content


def oneline_fasta(fasta_file: str):
    '''
    Accepts as an argument a fasta file name as a string and
    converts it to a fasta file without newline characters in
    the sequences. The oneline fasta file will be named
    oneline_{fasta_file}.
    '''
    with open(f"oneline_{fasta_file}", "w") as out:
        with open(fasta_file, "r") as fa:
            first_line = True
            for line in fa:
                if line[0] == ">":
                    if first_line == False:
                        out.write("\n")
                    out.write(line)
                else:
                    line = line.strip("\n")
                    out.write(line)
                    first_line = False

def rev_comp(seq: str) -> str:
    '''
    Accepts as an argument a string, checks if it contains valid nucleic
    acid characters, and returns the reverse complement of the sequence. 
    Works for both DNA and RNA. 
    '''
    upper_seq = seq.upper()
    complement = {"A":"T", "T":"A", "G":"C", "C":"G", "U":"A", "N":"N"}

    rcomp = "" # the reverse complement of seq. 

    for base in upper_seq:
        if base in DNA_bases or base in RNA_bases:
            rcomp += complement[base]
        else:
            raise Exception(f"Sequence contains invalid base: \"{base}\"")

    return rcomp[::-1]

if __name__ == "__main__":
    assert validate_base_seq("AATAGAT", False) == True
    assert validate_base_seq("AAUAGAU", True) == True
    assert validate_base_seq("TATUC",False) == False
    assert validate_base_seq("UCUGCU", False) == False

    assert gc_content("GTCA") == 0.5
    assert gc_content("ATTC") == 0.25

    assert convert_phred("E") == 36

    assert qual_score(["E", "A"]) == 34.0

    assert rev_comp("ATGCN") == "NGCAT"

    
