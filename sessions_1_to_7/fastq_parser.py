import numpy as np
import pandas as pd
from itertools import zip_longest

def parse_fastq(filepath):
    sequences =[]
    qualities = []

    with open(filepath,"r") as f:
        while True:
            header = f.readline().strip()
            if not header:
                break


            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()

            sequences.append(seq)
            qualities.append(qual)

    return sequences, qualities

def phred_scores(quality_string):
    return [ord(char)-33 for char in quality_string]

def base_comp_per_cycle(filepath):
    sequences =[]
    results = []

    with open(filepath,"r") as f:
        while True:
            header = f.readline().strip()
            if not header:
                break


            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()

            sequences.append(seq)
    
    padded = list(zip_longest(*sequences, fillvalue="N"))
    
    print(padded)

    for cycles in padded:
        counts = {
            "A":cycles.count("A"),
            "T":cycles.count("T"),
            "G":cycles.count("G"),
            "C":cycles.count("C"),
            "N":cycles.count("N")
        }
        results.append(counts)

    return results
    
def gc_computation(sequences):
    gcresults = []
    for reads in sequences:
        counts = {
            "G":reads.count("G"),
            "C":reads.count("G"),
        }
        total = len(reads)
        gcperc = (((counts["G"])+(counts["C"]))/total)*100
        gcresults.append(gcperc)

    return gcresults




