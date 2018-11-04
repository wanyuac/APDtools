#!/usr/bin/env python

'''
This script randomly draw N sequences from a multi-FASTA file. A FASTA
file of the same name as the input file will be saved in a destination
directory as the output of this script.

Usage: python sample_seqs.py cds.fna 1000 ./Output

Python 2 and 3 compatible.

Author: Yu Wan (wanyuac@gmail.com, GitHub: https://github.com/wanyuac)
First edition: 29 Nov 2015
Last edition: 3 Nov 2018
License: GNU GPL 2.1
''' 

from __future__ import print_function
import os
import sys
import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from argparse import ArgumentParser

'''
def draw_seqs(seqs, n):
    new_record = []
    selected_indices = random.sample(range(0, len(seqs)), n)
    i = 0
    for record in seqs:
        if i in selected_indices:
            new_record.append(record)
        i += 1
    return new_record
'''

def main():
    #files = sys.stdin.readlines()  # read a list of file names
    fasta = sys.argv[1]
    n = int(sys.argv[2])  # expected number of sequences per FASTA file
    outdir = sys.argv[3]
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    print("Sampling %d sequences from %s..." % (n, fasta))
    
    #subsampled.seqs = draw_seqs(list(SeqIO.parse(fasta, "fasta")), n)
    sampled_seqs = random.sample(list(SeqIO.parse(fasta, "fasta")), n)
        
    # Save the output file under the same name in the directory outdir.
    SeqIO.write(sampled_seqs, os.path.join(outdir, os.path.basename(fasta)), "fasta")

########## main program ##########
if __name__ == "__main__":
    main()