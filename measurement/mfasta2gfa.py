#!/usr/bin/env python

"""
This script is derived from fasta2gfa.py to take as input multi-FASTA files (e.g., those downloaded from the
GenBank assembly database). The FASTA files should be named in the format of '[strain name].fasta' and there
must be only a single dot in the file name. For example, strain1.exp1.fasta is an illegal name for this script.

Usage: python fasta2gfa -i *.fasta -o references/gfa -d ''

Author: Yu Wan (wanyuac@gmail.com, https://github.com/wanyuac)
First edition and the latest edition: 18/7/2018
Python 2 and 3 compatible
License: GNU GPL 2.1
"""

from __future__ import print_function
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from argparse import ArgumentParser


def parse_arguments():
    parser = ArgumentParser(description = "")
    parser.add_argument("-i", nargs = "+", type = str, required = True, help = "Input multi-FASTA files")
    parser.add_argument("-o", type = str, required = False, default = "gfa", help = "Output directory")
    return parser.parse_args()


def main():
    args = parse_arguments()
    if not os.path.exists(args.o):
        os.system("mkdir " + args.o)  # set up the output directory
        
    # Create a GFA file for each FASTA file
    for f in args.i:
        create_gfa(f, args.o)
    
    return


def create_gfa(mfasta, outdir):
    strain, _ = os.path.splitext(os.path.basename(mfasta))
    gfa = open(os.path.join(outdir, strain + ".gfa"), "w")  # initiate the output file
    print("H\tVN:Z:1.0", file = gfa)  # print an universal header line into the output
    with open(mfasta, "rU") as handle:
        for entry in SeqIO.parse(handle, "fasta"):
            bp = len(entry.seq)
            print("\nS\t{id}\t{seq}\tLN:i:{seq_len}\tRC:i:{read_cov}".format(id = entry.id, seq = entry.seq, \
                  seq_len = str(bp), read_cov = str(bp * 10)), file = gfa)  # assuming an arbitrary ten-fold coverage
            print("L\t{id}\t+\t{id}\t+\t0M".format(id = entry.id), file = gfa)  # Assuming that there is no overlap between two ends.
    gfa.close()
    
    return


if __name__ == "__main__":
    main()
