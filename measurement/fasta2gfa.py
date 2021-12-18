#!/usr/bin/env python

"""
This script converts every nucleotide sequence (FASTA format) into a circular topology in the GFA format that can be processed
and visualised using Bandage. It is particularly useful in measuring the shortest physical distances between genes in
finished-grade reference bacterial genomes which are often circular but stored as linear sequences in FASTA files.

This script expects every input file name is encoded by the formula [strain name]__[accession number/sequence name].[fasta/fna].
It assumes every sequence in the FASTA files are of circular topology in living organisms. Currently, this script only creates a
single GFA file for all sequences of each strain. For example, circular graphs of the chromosome and plasmids of the same strain
will be stored in the same GFA file for the simplicity of measuring physical distances between genes.

Usage: python fasta2gfa -i *__*.fasta -o references/gfa -d '__'

Author: Yu Wan (wanyuac@gmail.com, https://github.com/wanyuac)
First edition: 19/11/2016; the latest edition: 18/12/2021
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
from collections import defaultdict


def parse_arguments():
    parser = ArgumentParser(description = "")
    parser.add_argument("-i", nargs = "+", type = str, required = True, help = "Input FASTA files")
    parser.add_argument("-o", type = str, required = False, default = "gfa", help = "Output directory")
    parser.add_argument("-d", type = str, required = False, default = "__", \
                        help = "Delimiter between two fields in each encoded input file name")
    return parser.parse_args()


def main():
    args = parse_arguments()
    if not os.path.exists(args.o):
        os.system("mkdir " + args.o)  # set up the output directory
    file_sets = parse_filenames(args.i, args.d)
    make_gfa(file_sets, args.o)
    return


def parse_filenames(inputs, delimiter):
    # e.g., NJST258_1__CP006923.fasta => NJST258_1, CP006923
    file_sets = defaultdict(dict)
    for f in inputs:
        base = os.path.splitext(os.path.basename(f))[0]  # omit the filename extension
        strain, accession = base.split(delimiter)
        file_sets[strain][accession] = f
    return file_sets


def make_gfa(file_sets, outdir):
    for strain, files in file_sets.items():
        gfa = open(os.path.join(outdir, strain + ".gfa"), "w")  # make a GFA file for every strain
        print("H\tVN:Z:1.0", file = gfa)  # print a header line that contains a GFA version number to the current output file
        for accession, f in files.items():
            convert_fasta(accession, f, gfa)  # process each file
        gfa.close()
    return


def convert_fasta(accession, fasta, gfa):
    c = 0  # an additional counter for sequences in each FASTA file
    with open(fasta, "r") as handle:
        for genome in SeqIO.parse(handle, "fasta"):  # go through every record (assuming to be a complete circular genome) in the current FASTA file
            c += 1
            bp = len(genome.seq)
            if c == 1:
                seq_id = accession
            else:
                seq_id = accession + "." + str(c)  # in case a FASTA file contains more than a single sequence
            print("\nS\t{id}\t{seq}\tLN:i:{seq_len}\tRC:i:{read_cov}".format(id = seq_id, seq = genome.seq, \
                  seq_len = str(bp), read_cov = str(bp * 10)), file = gfa)  # assuming a ten-fold coverage
            print("L\t{id}\t+\t{id}\t+\t0M".format(id = seq_id), file = gfa)  # Assuming that there is no overlap between two ends.
    return

    
if __name__ == "__main__":
    main()
