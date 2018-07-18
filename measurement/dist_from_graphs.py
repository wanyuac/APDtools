#!/usr/bin/env python

"""
This script runs Ryan"s Bandage recursively to pull out physical distances between query sequences for a set of samples.
It submits a job through SLRUM for each sample. Note that every sample must have a single pair of graph/FASTA file and a
FASTA file of query sequences.

Usage:
python dist_from_graphs.py --graphs *_spades.fastg --queries *_kp.all_consensus_alleles.fasta --filter sample_alleles.tsv \
--suffix_graphs _spades.fastg --suffix_queries _kp.all_consensus_alleles.fasta --suffix_out spd \
--bandage /vlsci/VR0082/shared/ryan/Bandage/Bandage --mem 8192 --debug

Format of the filter TSV file:[sample name]\t[allele1,...,alleleN]. Specifically, allele names are comma-delimited. For instance:
    strain1\ta11,a12,a13\n
    strain2\ta21,a22\n
    strain3\ta31,a32,a33,a34,a35\n
    ...

Python version 2 and 3 compatible
Development history: 18-22/11/2016, 21-22/5/2017, ...; the latest edition: 8/1/2018

Copyright 2017-2018 Yu Wan <wanyuac@gmail.com>
Licensed under the Apache License, Version 2.0
"""

from __future__ import print_function
import os
import sys
import time
from collections import namedtuple
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

"""
Global constants: currently, this script is designed to work on the computational cluster Snowy of VLSCI. Nonetheless,
it is easy to adapt this script to other SLURM-running systems through configuring the following constants. Set BLAST_MODULE
to "" when BLAST is not running as a module on your computer.
"""

BLAST_MODULE = "BLAST+/2.2.31-GCC-4.9.2"


def parse_args():
    parser = ArgumentParser(description = "Submit jobs of Bandage to measure physical distances between genes")
    parser.add_argument("--graphs", nargs = "+", type = str, required = True, help = "Assembly graphs in FASTG or GFA formats")
    parser.add_argument("--queries", nargs = "+", type = str, required = True, \
                        help = "FASTA files containing query nucleotide sequences")
    parser.add_argument("--filter", type = str, required = False, default = None, \
                        help = "A tab-delimited file indicating which samples and alleles will be used")
    parser.add_argument("--ext_id", action = "store_true", required = False, help = "Flag it when your allele IDs contain any extended IDs, such as '.125'")
    parser.add_argument("--ext_id_delim", type = str, required = False, default = ".", \
                        help = "Delimiter leading extended IDs. It must not be an underscore.")
    parser.add_argument("--suffix_graphs", type = str, required = False, default = "_spades.fastg", \
                        help = "Suffix that will be used to extract sample names from file names")
    parser.add_argument("--suffix_queries", type = str, required = False, default = "_consensus_alleles.fasta", \
                        help = "Suffix that will be used to extract sample names from file names")
    parser.add_argument("--suffix_out", type = str, required = False, default = "spd", \
                        help = "Suffix that will be appended to sample names to form an output file name.")
    parser.add_argument("--outdir", type = str, required = False, default = "distances", help = "Output directory")
    parser.add_argument("--bandage", type = str, required = False, default = "Bandage", help = "Path to Bandage")
    parser.add_argument("--mem", type = str, required = False, default = "8192", help = "Memory (Mb) for every job.")
    parser.add_argument("--walltime", type = str, required = False, default = "1-0:0:00", help = "The walltime of every task")
    parser.add_argument("--par", type = str, required = False, default = "main", help = "Partition specifing which queue will the jobs belong to")
    parser.add_argument("--delay", type = int, required = False, default = 5, help = "Number of seconds before submitting the next job. Set to 0 to turn the delay off.")
    parser.add_argument("--other_args", type = str, required = False, default = "", \
                        help = "Other arguments that will be directly passed to Bandage.")
    parser.add_argument("--debug", action = "store_true", help = "Set it to prohibit submission of jobs.")
    return parser.parse_args()


def main():
    args = parse_args()

    # prepare the output directory
    outdir = os.path.realpath(args.outdir)
    if not os.path.exists(outdir):
        print("Making the output folder.")
        os.system("mkdir " + outdir)
    
    # process sequence files
    files = organise_files(graph_files = args.graphs, query_files = args.queries, sf_g = args.suffix_graphs, sf_q = args.suffix_queries)
    
    # Reformat sequence headers and filter query sequences
    if args.filter != None:  # Import filter information when it is provided and filter query sequences based on the filter
        filter = setup_filter(args.filter)
        if len(filter) > 0:  # when a filter is loaded. Otherwise, no changes apply to query sequences.
            files = filter_queries(files, filter, outdir, args.ext_id, args.ext_id_delim)    
    
    submit_jobs(files = files, bandage = args.bandage, mem = args.mem, walltime = args.walltime, par = args.par, \
                other_args = args.other_args, outdir = outdir, sf_out = args.suffix_out, run = not args.debug, \
                delay = args.delay)
    
    return


def setup_filter(tsv):
    """
    Read a TSV file that follows the format: [sample]\t[allele1,allele2,...,alleleN].
    Output: filter = {strain1 : [allele1, ..., alleleN], strain2 : [...], ...}
    Specifically, query sequence files will be filtered in accordance with the filter when it is configured. Otherwise, an
    all-by-all pairwise distance measurement will be carried out by Bandage for every sample. Note that in the first scenario,
    this script tries to convert allele names in the TSV file back to what are stored in the query FASTA file because SRST2
    changes allele names slightly when it is compiling allele calls into a single table.
    """
    print("Reading filter information.")
    filter = {}
    with open(tsv, "rU") as f:
        for line in f:
            sample, alleles = line.rstrip("\n").split("\t")
            alleles = alleles.split(",")  # Allele IDs may have extended identifiers attached. For instance, AadA1-pm_1597.1287, where "1287" is the extended ID separated by a full stop.
            if len(alleles) > 1:                    
                filter[sample] = reformat(alleles)  # transform allele identifiers into the same form in an SRST2-compatible genetic database
            else:
                print("Sample " + sample + " is skipped because it has only a single query sequence.")  # Otherwise, Bandage makes this job fail.
    
    return(filter)


def reformat(ids):
    """
    Reformat allele identifiers in a list to be the same as those in an SRST2-compatible genetic database.
    For example, Aac6-Ib_1677 will be converted into Aac6-Ib__1677. This function deals with the problem that SRST2 shortens
    allele names to make allele calls. The only argument of this function is a list of allele names from SRST2's compiled
    allele calls.
    """
    new_ids = []
    for a in ids:
        indices = find_char(a, "_")  # get all indices of the character "_", which separates an allele name with its sequence identifier
        if len(indices) > 0:
            i = indices[-1]  # only substitute the last "_" character with "__"
            new_ids.append(a[ : i] + "__" + a[i + 1 : ])
        else:
            new_ids.append(a)
            
    return new_ids


def find_char(string, char):
    # Returns all indices of a character in a string. Note that a string starts at the position zero.
    indices = []
    for i, c in enumerate(string):
        if c == char:
            indices.append(i)
            
    return indices


def filter_queries(files, filter, outdir, is_ext, ext_delim):
    # configure a folder for filtered query sequences
    new_fasta_dir = os.path.join(outdir, "filtered_queries")
    if not os.path.exists(new_fasta_dir):
        os.system("mkdir " + new_fasta_dir)
        
    print("Filtering query sequences based on the filter.")
    targeted_samples = list(filter.keys())  # We only need to measure the distances in samples that habours query sequences.
    new_files = {}
    File_set = namedtuple("File_set", ["graph", "query"])  # for filtered and reformatted query sequences of each sample
    for sample, file_set in files.items():
        if sample in targeted_samples:  # Samples on the file list can be less than those in the filter.
            written_count = 0
            targeted_alleles_mapping = rm_extended_IDs(filter[sample], is_ext, ext_delim)  # remove extended IDs from allele IDs
            targeted_alleles = list(targeted_alleles_mapping.keys())
            new_fasta = os.path.join(new_fasta_dir, sample + "__querySeqs.fasta")
            f = open(new_fasta, "w")
            with open(file_set.query, "rU") as query_file:
                for seq in SeqIO.parse(query_file, "fasta"):
                    allele_id = extract_allele_id(seq.id)
                    if allele_id in targeted_alleles:
                        seq.id = targeted_alleles_mapping[allele_id].replace("__", "_")  # retrieve the original allele ID to make a new sequence name
                        seq.description = ""
                        SeqIO.write(seq, f, "fasta")
                        written_count += 1
            f.close()
            new_files[sample] = File_set(graph = file_set.graph, query = new_fasta)
            
            # check if all targeted alleles have been accessed
            n = len(filter[sample])  # number of unfiltered alleles in the current sample
            if written_count < n:  # which should not happen, but sometimes R replaces dashes with full stops
                print("Warning: there are %i alleles in the sample %s failed to find their corresponding query sequences." % \
                      (n - written_count, sample))
    
    return new_files


def rm_extended_IDs(ids, is_ext, delim):  # expect ids to be a list
    new_ids = {}  # {new name : original name}
    if is_ext:
        for id in ids:
            new_id = id.split(delim)[0]
            new_ids[new_id] = id
    else:
        for id in ids:
            new_ids[id] = id  # just copy this ID
            
    return new_ids


def extract_allele_id(header):
    """
    This function extracts the allele ID from a header line. For example, the allele ID "TEM-198__1035" is extracted from
    the header "205__TEM-1D_Bla__TEM-198__1035.consensus ERR720250" or "205__TEM-1D_Bla__TEM-198__1035.consensus|ERR720250".    
    """
    db_entry = header.split(".")[0]  # restores the sequence header in the SRST2-formatted genetic database
    try:
        allele_id = "__".join(db_entry.split("__")[-2 : ])  # take the last two elements. For instance, "OqxB" + "__" + "49" => "OqxB__49".
    except IndexError:
        print("Error: the header line of the query sequence " + db_entry + " is incompatible with SRST2.")
        raise
    
    return allele_id


def organise_files(graph_files, query_files, sf_g, sf_q):
    """
    This function returns a one dimensional dictionary in which each element is a named tuple, namely,
    {sample name : (graph = a FASTG/GFA/FASTA file, query = a query FASTA file)].
    """
    print("Organising graph and query files.")
    graphs = make_dict(graph_files, sf_g)  # first, import all graph files
    queries = make_dict(query_files, sf_q)  # import query files
    
    # merge two dictionaries into a list to named tuples
    # So graph_files and query_files do not need to completely match.
    samples = lists_intersect(list(graphs.keys()), list(queries.keys()))  # find out common samples
    File_set = namedtuple("File_set", ["graph", "query"])
    files = {}  # sample name: namedtuple
    for sample in samples:
        files[sample] = File_set(graph = graphs[sample], query = queries[sample])
    
    return files


def make_dict(files, suffix):
    d = {}
    for f in files:
        try:
            sample = os.path.basename(f).rstrip(suffix)  # E.g., ERR562360 will be extracted from the string "~/assemblies/ERR562360_spades.fastg".
            d[sample] = f
        except KeyError:
            print("Error: invalid sample name " + sample + " among graph files.")
            raise
        
    return d


def lists_intersect(a, b):
    # takes an intersection between two lists
    return list(set(a) & set(b))

  
def submit_jobs(files, bandage, mem, walltime, par, other_args, outdir, sf_out, run, delay):
    print("Submitting jobs of Bandage for " + str(len(files)) + " samples.")
    to_delay = delay > 0
    
    for sample, file_set in files.items():
        cmd = "#!/bin/bash"
        if par != "":  # when the user specified a partition
            cmd += "\n#SBATCH -p " + par  # Otherwise, the "-p" will not be added into the command line, making the SLURM to use the default partition.
        cmd += "\n#SBATCH --job-name=" + sample  # keep the job name short to make it easier to debug
        cmd += "\n#SBATCH --ntasks=1"
        cmd += "\n#SBATCH --nodes=1"
        cmd += "\n#SBATCH --ntasks-per-node=1"
        cmd += "\n#SBATCH --cpus-per-task=1"
        cmd += "\n#SBATCH --mem-per-cpu=" + mem
        cmd += "\n#SBATCH --time=" + walltime
        if BLAST_MODULE != "":
            cmd += "\nmodule load " + BLAST_MODULE
        cmd += "\n{Bandage} distance {graph_file} {query_file} {other_options} > {result}\n".format(
        Bandage = bandage,
        graph_file = file_set.graph,
        query_file = file_set.query,
        other_options = other_args,
        result = os.path.join(outdir, sample + "__" + sf_out + ".tsv"))
        
        slurm_filename = sample + "__" + sf_out + ".slurm"  # put the sf_out into the filename to avoid name clash
        with open(slurm_filename, "w") as slurm_script:
            slurm_script.write(cmd)
            
        if run:
            os.system("sbatch " + slurm_filename)  # submit this job
            if (to_delay):
                time.sleep(delay)  # It is recommended to take a break before submitting the next job to give the SLURM enough time for sheduling resources.
                
    return


if __name__ == "__main__":
    main()
