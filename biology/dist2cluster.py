#!/usr/bin/env python

"""
Measuring the shortest distance between target alleles and a cluster of other alleles in every genome assembly.

Algorithm:
1. Read a cluster-definition file to acquire cluster IDs, alleles per cluster and strains harbouring each cluster.
2. Import a sequence database of alleles that are members of clusters.
3. Import target allele sequences.
4. Create a dictionary for assembly graphs.
5. For each cluster:
    1) write a FASTA file into the output directory to include relevant cluster allele sequences and all target sequences
    2) run Bandage to measure the shortest-path physical distances
    
Output directories
    outdir/
        queries/
        raw/

Examples:
python dist2cluster.py -t intI.fna -c clusters.tsv -d cluster_alleles.fna -a assemblies/*__spades.fastg \
-p query_paths.tsv -o output/intI_dists.tsv -s "__spades" -k

Notes:
    1. This algorithm runs Bandage to measure the shortest-path distance between best query paths (hence the option
    --allquerypaths is not used here). This is the same as the script dist_from_graphs.py.
    2. It does not work under different alignment filter for targets and cluster members so far.
    3. Copies of the same allele in a genome may cause errors in the distance measurements.

Dependencies: Python version 2 and 3 (recommended) compatible, BioPython, BLAST and pandas (http://pandas.pydata.org/)

Copyright 2017 Yu Wan <wanyuac@gmail.com>
Licensed under the Apache License, Version 2.0
Created on 10 Oct 2017, the lastest edition: 18 Oct 2017
"""

from __future__ import print_function
import os
import sys
import shlex  # parses a command line into a list
import pandas as pd  # requires users to download and install on their computers
import multiprocessing as mp
from argparse import ArgumentParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from distutils import spawn
from collections import namedtuple, defaultdict
from subprocess import check_call, STDOUT, PIPE  # to launch concurrent processes of Bandage
from multiprocessing.pool import ThreadPool  # to limit the number of the concurrent processes


def parse_arguments():
    parser = ArgumentParser(description = "Measuring distances from target alleles to allele clusters in genome assemblies.")
    parser.add_argument("-t", "--targets", dest = "targets", type = str, required = True, \
                        help = "A FASTA file for sequences of target alleles")
    parser.add_argument("-c", "--clusters", dest = "clusters", type = str, required = True, \
                        help = "A tab-delimited cluster-definition file")
    parser.add_argument("-d", "--db", dest = "db", type = str, required = True, \
                        help = "A FASTA file containing all allele sequences listed in the cluster-definition file as a database.")
    parser.add_argument("-a", "--assemblies", dest = "assemblies", nargs = "+", type = str, required = True, \
                        help = "A list of assembly graphs in FASTA/FASTG/GFA format")
    parser.add_argument("-p", "--paths", dest = "paths", type = str, required = True, help = "A tab-delimited file for best query paths of alleles in each assembly")
    parser.add_argument("-o", "--outdir", dest = "outdir", type = str, required = False, default = "output", \
                        help = "Output directory (without the ending backslash)")
    parser.add_argument("-s", "--suffix", dest = "suffix", type = str, required = True, \
                        help = "Suffix with which strain names are extracted from file names of genome assemblies.")
    parser.add_argument("-k", "--skip", dest = "skip", action = "store_true", \
                        help = "Set it to avoid overwriting existing output files.")
    parser.add_argument("-b", "--bandage", dest = "bandage", type = str, required = False, default = "Bandage", \
                        help = "Path to Bandage")
    parser.add_argument("-u", "--threads", dest = "threads", type = int, required = False, default = 1, \
                        help = "Number of threads to run Bandage simultanesously (default: single thread)")
    parser.add_argument("-m", "--compiler", dest = "compiler", type = str, required = False, default = "compile_dists.py", \
                        help = "Path to compile_dists.py")
    parser.add_argument("-n", "--max_node_num", dest = "max_node_num", type = int, required = False, default = 0, \
                        help = "Maximal number of nodes for reliable measurements of physical distances. Set it to zero to disable this filter.")
    parser.add_argument("-l", "--max_dist", dest = "max_dist", type = int, required = False, default = 0, \
                        help = "Maximal physical distance to be considered as reliable. Set it to zero to disable this filter.")
    parser.add_argument("-r", "--other_args", dest = "other_args", type = str, required = False, default = "--ifilter 95 --minpatcov 0.5 --evfilter 1e-3", \
                        help = "Other arguments directly passed to Bandage.")
    return parser.parse_args()


def main():
    args = parse_arguments()
    
    # Start up
    check_for_blast()
    check_file_exists([args.targets, args.clusters, args.db, args.bandage] + args.assemblies)
    if args.threads > mp.cpu_count():  # check the number of threads
        print("Warning: %i threads are specified but there are only %i CPUs available." % (args.threads, mp.cpu_count()))
        sys.exit("Please check your thread specification and run this script again.")
    elif args.threads <= 0:
        sys.exit("Argument error: the number of threads must be a positive integer.")
        
    outdir = args.outdir
    if outdir == "":  # when the user only supplies a file name without a relevant path
        outdir = "."
    else:
        check_dir_exists(outdir)  # set up an output directory
    outdir_raw = check_dir_exists(os.path.join(outdir, "raw"))  # the directory for raw outputs from Bandage
    outdir_queries = check_dir_exists(os.path.join(outdir, "queries"))  # outdir/queries, for individual FASTA files of query sequences
        
    # Import and parse the cluster definitions
    cls = import_cluster_def(args.clusters)  # cls is a dictionary of namedtuples {cluster id : [alleles, strains]}
    
    # Import the allele database for clusters
    db = read_fasta(args.db)  # {allele : sequence}
    
    # Import target sequences
    targets = read_fasta(args.targets)  # {allele : sequence}; use list(targets.keys()) to get target allele names
    target_sizes = get_allele_sizes(targets)  # a dictionary of base-pair counts
    
    # Create a dictionary of assembly graphs
    assemblies = store_assembly_names(args.assemblies, args.suffix)
    
    # Import best and unique query paths
    query_paths = import_query_paths(args.paths)
    
    """
    Use multiprocessing to run Bandage and catch raw outputs from each subprocess
    In total, there are n(cluster) * n(strains) Bandage processes to be submitted.
    It is recommended to use a smaller number (<=16) of cores to launch jobs so as
    to avoid the speed bottleneck of the hard driver in writing data.
    """
    print("Using " + str(args.threads) + " cores to run Bandage jobs.")
    Cmd = namedtuple("Cmd", ["cmd", "out"])  # element of a list for parallel computing
    commands = ()
    
    # Compose a command line for each strain
    raw_outputs = []
    for cid, cl in cls.items():
        query_fasta = write_queries(cid, cl, db, targets, outdir_queries, args.skip)  # write a FASTA file of query sequences for each cluster
        for s in cl.strains:
            bandage_out = os.path.join(outdir_raw, "__". join([s, cid, "dists.tsv"]))
            cmd_fields = [args.bandage, "distance", assemblies[s], query_fasta] + shlex.split(args.other_args)
            if os.path.isfile(bandage_out) and args.skip:
                print("Skipped: not rewrite the result file " + bandage_out)
            else:
                print("Command to run: " + " ".join(cmd_fields))
                commands += (Cmd(cmd = cmd_fields, out = bandage_out), )  # append an element into the tuple
            raw_outputs.append(bandage_out)
    
    # Run Bandage jobs concurrently
    if len(commands) > 0:  # else, skip this step
        tp = ThreadPool(processes = args.threads)  # create a ThreadPool object that runs core_num worker threads
        tp.map_async(run_cmd, commands)  # The results object belongs to the class multiprocessing.pool.MapResult.
        tp.close()
        tp.join()
    
    # Compile Bandage outputs into a single file (dists.tsv under the outdir)
    check_call(args = ["python", args.compiler, "-sf", "__dists.tsv", "-d", "__", "-od", outdir, "-m"] + raw_outputs, \
               shell = False, stdout = PIPE, stderr = STDOUT)  # Cluster ID (cid) in each input file name serves as the "source" for the script compile_dists.py.
    
    # Preprocess the merged distance table
    ds = pd.read_csv(os.path.join(outdir, "dists.tsv"), sep = "\t", \
                     names = ["query1", "query2", "strain", "distance", "node_number", "cluster", "query1_path", "query2_path", "orientation", "distance_path"], \
                     header = 0)  # import the compiled distance measurements as a data frame
    ds = ds[["cluster", "strain", "query1", "query2", "query1_path", "query2_path", "distance", "node_number", "distance_path"]]  # rearrange columns through a list of column names and drop the orientation column
    
    # filter for measurement reliability
    if args.max_node_num > 0:
        ds = ds.loc[ds["node_number"] <= args.max_node_num]
    is_empty(ds)
    if args.max_dist > 0:
        ds = ds.loc[ds["distance"] <= args.max_dist]
    is_empty(ds)
    
    # filter for target-cluster distances
    target_alleles = list(targets.keys())  # names of target alleles
    ds = ds[~(ds["query1"].isin(target_alleles) & ds["query2"].isin(target_alleles))]  # remove distances between target alleles
    ds = reorder_queries(ds, target_alleles)  # to make target alleles alway appear in the query1 column so further programming becomes easier
    
    """
    Remove distances between cluster alleles. Since target alleles are only present in the query1 column after reorder_queries, any non-target alleles
    in the query1 column indicate distances to other non-target alleles. After this filter, the query1 column only contains names of the target alleles.
    """
    ds = ds[ds["query1"].isin(target_alleles)]
    is_empty(ds)  # move on if there are rows left
    ds.columns = ["cluster", "strain", "target", "member", "target_path", "member_path", "distance", "node_number", "distance_path"]  # rename columns
    ds = ds.reset_index(drop = True)  # reset indices and do not include indices as a new column. Otherwise, lock_member_alleles raises an error of list index overflow.
    ds = lock_member_alleles(ds, query_paths)  # Pinpoint best query paths of cluster-member alleles (that is, removing distances involving other query paths)
    
    # Make a strain-level summary for every cluster
    report_str = strain_level_summary(ds)
    
    # Finally, save results as text files
    ds.to_csv(os.path.join(outdir, "dists_filtered.tsv"), sep = "\t", encoding = "utf-8", index = False)
    report_str.to_csv(os.path.join(outdir, "dists_summary.tsv"), sep = "\t", encoding = "utf-8", index = False)
    target_sizes.to_csv(os.path.join(outdir, "target_sizes.tsv"), sep = "\t", encoding = "utf-8", index = False)
    
    return


def strain_level_summary(ds):
    strains = set(ds["strain"])  # de-duplicate the list of strain names
    clusters = set(ds["cluster"])
    
    # Make a data frame from a list of dictionaries
    report = []
    for c in clusters:
        for s in strains:
            ds_sub = ds.loc[(ds["cluster"] == c) & (ds["strain"] == s)]
            target_alleles = set(ds_sub["target"])  # Usually there is only a single target.
            for a in target_alleles:
                ds_sub1 = ds_sub.loc[ds_sub["target"] == a]
                ds_sub1 = ds_sub1.reset_index(drop = True)  # in order to use i_max and i_min below
                dists = list(ds_sub1["distance"])  # "Int64Index" object is not callable without being converted into a list.
                i_max = dists.index(max(dists))  # get the index of the first maximum distance
                i_min = dists.index(min(dists))  # the index of the first minimum distance
                member_min, mem_min_path, d_min, d_min_nodes = ds_sub1.loc[i_min, ["member", "member_path", "distance", "node_number"]]
                member_max, mem_max_path, d_max, d_max_nodes = ds_sub1.loc[i_max, ["member", "member_path", "distance", "node_number"]]
                
                """
                Push a dictionary into the list "report".
                Only take the first target path because it is the same in every row of ds_sub1.
                """
                report.append({"cluster" : c, "strain" : s, "target" : a, \
                               "member_min" : member_min, "member_max" : member_max, \
                               "d_min" : d_min, "d_max" : d_max, \
                               "d_count" : len(ds_sub1.index), \
                               "target_path" : ds_sub1.loc[0, "target_path"], \
                               "mem_min_path" : mem_min_path, "mem_max_path" : mem_max_path, \
                               "d_min_nodes" : d_min_nodes, "d_max_nodes" : d_max_nodes})
                out = pd.DataFrame(report)  # An error arises when using the command report = pd.DataFrame(report).
                out = out[["cluster", "strain", "target", "member_min", "member_max", "d_min", "d_max", \
                           "d_count", "target_path", "mem_min_path", "mem_max_path", "d_min_nodes", "d_max_nodes"]]
    return out


def lock_member_alleles(df, paths):  # df: pass-by-reference
    """
    Remove rows of df where member allele paths (member_path) are not found in a list of given paths.
    """
    keep = [False] * len(df.index)  # a logical list (with the length equalling the row count) for rows to keep
    for i, row in df.iterrows():
        try:
            keep[i] = row["member_path"] == paths[row["member"]][row["strain"]]
        except IndexError:
            print("Error: the key [%s][%s] is not found in the dictionary of paths." % (row["member"], row["strain"]))
            raise
    df = df[keep]
    return df


def reorder_queries(df, targets):  # df is a reference to a data frame (not pass-by-value)
    """
    Reorder query1 and query2 names and paths to put names and paths of the target alleles into the query1 and query1_path columns
    Notice Python passes a reference of df to this function instead of its value, hence any modifications to df modifies the data
    frame it points to. Use d = df.copy() if you do not want the function to directly modify your argument data frame, although this
    is not necessary for my purpose of swapping the query names. Using a reference is more memory-efficient than copying the value,
    which is beneficial when we have a large data frame to process.
    """
    for i, row in df.iterrows():
        q1, q2, q1_path, q2_path = row[["query1", "query2", "query1_path", "query2_path"]]
        if q2 in targets:
            # Because I already excluded rows where both q1 and q2 are in the targets list, q1 must not in targets here.
            df.loc[i, ["query1", "query2", "query1_path", "query2_path"]] = [q2, q1, q2_path, q1_path]
    return df


def run_cmd(cmd):
    try:
        bandage_out = open(cmd.out, "w")
        bandage_err = open(cmd.out + ".err", "w")
        check_call(args = cmd.cmd, shell = False, stdout = bandage_out, stderr = bandage_err)  # launch a Bandage session as a thread
        bandage_out.close()
        bandage_err.close()
    except IOError as e:
        sys.exit("I/O error on '%s': %s" % (e.filename, e.strerror))
    except CalledProcessError as e:
        sys.exit("The command line %s failed, returned code %d." % (" ".join(cmd.cmd), e.returncode))
    except OSError as e:
        sys.exit("Failed to run the shell: %s" % (str(e)))
    return
    

def import_query_paths(path_file):
    '''
    A path file is comprised of three columns separated by tab keys: query, sample, path. The file must not contain
    column names. Each path anchors the best hit of an allele in a given assembly, ruling out alternative hits when
    a low query coverage is allowed for target sequences.
    '''
    paths = defaultdict(dict)
    with open(path_file, "rU") as f:
        lines = f.read().splitlines()
    for line in lines:
        query, sample, query_path = line.split("\t")
        paths[query][sample] = query_path
    return paths


def write_queries(cid, cl, db, targets, outdir_queries, skip):
    """
    Create a FASTA file of query sequences for a given cluster.
    cid: cluster ID; cl: a Cluster object; db: an allele database; targets: target sequences;
    outdir_queries: output directory for FASTA files ([parental outdir]/queries)
    Output file: [parental outdir]/queries/[cid]__queries.fna
    This function returns the name of the output FASTA file.
    """
    out = os.path.join(outdir_queries, cid + "__queries.fna")
    if not (os.path.isfile(out) and skip):
        f = open(out, "w")
        
        # first, write target sequences
        for name, seq in targets.items():
            f.write(">" + name + "\n" + seq + "\n")
        
        # next, write sequences of cluster alleles
        for allele in cl.alleles:
            f.write(">" + allele + "\n" + db[allele] + "\n")
        
        f.close()
    return out
    

def import_cluster_def(filename):
    """
    Importing cluster definitions from a tab-delimited file of three columns.
    Returns a dictionary of namedtuples.
    """
    clusters = {}
    Cluster = namedtuple("Cluster", ["alleles", "strains"])
    with open(filename, "rU") as f:
        lines = f.read().splitlines()
    for line in lines:
        cid, a, s = line.split("\t")  # cluster ID, allele name and strain name
        clusters[cid] = Cluster(alleles = a.split(","), strains = s.split(","))
    return clusters


def read_fasta(filename):
    """ Returns a dictionary of sequences """
    db = {}
    records = list(SeqIO.parse(filename, "fasta"))
    for rec in records:
        db[rec.name] = str(rec.seq)  # Do not use rec.id as the key because it is the whole header line. rec.id = rec.name + rec.description.
    return db


def get_allele_sizes(d):
    """ Measure the length of each sequence in the sequence dictionary d """
    s = []
    for allele, seq in d.items():
        s.append({"Allele" : allele, "Size" : len(seq)})
    return pd.DataFrame(s)


def store_assembly_names(filenames, suffix):
    assemblies = {}
    for f in filenames:
        strain = os.path.basename(f)  # Otherwise, the common directory name goes into every key with assembly names.
        strain = strain.replace(suffix, "")  # chop off the suffix
        assemblies[strain] = f
    print("Links to %d assemblies have been imported from %d filenames." % (len(assemblies), len(filenames)))
    return assemblies


def check_for_blast():
    makeblastdb_path = spawn.find_executable("makeblastdb")
    blastn_path = spawn.find_executable("blastn")
    blast_installed = (makeblastdb_path != None and blastn_path != None)
    if not blast_installed:
        sys.exit("Error: could not find BLAST program")
    return


def check_file_exists(files):
    for filename in files:
        if not os.path.isfile(filename):
            sys.exit("Error: could not find " + filename)
    return


def check_dir_exists(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    return dirname


def is_empty(df):
    """ Exit if a data frame is empty """
    if df.size == 0:
        print("No further outputs are generated as there is no distance between target alleles and clusters available for any strains.")
        sys.exit(0)  # normal exit
    

if __name__ == "__main__":
    main()
