#!/usr/bin/env python

"""
This script compiles distance measurements produced by Bandage. It expects input file names to observe the format:
[sample name][delimiter][source][suffix such as .tsv]. For instance, a valid input file may have the name:
NJST258__graph.tsv, in which a double underscores is the delimiter within the file name.

Usage example:
    python compile_dists.py -m *__*_spd.tsv -sf _spd.tsv -d "__" -od compiled -os shortest_path_distances.tsv \
    -of failed_query_paths.tsv -oe bandage_errors.tsv > compile_dists.log

Outputs:
    1. A tab-delimited table, e.g., "shortest_path_distances.tsv", about concatenated distance measurements.
    2. A tab-delimited table, e.g., "failed_query_paths.tsv" stored under the current working directory.

Python version 2 and 3 compatible
Copyright 2017 Yu Wan <wanyuac@gmail.com, https://github.com/wanyuac>
Licensed under the Apache License, Version 2.0
Development history: 21-22/11/2016, 22/5/2017; the latest edition: 14/12/2017
"""

from __future__ import print_function
from argparse import ArgumentParser
from collections import namedtuple
import os
import sys

# Constants
DEFAULT_OS = "distances.tsv"
DEFAULT_OF = "failed_query_paths.tsv"
DEFAULT_OE = "errors.tsv"
COLUMNS = ["query1", "query2", "sample", "distance", "node_number", "source", "query1_path", "query2_path", "orientation", "distance_path"]

def parse_arguments():
    # six arguments
    parser = ArgumentParser(description = "Compile distance measurements from Bandage.")
    parser.add_argument("-m", "--measurements", dest = "measurements", nargs = "+", type = str, required = True, \
                        help = "Output files of Bandage in the tab-delimited format.")
    parser.add_argument("-sf", "--suffix", dest = "suffix", type = str, required = False, default = ".tsv", \
                        help = "The suffix you need to drop from input file names before extracting sample and source information.")
    parser.add_argument("-d", "--delimiter", dest = "delimiter", type = str, required = False, default = "__", \
                        help = "The delimiter between the sample name and the suffix in each input file name")
    parser.add_argument("-od", "--output_dir", dest = "output_dir", type = str, required = False, default = ".", \
                        help = "Output directory")
    parser.add_argument("-os", "--output_dists", dest = "output_dists", type = str, required = False, default = DEFAULT_OS, \
                        help = "File name for successful distance measurements")
    parser.add_argument("-of", "--output_failures", dest = "output_failures", type = str, required = False, default = DEFAULT_OF, \
                        help = "File name for failed query paths")
    parser.add_argument("-oe", "--output_errors", dest = "output_errors", type = str, required = False, default = DEFAULT_OE, \
                        help = "File name for Bandage errors")
    return parser.parse_args()

def main():
    args = parse_arguments()
    files = []  # a list of Bandage_output objects
    for input_file in args.measurements:
        files.append(Bandage_output(input_file, args.suffix, args.delimiter))  # creates Bandage_output objects
    
    # initialise three kinds of output files
    if not os.path.exists(args.output_dir):
        os.system("mkdir " + args.output_dir)
    output_dists = initialise_output_file(os.path.join(args.output_dir, args.output_dists), DEFAULT_OS, COLUMNS)
    output_failures = initialise_output_file(os.path.join(args.output_dir, args.output_failures), DEFAULT_OF, ["sample", "query"])
    output_errors = initialise_output_file(os.path.join(args.output_dir, args.output_errors), DEFAULT_OE, ["sample", "error"])
    
    # process every object
    for record in files:  # each record is a reference to a Bandage_output object
        record.parse_record()  # parse Bandage's output files
        record.print_attr(output_dists, output_failures, output_errors)
    
    print("Completion: " + str(len(files)) + " files were combined successfully.")
    return

def initialise_output_file(output_file, default, column_names):
    try:
        f = open(output_file, "w")  # erase existing contents
        file_name = output_file
    except OSError:
        print("The name of the output file " + output_file + " is not accessible.\nReset the output to the default file, " + \
              default + " under the current working directory.")
        f = open(default, "w")
        file_name = default
    print("\t".join(column_names), file = f)  # write the header line into the file
    f.close()
    return file_name

class Bandage_output:
    'Stores and parses information in a single output file from Bandage.'
    
    # private class-level constants
    __NO_QUERY_PATHS = "No query paths found for "
    __FAILURE_HEADER_LEN = len(__NO_QUERY_PATHS)
    __BANDAGE_ERROR = "Bandage error: "
    __BANDAGE_ERROR_LEN = len(__BANDAGE_ERROR)
    __DISTANCE_HEADER = "Query 1 name"  # define the first line of the effective information; skip contents before this line.
    
    def __init__(self, result_file, sf, delim):
        # configures basic attributes of an object
        self.file = result_file  # about which result file is linked to this object
        self.sample, self.source = os.path.basename(self.file).rstrip(sf).split(delim)
        self.dists = []
        self.failures = []
        self.error = ""
        return
    
    def parse_record(self):  # parse Bandage's output file
        with open(self.file, "rU") as record:
            lines = record.read().splitlines()  # import all lines from the distance file
        i = -1  # index of elements in the list "lines"
        for line in lines:  # skip the header line in the distance file
            i += 1
            if line.startswith(Bandage_output.__NO_QUERY_PATHS):
                self.failures.append(line[Bandage_output.__FAILURE_HEADER_LEN : ])  # extract the allele ID from the current line
            elif line.startswith(Bandage_output.__BANDAGE_ERROR):
                self.error = line[Bandage_output.__BANDAGE_ERROR_LEN : ]  # get the error message; Only a single error is present becasue Bandage terminates here.
            elif line.startswith(Bandage_output.__DISTANCE_HEADER) and len(lines) - i > 1:
                self.dists = self.__parse_distance_field(lines[i + 1 : ])  # extract distances from the distance field
                break
            else:
                print("Warning: the %i-th line of the file %s is unrecognisable or there is no distance available." % \
                      (i + 1, self.file))
        return
    
    def print_attr(self, file_dists = None, file_failures = None, file_errors = None):
        # print attributes into files
        if file_dists != None and self.dists != []:
            with open(file_dists, "a") as f_dists:
                for line in self.dists:
                    print("\t".join([line.query1, line.query2, line.sample, line.distance, line.node_number, line.source, \
                                     line.query1_path, line.query2_path, line.orientation, line.distance_path]), \
                          file = f_dists)
                
        if file_failures != None and self.failures != []:
            with open(file_failures, "a") as f_failures:
                for allele in self.failures:
                    print("\t".join([self.sample, allele]), file = f_failures)
        
        if file_errors != None and self.error != "":
            with open(file_errors, "a") as f_errors:
                print("\t".join([self.sample, self.error]), file = f_errors)
            
        return

    def __parse_distance_field(self, distance_field):  # a private method
        dists = []  # a list of named tuples; each named tuple stores information of a single line
        Line = namedtuple("Line", COLUMNS)
        for line in distance_field:  # go through every row (line) of the distance table
            fields = line.split("\t")
            dists.append(Line(query1 = fields[0],
                              query2 = fields[1],
                              sample = self.sample,
                              distance = fields[5],
                              node_number = str(self.__count_nodes(fields[6])),
                              source = self.source,
                              query1_path = fields[2],
                              query2_path = fields[3],
                              orientation = fields[4],
                              distance_path = fields[6]))
        return dists
    
    def __count_nodes(self, dist_path):  # counts how many nodes are there in a path
        return dist_path.count("+") + dist_path.count("-")
    
if __name__ == "__main__":
    main()
