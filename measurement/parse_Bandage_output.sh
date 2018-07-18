#!/bin/bash
#
# Copyright (C) 2017 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public License, version 3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
# First edition and the latest edition: 14 Dec 2017

# Help information ###############
display_usage() {
    echo "This script utilises compile_dists.py to split each output file of the Bandage distance command into
    three files (distances, failures and errors).
    Usage:
    ./parse_Bandage_output.sh --outdir=[output directory] --compiler=[path to compile_dists.py] --suffix='.tsv' [input directory]/*.tsv
    This script assumes every input file names follow the same format: [strain name]__[source][filename extension, such as '.tsv'].
    "
}

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	display_usage
	exit
fi

# Defaults ###############
outdir=$PWD
suffix=".tsv"

# get the directory of this script itself
script_dir=$(cd `dirname $0` && pwd)  # where compile_dists.py is expected to be installed
compiler="${script_dir}/compile_dists.py"

# Read arguments ###############
raw_files=()
for i in "$@"; do  # loop through every argument
    case $i in
        --outdir=*)
        outdir="${i#*=}"
        ;;
		--compiler=*)
		compiler="${i#*=}"
		;;
        --suffix=*)
        suffix="${i#*=}"  # filename suffix to be chopped off from every input file name
        ;;
        *)
        raw_files+=( "${i#*=}" )  # "0" or "1"  # load every file name
        ;;
    esac
done

# Run compile_dists.py for every strain ###############
if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi

for f in "${raw_files[@]}"; do
    base_name=`basename $f $suffix`  # e.g., ./distances/strain1__graph.tsv => strain1__graph
    base_name_fields=(`echo $base_name | sed "s/__/\n/g"`)
    strain=`echo ${base_name_fields[0]}`
    source=`echo ${base_name_fields[1]}`
    echo "Processing the Bandage output for the strain ${strain}."
    python $compiler --measurements $f --suffix "$suffix" --delimiter "__" --output_dir $outdir --output_dists ${strain}__dists.tsv --output_failures ${strain}__missingQueryPaths.tsv --output_errors ${strain}__errors.tsv > ${strain}__compileDists.log
done

echo "All of ${#raw_files[@]} files have been parsed successfully."
