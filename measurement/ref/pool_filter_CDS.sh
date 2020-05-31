#!/bin/bash
# Pool and filter CDS of every strain.
#
# Example
#    pool_filter_CDS.sh --dir_in="./CDS" --dir_loci="./Loci" --dir_out="./CDS/Subset" [strain1] [strain2] ...
# Copyright (C) 2017 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public License, version 3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
# The first edition: 12 Dec 2017; the latest edition: 10 Apr 2018

# Define defaults
dir_in="."  # Directory of input FASTA files.
dir_out="./Subset"
dir_script="."  # where extract_fasta_loci.py stores
strains=()
suffix="__loci.txt"

# Read arguments
for i in "$@"; do  # loop through every argument (e.g. "--gbk=*.gbk" is treated as a single word by the $@ operator)
    case $i in
        --dir_in=*)
        dir_in="${i#*=}"
        ;;
        --dir_loci=*)
        dir_loci="${i#*=}"
        ;;
        --dir_out=*)
        dir_out="${i#*=}"  # output directory without the forward slash
        ;;
        --dir_script=*)
        dir_script="${i#*=}"
        ;;
        --suffix=*)
        suffix="${i#*=}"
        ;;
        *)
        strains+=(${i})  # Other arguments are treated as GenBank filenames.
        ;;
    esac
done

# Make output directories
if [ ! -d $dir_out ]; then
    mkdir -p $dir_out
fi

# Pool FASTA files of each strain into a single multi-FASTA file
for s in "${strains[@]}"; do
    cat ${dir_in}/${s}__*.fna | python ${dir_script}/extract_fasta_loci.py "$(cat ${dir_loci}/${s}${suffix})" > ${dir_out}/${s}.fna
done
