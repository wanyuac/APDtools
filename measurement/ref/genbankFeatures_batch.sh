#!/bin/bash
#
# This script coverts GenBank files into FAA (protein sequence) or FFn (nucleotide sequence) files.
#
# Usage: genbankFeatures.sh --outdir=[output directory] --script_dir=[directory containing the script genbankFeatures.py] \
#                                     --protein --keep_pseudo --ext=[filename extension for GenBank files] \
#                                     --out_ext=[filename extension of output files] [input GenBank files]
# For example: bash genbankFeatures.sh ./references/gbk/to_convert ./references/faa protein
#
# It is derived from my script gbk2faa.sh.
# Copyright (C) 2016 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public License (GPL) version 3
# First edition: 27 Jun 2016; the latest edition: 3 Nov 2018

# default arguments
outdir="."
script_dir="."
nucl=true  # By default, the script pulls out nucleotide sequences.
skip_pseudo=true  # By default, the script skips pseudo genes.
gbk=()
ext=".gbk"
out_ext=".fasta"  # filename extension for output files

# read options
for i in "$@"; do
    case $i in
        --outdir=*)
        outdir="${i#*=}"  # output directory without the forward slash
        ;;
        --script_dir=*)
        script_dir="${i#*=}"
        ;;
		--protein)
		nucl=false
		;;
		--keep_pseudo)
		skip_pseudo=false
		;;
		--ext=*)
		ext="${i#*=}"
		;;
        --out_ext=*)
        out_ext="${i#*=}"
        ;;
        *)
        gbk+=(${i})  # Other arguments are treated as GenBank filenames.
    esac
done

# loop through the "ls" output
for i in "${gbk[@]}"; do
	name=$(basename $i ${ext})  # remove the path as well as the file name extension
	if $nucl; then  # for string comparison
		if $skip_pseudo; then
			python ${script_dir}/genbankFeatures.py --genbank $i --output ${outdir}/${name}.${out_ext} --features CDS --type nucleotide --ext $out_ext --skippseudo &
		else
			python ${script_dir}/genbankFeatures.py --genbank $i --output ${outdir}/${name}.${out_ext} --features CDS --type nucleotide --ext $out_ext --markpseudo &
		fi
	else
		if $skip_pseudo; then
			python ${script_dir}/genbankFeatures.py --genbank $i --output ${outdir}/${name}.${out_ext} --features CDS --type protein --ext $out_ext --skippseudo &
		else
			python ${script_dir}/genbankFeatures.py --genbank $i --output ${outdir}/${name}.${out_ext} --features CDS --type protein --ext $out_ext --markpseudo &
		fi
	fi
done
