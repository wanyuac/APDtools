#!/bin/bash
# Tabulate gene features for every GenBank file and pull them together afterwards
# Recommend to run this bash script under an interactive session through the screen command when there are too many files to process
#
# Command:
#   ./gbk2table.sh --script=[script that extracts features] --outdir=[output directory without the forward slash] [GenBank files]
# Name format of GenBank files:
#   [strain name]__[accession].[filename extension]
#   e.g. 2011C-3493__CP003290.gbk for the chromosomal sequence of the Escherichia coli strain 2011C-3493
#
# Copyright (C) 2017-2018 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public License, version 3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
# First edition and the latest edition: 12 Dec 2017

# Define default or initial values
script="/vlsci/SG0006/shared/code/holtlab/genbank2table.py"  # the default script for feature extraction
outdir=$PWD
gbk=()
strains=()  # an empty arrary
default_header="contig_id\tcontig_name\tfeature_type\tstart\tend\tlength\tstrand\tlocus_tag\tproduct\tpseudo"  # header line from the genbank2table.py script
new_header="Accession\tContig_name\tFeature_type\tStart\tEnd\tLength\tStrand\tLocus_tag\tProduct\tPseudo"  # header of final TSV files

# Read arguments from the command line
for i in "$@"; do  # loop through every argument (e.g. "--gbk=*.gbk" is treated as a single word by the $@ operator)
    case $i in
		--script=*)
		script="${i#*=}"  # script that extracts features
		;;
        --outdir=*)
        outdir="${i#*=}"  # output directory without the forward slash
        ;;
        *)
        gbk+=(${i})  # Other arguments are treated as GenBank filenames.
    esac
done

# Check the output directory
if [ ! -d $outdir ]; then
    mkdir $outdir
fi
cd $outdir

# Process every strain when there are GenBank files found
if [ ${#gbk[@]} -gt 0 ]; then
    # Extract feature information from the GenBank file
    for f in "${gbk[@]}"; do
        filename=`basename ${f}`
        base_name="${filename%.*}"  # base of the filename
        base_name_fields=(`echo $base_name | sed "s/__/\n/g"`)  # read the split string as an arrary. Do not use the tr command as it allows a partial match "_" of "__".
        strain_name=`echo ${base_name_fields[0]}`
        strains+=(${strain_name})  # push this strain name into the arrary "strains"
        python $script --genbank $f --features CDS --qualifiers 'locus_tag,product,pseudo' > ${outdir}/${base_name}.txt  # each with a header line
    done

    # Concatenate feature files of each strain into a single one
    # deduplicate strain names in the array strains
    strains=(`echo "${strains[@]}" | tr " " "\n" | sort -u | tr "\n" " "`)  # Reference: https://stackoverflow.com/questions/13648410/how-can-i-get-unique-values-from-an-array-in-bash
    for s in "${strains[@]}"; do
        printf "${new_header}\n" > ${outdir}/${s}.tsv  # make the header line; change the filename extension to make it easier to sort
        cat ${outdir}/${s}__*.txt | grep -P -v "$default_header" >> ${outdir}/${s}.tsv  # skip the header lines in each file
    done
else
    echo "This script is stopped as there is no GenBank files found."
fi
