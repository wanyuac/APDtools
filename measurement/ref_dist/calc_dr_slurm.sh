#!/bin/bash
# This bash script calls "calc_dr.R" to calculate reference distances for a list of strain names.
#
# This script assumes all input filenames ended with the extension "tsv".

# usage:
#   cat strain_names.txt | calc_dr.sh
#	ls -1 Feature_table/*.tsv | xargs -I {} basename {} ".tsv" | calc_dr.sh
#
# Copyright (C) 2017 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public License, version 3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
# The first edition: 12 Dec 2017; the latest edition: 9 Apr 2018

# defaults
in_dir="."
script="./calc_dr.R"
partition="sysgen,sysgen-long,main"
R_module=""
run=true
strains=""  # By default, the script extracts strain names from filenames.

# Read arguments
for i in "$@"; do  # loop through every argument
    case $i in
        --in_dir=*)  # Should not be ended with a forward slash.
		in_dir="${i#*=}"
		;;
        --strains=*)  # comma-delimited
        strains="${i#*=}"
        ;;
		--contig_sizes=*)
		contig_sizes="${i#*=}"
		;;
		--script=*)  # path to call calc_dr.R
		script="${i#*=}"
		;;
		--partition=*)
		partition="${i#*=}"
		;;
		--R_module=*)  # optional
		R_module="${i#*=}"
		;;
		--debug*)
		run=false
		;;
		*)  # Do nothing otherwise.
		;;
    esac
done

# Extract strain names
strain_names=()
if [ -z "$strains" ]; then
    IFS=$'\n'; for f in `ls -1 ${in_dir}/*.tsv`; do  # by default, extract strain names from all filenames
	    bn=`basename $f '.tsv'`
        strain_names+=("$bn")
    done
else
    strain_names=( `echo $strains | tr "," "\n"` )
fi

# Measure shortest-path distances for each set of query sequences
for s in "${strain_names[@]}"; do
	printf "\nMaking SLURM scripts for the strain ${s}:\n\n"
	cmd="#!/bin/bash\n\n"  # The echo command does not interpret the newline character.
	cmd+="#SBATCH -p ${partition}\n"
	cmd+="#SBATCH --job-name='calc dr: ${s}'\n"
	cmd+="#SBATCH --ntasks=1\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=1\n"
	cmd+="#SBATCH --cpus-per-task=1\n"
	cmd+="#SBATCH --mem-per-cpu=4096\n"
	cmd+="#SBATCH --time=0-24:0:00\n\n"
	if [[ !  -z  $R_module  ]]; then
		cmd+="module load ${R_module}\n\n"
	fi
	cmd+="cd ${PWD}\n\n"
	cmd+="Rscript ${script} --strain $s --features ${in_dir}/${s}.tsv --contig_sizes ${contig_sizes} --nofuzzy_start\n"
	echo -e $cmd
	if [ "$run" = true ]; then
		echo -e $cmd | sbatch 
	fi
	sleep 1
done

echo -e "\nAll jobs have been submitted successfully."
