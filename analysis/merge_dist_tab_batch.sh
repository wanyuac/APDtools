#!/bin/bash
#
# Submit SLURM jobs of merge_dist_tab.R for multiple strains. Strain names are comma-delimited
# and provided as a single string for the first (and the only) argument. All jobs will be running
# under the current working directory by default.
# Copyright 2017 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# First edition: 17-18 Dec 2017; the latest edition: 30 Mar 2018

# Split $1 into an array or display help information ###############
display_usage() {
    echo "Submitting SLURM jobs of merge_dist_tab.R for multiple strains.
    Usage:
        ./merge_dist_tab_batch.sh --strains='strain1,strain2,strain3' --dr='/vlsci/SG0006/shared/wan/Ec/Dist/Sim/Dr' --d='/scratch/sysgen/wan/Ds/Parsed/Contig' --cores=8
	"
}

# Defaults ###############
dr_suffix="__dr.tsv"  # filename suffix for reference distances
d_suffix="__dists.tsv"  # filename suffix for parsed Bandage outputs
partition="sysgen-long,sysgen,main"
n_cores=8
mem_per_cpu=8192
walltime="6-0:0:00"
max_d=300000
max_nodes=10
qs=""
sampling_rate_d=1  # That is no sampling by default.
sampling_rate_dr=1
outdir=$PWD  # output directory
R_module="R/3.3.3-vlsci_intel-2015.08.25"  # replace it to your own module or comment it out
script_dir=$(cd `dirname $0` && pwd)  # assuming merge_dist_tab.R is stored under the same directory as this script
n=0  # strain count
run=true

# Read arguments
for i in "$@"; do  # loop through every argument
    case $i in
        --strains=*)
        strains="${i#*=}"
        strains=( `echo $strains | tr "," "\n"` )
        n=${#strains[@]}  # number of strains
        ;;
		--dr=*)
		dir_dr="${i#*=}"
		;;
        --d=*)
        dir_d="${i#*=}"  # directory of distances measured using Bandage
        ;;
        --cores=*)
        n_core="${i#*=}"
        ;;
        --mem_per_cpu=*)
        mem_per_cpu="${i#*=}"
        ;;
		--sampling_rate_d=*)
        sampling_rate_d="${i#*=}"
        ;;
        --sampling_rate_dr=*)
        sampling_rate_dr="${i#*=}"
        ;;
        --max_nodes=*)
        max_nodes="${i#*=}"
        ;;
        --max_d=*)
        max_d="${i#*=}"
        ;;
        --qs=*)
        qs="${i#*=}"
        ;;
        --walltime=*)
        walltime="${i#*=}"
        ;;
        --dr_suffix=*)
        dr_suffix="${i#*=}"
        ;;
        --d_suffix=*)
        d_suffix="${i#*=}"
        ;;
        --partition=*)
        partition="${i#*=}"
        ;;
        --debug*)
		run=false
		;;
        *)  # Do nothing otherwise.
        ;;
    esac
done

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	display_usage
	exit
fi

# Submitting jobs through SLURM ###############
for s in "${strains[@]}"; do
    cmd='#!/bin/bash\n'  # Use '' instead of "", otherwise this line will be executed immediately.
    cmd+='\n#SBATCH --partition='$partition  # Similarly, bash executes this command immediately rather than adding it to cmd when double quotes are used.
    cmd+='\n#SBATCH --job-name='$s
    cmd+='\n#SBATCH --ntasks=1'
    cmd+='\n#SBATCH --cpus-per-task='$n_core
    cmd+='\n#SBATCH --mem-per-cpu='$mem_per_cpu
    cmd+='\n#SBATCH --time='$walltime
    cmd+="\n\nmodule load ${R_module}"
    cmd+="\ncd ${outdir}\n\n"
    cmd+="Rscript ${script_dir}/merge_dist_tab.R --dr ${dir_dr}/${s}${dr_suffix} --d ${dir_d}/${s}${d_suffix} "
    cmd+="--max_d ${max_d} --max_nodes ${max_nodes} "
    if [ ! -z "$qs" ]; then
        cmd+="--qs ${qs} "
    fi
    cmd+="--sampling_rate_d ${sampling_rate_d} --sampling_rate_dr ${sampling_rate_dr} --cores ${n_core} --output ${s}__Dm.tsv\n"
    echo -e $cmd  # display the commands
    if [ "$run" = true ]; then
        echo -e $cmd | sbatch  # actually execute the SLURM script
    fi
    sleep 1
done
