#!/bin/bash
# SPD comparisons
#
# This script compare SPDs measured in de novo assemblies to true SPDs (measured in circularised complete genomes)
# for every pair of query sequences (de-duplicated CDS).
#
# Run this script in a screen session
#
# Usage
#	bash pipeline_SPDs.sh [strain name] [stage] [comma-delimited graph names] [no. of iterations]
# 	bash pipeline_SPDs.sh HS11286 1 'c,l,k69' 15
#	But do not mix 'ctg' with others: bash pipeline_SPDs.sh HS11286 1 'ctg,l,c' 15 will not work correctly.
# Use the command 'sacct | less' to inspect status of every job. Hence the email notification is not necessary.
#
# Dependencies: SLURM (queuing system), Python, BioPython, R, nucleotide BLAST.
#
# Copyright 2018 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# First version: 2 July 2016; the latest edition: 4 Nov 2018


# ARGUMENT SETTINGS ###############
# Environmental constants
WALLTIME_SHORT="0-04:00:00"  # walltime for short jobs
WALLTIME_LONG="0-08:00:00"  # walltime for long jobs
MEM_SHORT=16384  # RAM (Mb) for short jobs
MEM_LONG=32768  # RAM (Mb) for long jobs

# Default arguments
out_dir="./SPD"
iter_num=5
query_suffix="fna"
graph_suffix="fasta"  # gfa, fna, etc.
query_num=1000
bandage_dir="./Bandage"
blast_module=""
python_module=""
biopython_module=""
R_module=""
cores_short=1
cores_long=4
par_short="main"
par_long="main"
strains=""
bandage_args="--evfilter 1e-5 --ifilter 95 --minpatcov 0.95"
out_suffix="graph"
account=""
n_max=10
d_max="300e3"

# Read 19 user specific arguments
for i in "$@"; do  # loop through every argument
    case $i in
        --out_dir=*)  # the working directory
        out_dir="${i#*=}"
        ;;
		--out_suffix=*)  # suffix for output sub-directories and SPD files
		out_suffix="${i#*=}"
		;;
		--script_dir=*)  # path to physDist
		script_dir="${i#*=}"
		;;
		--bandage_dir=*)  # path to call Bandage
		bandage_dir="${i#*=}"
		;;
		--bandage_args=*)  # arguments for Bandage to measure SPDs in assembly graphs and contigs
		bandage_args="${i#*=}"
		;;
		--blast_module=*)  # module name for nucleotide BLAST, required by Bandage
		blast_module="${i#*=}"
		;;
		--python_module=*)  # module name for Python
		python_module="${i#*=}"
		;;
        --biopython_module=*)  # BioPython
        biopython_module="${i#*=}"
        ;;
        --R_module=*)  # name for the R module
        R_module="${i#*=}"
        ;;
		--ref_dir=*)  # directory for reference sequences (*.gfa)
		ref_dir="${i#*=}"
		;;
		--query_dir=*)  # directory for query sequences
		query_dir="${i#*=}"
		;;
		--query_suffix=*)  # filename extension (no dot) for input query files (CDS extracted from GenBank files)
		query_suffix="${i#*=}"
		;;
		--query_num=*)  # number of query sequences to be sampled from the total in each iteration
		query_num="${i#*=}"
		;;
        --iter_num=*)  # number of iterations
		iter_num="${i#*=}"
		;;
		--graph_dir=*)  # contigs (*.fasta), assembly graphs (*.gfa), etc.
		graph_dir="${i#*=}"
		;;
		--graph_suffix=*)  # filename extension for input subject files, no dot sign
		graph_suffix="${i#*=}"
		;;
		--strains=*)  # comma-delimited strain names
		strains="${i#*=}"
		;;
		--cores_short=*)  # number of computational cores to be used for each short SLURM job
		cores_short="${i#*=}"
		;;
        --cores_long=*)  # number of computational cores to be used for each short SLURM job
		cores_long="${i#*=}"
		;;
		--par_short=*)  # SLURM partition for short jobs to join
		par_short="${i#*=}"
		;;
		--par_long=*)  # name of SLURM partition for long jobs to join
		par_long="${i#*=}"
		;;
        --account=*)  # project name to join a partition
		account="${i#*=}"
		;;
        --n_max=*)  # maximum node number of SPDs
        n_max="${i#*=}"
        ;;
        --d_max=*)  # maximum distance
        d_max="${i#*=}"
        ;;
        *)  # do nothing by default
        ;;
    esac
done

if [ ! -z "$strains" ]; then
	strains=(`echo $strains | tr "," "\n"`)  # parse strain names
else
	echo "Argument error: strain names must be specified."
	exit 1
fi
iterations=(`seq 1 $iter_num`)  # indices of iterations

# ENVIRONMENTAL CONFIGURATION ###############
if [ ! -z "$python_module" ]; then  # when modules are set
    load_python=true
    module load "$python_module"  # Load a module for the current bash environment, where the module cannot be seen from outside.
	# [debug] module list  # print module names
else
    load_python=false
fi

if [ ! -z "$biopython_module" ]; then
    module load "$biopython_module"  # for sample_seqs.py
fi

if [ ! -z "$R_module" ]; then
    load_R=true
    module load $R_module
else
    load_R=false
fi

if [ ! -z "$account" ]; then
    add_account=true
else
    add_account=false
fi

# Set up the output directory
if [ ! -d "$out_dir" ]; then
	echo "Making output directories ${out_dir}."
	mkdir $out_dir
    mkdir "${out_dir}/Merged_${out_suffix}"
	for s in ${strains[@]}; do
		mkdir ${out_dir}/${s}
		for i in ${iterations[@]}; do  # A double quote retrieves the value of the variable.
			mkdir -p ${out_dir}/${s}/${i}/Parsed
		done
	done
fi


# STAGE 1: SAMPLE FROM CDS PER STRAIN ###############
for s in ${strains[@]}; do
	checkpoint="${out_dir}/${s}/stage_1.success"
	if [ ! -f "$checkpoint" ]; then
		for i in ${iterations[@]}; do
			echo "Sampling ${query_num} queries for strain ${s}, round ${i}."
			python ${script_dir}/measurement/sample_seqs.py "${query_dir}/${s}.${query_suffix}" "$query_num" "${out_dir}/${s}/${i}"
		done
		touch "$checkpoint"
	else
		echo "CDS of strain ${s} have been sampled."
	fi
done


# STAGE 2: MEASURE REAL SPDS IN COMPLETE GENOMES ###################
# 2 GB RAM per CPU (core)
for s in ${strains[@]}; do
	checkpoint="${out_dir}/${s}/stage_2.success"
	if [ ! -f "$checkpoint" ]; then
		script="${out_dir}/${s}/run_Bandage_Ref.slurm"
		if [ ! -f "$script" ]; then
			echo "Creating a SLURM script for measuring SPDs in reference genomes of the strain ${s}."
			printf "#!/bin/bash\n\n" > $script  # The echo command does not interpret the newline character.
			printf "#SBATCH --partition=${par_short}\n" >> $script
			printf "#SBATCH --job-name='Bandage:${s}'\n" >> $script
            if [ "$add_account" = true ]; then
                printf "#SBATCH --account='${account}'\n" >> $script
            fi
			printf "#SBATCH --nodes=1\n#SBATCH --ntasks=1\n" >> $script
			printf "#SBATCH --cpus-per-task=${cores_short}\n" >> $script
			printf "#SBATCH --mem=${MEM_SHORT}\n" >> $script
			printf "#SBATCH --time=${WALLTIME_SHORT}\n\n" >> $script
            if [ ! -z "$blast_module" ]; then
                printf "module load ${blast_module}\n\n" >> $script
            fi
			for i in ${iterations[@]}; do
				printf "${bandage_dir}/Bandage distance ${ref_dir}/${s}.gfa ${out_dir}/${s}/${i}/${s}.${query_suffix} --evfilter 1e-5 --ifilter 100 --minpatcov 1 > ${out_dir}/${s}/${i}/${s}__dr.tsv &\n" >> $script
			done
			printf "\nwait\n" >> $script
		fi
		echo "Running Bandage to measure SPDs in complete genomes of strain ${s}."
		sbatch --wait $script  # until this job finishes, then move to the next strain
		touch "$checkpoint"
	else
		echo "Reference SPDs have been measured in complete genomes for strain ${s}."
	fi
done


# STAGE 3: MEASURE SPDS IN DE NOVO ASSEMBLIES ##################
# 2 GB RAM per CPU (core)
for s in ${strains[@]}; do
	checkpoint="${out_dir}/${s}/stage_3_${out_suffix}.success"
	if [ ! -f "$checkpoint" ]; then
		script="${out_dir}/${s}/run_Bandage_${out_suffix}.slurm"
		if [ ! -f "$script" ]; then
			echo "Creating a SLURM script for measuring SPDs in de novo assemblies of the strain ${s}."
			printf "#!/bin/bash\n\n" > $script  # The echo command does not interpret the newline character.
			printf "#SBATCH --partition=${par_short}\n" >> $script
			printf "#SBATCH --job-name='Bandage:${s}'\n" >> $script
            if [ "$add_account" = true ]; then
                printf "#SBATCH --account='${account}'\n" >> $script
            fi
			printf "#SBATCH --nodes=1\n#SBATCH --ntasks=1\n" >> $script
			printf "#SBATCH --cpus-per-task=${cores_short}\n" >> $script
			printf "#SBATCH --mem=${MEM_SHORT}\n" >> $script
			printf "#SBATCH --time=${WALLTIME_SHORT}\n\n" >> $script
			if [ ! -z "$blast_module" ]; then
                printf "module load ${blast_module}\n\n" >> $script
            fi
			for i in ${iterations[@]}; do
				printf "${bandage_dir}/Bandage distance ${graph_dir}/${s}.${graph_suffix} ${out_dir}/${s}/${i}/${s}.${query_suffix} ${bandage_args} > ${out_dir}/${s}/${i}/${s}__${out_suffix}.tsv &\n" >> $script
			done
			printf "\nwait\n" >> $script
		fi
		echo "Running Bandage to measure SPDs in de novo assemblies of strain ${s}."
		sbatch --wait $script
		touch "$checkpoint"
	else
		echo "SPDs have been measured for the strain ${s}."
	fi
done


# STAGE 4: PARSE BANDAGE OUTPUTS FOR REFERENCE SPDS ###############
for s in ${strains[@]}; do
	checkpoint="${out_dir}/${s}/stage_4.success"
	if [ ! -f "$checkpoint" ]; then
		script="${out_dir}/${s}/parse_dr.slurm"
		if [ ! -f "$script" ]; then
			echo "Creating a SLURM script for parsing reference SPDs of the strain ${s}."
			printf "#!/bin/bash\n\n" > $script  # The echo command does not interpret the newline character.
			printf "#SBATCH --partition=${par_short}\n" >> $script
			printf "#SBATCH --job-name='ParseDr:${s}'\n" >> $script
            if [ "$add_account" = true ]; then
                printf "#SBATCH --account='${account}'\n" >> $script
            fi
			printf "#SBATCH --nodes=1\n#SBATCH --ntasks=1\n" >> $script
			printf "#SBATCH --cpus-per-task=${cores_short}\n" >> $script
			printf "#SBATCH --mem=${MEM_SHORT}\n" >> $script
			printf "#SBATCH --time=${WALLTIME_SHORT}\n\n" >> $script
            if [ "$load_python" = true ]; then
                printf "module load ${python_module}\n\n" >> $script
            fi
			for i in ${iterations[@]}; do
				printf "python ${script_dir}/measurement/compile_dists.py -m ${out_dir}/${s}/${i}/${s}__dr.tsv -sf '.tsv' -d '__' -od ${out_dir}/${s}/${i}/Parsed -os dr.tsv -of dr_failures.tsv -oe dr_errors.tsv > ${out_dir}/${s}/${i}/Parsed/parse_dr.log &\n" >> $script
			done
			printf "\nwait\n" >> $script
		fi
		echo "Parsing reference SPDs of the strain ${s}."
		sbatch --wait $script
		touch "$checkpoint"
	else
		echo "Reference SPDs have been parsed for the strain ${s}."
	fi
done


# STAGE 5: PARSE BANDAGE OUTPUTS FOR QUERY SPDS ###############
for s in ${strains[@]}; do
	checkpoint="${out_dir}/${s}/stage_5_${out_suffix}.success"
	if [ ! -f "$checkpoint" ]; then
		script="${out_dir}/${s}/parse_SPDs_${out_suffix}.slurm"
		if [ ! -f "$script" ]; then
			echo "Creating a SLURM script for parsing SPDs of the strain ${s}."
			printf "#!/bin/bash\n\n" > $script  # The echo command does not interpret the newline character.
			printf "#SBATCH --partition=${par_short}\n" >> $script
			printf "#SBATCH --job-name='ParseSPD:${s}'\n" >> $script
            if [ "$add_account" = true ]; then
                printf "#SBATCH --account='${account}'\n" >> $script
            fi
			printf "#SBATCH --nodes=1\n#SBATCH --ntasks=1\n" >> $script
			printf "#SBATCH --cpus-per-task=${cores_short}\n" >> $script
			printf "#SBATCH --mem=${MEM_SHORT}\n" >> $script
			printf "#SBATCH --time=${WALLTIME_SHORT}\n\n" >> $script
			if [ "$load_python" = true ]; then
                printf "module load ${python_module}\n\n" >> $script
            fi
			for i in ${iterations[@]}; do
				printf "python ${script_dir}/measurement/compile_dists.py -m ${out_dir}/${s}/${i}/${s}__${out_suffix}.tsv -sf '.tsv' -d '__' -od ${out_dir}/${s}/${i}/Parsed -os d__${out_suffix}.tsv -of d_failures__${out_suffix}.tsv -oe d_errors__${out_suffix}.tsv > ${out_dir}/${s}/${i}/Parsed/parse_d__${out_suffix}.log &\n" >> $script
			done
			printf "\nwait\n" >> $script
		fi
		echo "Parsing SPDs of the strain ${s}."
		sbatch --wait $script  # Prepare and launch the next script when all jobs are finished.
		touch "$checkpoint"
	else
		echo "SPDs have been parsed for the strain ${s}."
	fi
done


# STAGE 6: MERGE TABLES OF REFERENCE SPDS AND THOSE FROM DE NOVO ASSEMBLIES PER STRAIN ###################
checkpoint="${out_dir}/pipeline.success"
script="${out_dir}/merge_SPDs_${out_suffix}.slurm"
if [ ! -f "$checkpoint" ]; then
	if [ ! -f "$script" ]; then
		echo "Creating a SLURM script for merging SPDs of the strain ${s}."
		printf "#!/bin/bash\n\n" > $script  # The echo command does not interpret the newline character.
		printf "#SBATCH --partition=${par_long}\n" >> $script
		printf "#SBATCH --job-name='MergeSPD'\n" >> $script
        if [ "$add_account" = true ]; then
            printf "#SBATCH --account='${account}'\n" >> $script
        fi
		printf "#SBATCH --nodes=1\n#SBATCH --ntasks=1\n" >> $script
		printf "#SBATCH --cpus-per-task=${cores_long}\n" >> $script
		printf "#SBATCH --mem=${MEM_LONG}\n" >> $script
		printf "#SBATCH --time=${WALLTIME_LONG}\n\n" >> $script
		if [ "$load_R" = true ]; then
            printf "module load ${R_module}\n\n" >> $script
        fi
        for s in ${strains[@]}; do
            dr_tables="${out_dir}/${s}/1/Parsed/dr.tsv"
            d_tables="${out_dir}/${s}/1/Parsed/d__${out_suffix}.tsv"
            if [ "$iter_num" -gt "1" ]; then
                for i in `seq 2 $iter_num`; do
                    dr_tables="${dr_tables},${out_dir}/${s}/${i}/Parsed/dr.tsv"
                    d_tables="${d_tables},${out_dir}/${s}/${i}/Parsed/d__${out_suffix}.tsv"
                done
            fi
            printf "Rscript --vanilla ${script_dir}/analysis/merge_dr_ds.R --dr ${dr_tables} --d ${d_tables} --d_max ${d_max} --n_max ${n_max} --out ${out_dir}/Merged_${out_suffix}/${s}__Dm.tsv &\n" >> $script
        done
		printf "\nwait\n" >> $script
	fi
	echo "Merging reference and query distances for all strains."
	sbatch --wait $script
	touch "$checkpoint"
fi


echo "Congratulations! The pipeline has been run through successfully!"
