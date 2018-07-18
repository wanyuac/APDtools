#!/bin/bash
# This is a handy little bash script that submits jobs through SLURM so that you do not need to rerun dist_from_graphs.py to get slurm files
# submitted. This is particularly useful when you finished running dist_from_graphs.py under a debugging mode (flagged by the --debug) option.
# Usage example: bash submit_jobs.sh
# Licence: GNU GPL 2.1
# Author: Yu Wan (wanyuac@gmail.com)
# Development history: 21/11/2016

for script in $(ls -1 *.slurm); do
	sbatch $script
done
