# Determining Thresholds for Obtaining Accurate Measurement of the Shortest-path Distances via a Simulation-Validation Approach

This directory offers scripts helping users to determine appropriate thresholds for filtering out shortest-path distances (SPDs) that are believed to have a high inaccuracy rate. Users may refer to our GeneMates paper for the method, which consists of read simulation using [readSimulator](https://github.com/wanyuac/readSimulator), *de novo* genome assembly by [SPAdes](http://cab.spbu.ru/software/spades/), SPD measurement with Bandage, and accuracy evaluation using scripts from this directory. At present, the evaluation concerns two aspects:

- Accuracy as a function of node numbers in paths in which SPDs are measured.
- Accuracy as a function of the upper bound of SPDs.

This page explains functionality of scripts used for this evaluation. Here, an accurate SPD is referred to as an SPD that deviates from its corresponding true distance by no more than _e_ bp (by the absolute value of the error), where _e_ is a predefined tolerance of errors. (See the GeneMates paper for details)

Note that users do not have to use the method demonstrated here for their studies, and we welcome suggestions about other methods. Since some scripts (such as ``SPDcomp.sh`) is platform dependent, users are encouraged to modify the code for their computer systems.

<br/>

## 1. Evaluating accuracy as a function of node numbers of SPDs

Empirically, as we have demonstrated in the GeneMates paper, we found that SPDs measured in two connected nodes generally reach accuracy of more than 90%. Script `accuracy_vs_nodes.R` is designed to perform the calculation of accuracy given a series of node numbers and plot the accuracy as a line graph (accuracy ~ SPD, given a node number). This script takes as input a distance table generated by `merge_dist_tab.R`. This table matches SPDs (_d_) with their true distances (_dr_) calculated in complete reference genomes.

<br/>

## 2. Evaluating accuracy as a function of the upper bound of SPDs

Script `accuracy_vs_dist.R` is designed for this evaluation and it also takes as input the distance table generated by `merge_dist_tab.R`. Outputs include a list of accuracy given a maximum SPD and graphs visualising the accuracy.

<br/>

## 3. Other scripts

- `merge_dist_tab_batch.sh`: Generates and submits SLURM jobs of `merge_dist_tab.R` for multiple reference genomes.
- `SPDcomp.sh`: This is a pipeline measuring reference physical distances (the true distances) and SPDs using Bandage and compiling the measurements using script `merge_dr_ds.R`.
- `d_vs_dr.R`: This is another but simpler script comparing SPDs (measured from assembly graphs) to their true distances given a maximum SPD and a maximum node number. This script does not produce any figure but two tab-delimited text files.