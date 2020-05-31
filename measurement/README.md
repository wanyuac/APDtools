# Measurement of Shortest-path Distances in Genome Assemblies

Sections of this page describes acquisition of shortest-path distances (SPDs) from different inputs. [Bandage (distance)](https://github.com/wanyuac/Bandage/tree/distance) is required for measuring the distances in draft genome assembly graphs.

<br/>

## 1. A generic approach for any kind of genome assemblies

Steps for measuring SPDs in an assembly graph are described as follows. In addition to GFA files, particularly, a FASTA file stores an assembly graph that has no connection between nodes (that is, contigs in the FASTA file).



### 1.1. Preparing a table of target alleles per sample genome

R function `mkFilterTSV` in GeneMates creates a two-column, tab-delimited table of target alleles from an allelic presence-absence matrix for subsequent distance measurement. The table follows the format:

```bash
[Genome name]\t[comma-delimited vector of allele names]
```

Allele names much be the same as those in an input multi-FASTA file of allele sequences. Users may use their own scripts to generate a target table in the same format.



### 1.2. Measuring SPDs in assembly graphs

Script `dist_from_graphs.py` runs `Bandage (distance)` to measure the SPDs. An example command is:

```bash
python ./APDtools/measurement/dist_from_graphs.py --graphs ./assembly/*.gfa --queries ./alleles/*_all_consensus.fna --filter target_alleles_perGenome.tsv --ext_id --suffix_graphs .gfa --suffix_queries _all_consensus.fna --suffix_out dataset1 --outdir output --bandage ./Bandage/Bandage --mem 4096 --walltime "0-0:30:0" --par group1 --other_args "--evfilter 1e-10 --ifilter 95 --minpatcov 0.95" > dataset1_SPDs.log
```

FASTA files are also legit arguments for the `graphs` parameter. This script only supports the SLURM job scheduler on a computer cluster. A simplified, single-job version of this script is to be created for a broader range of scenarios. (Please feel free to make your own versions)



### 1.3. Parsing and compiling Bandage output files

Raw Bandage outputs, whose format is not program friendly, needs to be parsed into pure tab-delimited files and compiled into a single table for subsequent analysis. Script `compile_dists.py` is developed for this purpose. Example command:

```bash
python ./APDtools/measurement/compile_dists.py -m output/*.tsv -sf .tsv -d "__" -od output/compiled -os SPDs.tsv -of failed_query_paths.tsv -oe errors.tsv > log/parse_dists.log
```



### 1.4. APD prioritisation

This step merges APDs from different sources in accordance with user's predefined weights of the sources (such as draft assembly graphs, contigs, and complete genome assemblies). The weight is proportional to user's confidence to the accuracy of SPDs from a specific source. Our empirical study has shown that this step improves the overall accuracy of SPD measurements. 

```bash
# Uses two cores for computation
Rscript --vanilla ./APDtools/measurement/prioritise_dists.R -d output/compiled/APDs.tsv -w output/compiled/source_weights.tsv -o output/compiled/SPDs_merg.tsv -c 2
```

The resulting `SPDs_merg.tsv` can be used as an input (specified by argument `phys.dists`) of the function `findPhysLink` in GeneMates, after filtering out unreliable SPDs based on user's predefined thresholds. (See our GeneMates paper and descriptions in the directory [`accuracy`](https://github.com/wanyuac/APDtools/tree/master/accuracy)) Since our empirical study also has shown that the reliability of SPDs usually follow the order: complete genomes > contigs > draft assembly graphs, the tab-delimited weight file `source_weights.tsv` has the content:

```bash
source	weight
graph	1
contig	2
complete	3
```

Users are advised to choose appropriate weights based on their data.

<br/>

## 2. For complete genome assemblies

It is straightforward to calculate SPDs from coordinates of genetic features in a linear genome so Bandage is not necessary. Nevertheless, since FASTA files store DNA and amino acid sequences in a linear manner, direct calculation of the SPD in a circular genome may result in a large error when two alleles are far apart. Specifically, every SPD in a circular genome must not exceed 1/2 of the genome length. To address this challenge, APDtools offers two approaches to the calculation of SPDs in complete assemblies of circular genomes.



### 2.1. Using script `calc_dr.R`

This script (in subdirectory `ref`) reads a tab-delimited coordinate table of gene features and calculates SPDs under the assumption of circular genome topology. Although it offers a fast and low-cost way for the acquisition of SPDs from complete genomes, it requires an argument of genome size, which is not convenient when a user has a number of genomes to analyse.

The table of gene features can be generated using scripts `gbk2table.sh` in `./ref` and [`gbk2tsv.py`](https://github.com/wanyuac/BINF_toolkit/blob/master/gbk2tsv.py) in BINF_toolkit. (Some manual reformatting may apply to their inputs) Script `calc_dr_slurm.sh` submits jobs of `calc_dr.R` to an SLURM job scheduler on a computer cluster for processing a number of reference genomes in parallel.



### 2.2. Using Bandage (distance)

An alternative approach to the calculation is using `Bandage (distance)`. The first step is to convert a FASTA file (may contain multiple contigs) to a GFA file (the format for assembly graphs) with the script `mfasta2gfa.py`. This step is only correct when every sequence in the FASTA file comes from a circular genome. Then the procedure of distance measurement is the same as measuring SPDs in any assembly graph. (See Section 1)

