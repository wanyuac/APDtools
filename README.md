# Measurement and Analysis of Allelic Physical Distances

APDtools is a stand-alone, optional component of the [GeneMates](https://github.com/wanyuac/GeneMates) package. It offers helper scripts for the measurement and analysis of allelic physical distances (APDs). APDtools implements the following functionality:

- Measurement of APDs in genome assemblies; (`measurement`)
  - Measurement of APDs from complete genome assemblies; (`measurement/ref`)
- Accuracy analysis of shortest-path distances (SPDs) by comparing SPDs to true APDs obtained from complete genome assemblies; (`accuracy`)
- Establishment of relationships between APDs and other biological data for interpretation. (`biology`)

Please see `README.md` of each subdirectory for details.

<br/>

**Dependencies**

- R
- Linux bash
- Python (version 3 is recommended)

<br/>

**Terminology**

- APD: The physical distance (in bp) between two alleles of one or two genes in a genome assembly.
- SPD: A particular kind of APDs. An SPD equals the number of base pairs between two target alleles following the shortest path in an assembly graph.
  - SPD = 0 when two alleles overlap.
  - The SPD is the true APD in complete (finished-grade) genome assemblies.
  - In draft assembly graphs, we consider SPDs as approximations of true APDs.

<br/>

**Citation**

Wan, Y., Wick, R. R., Zobel, J., Ingle, D. J., Inouye, M., & Holt, K. E. (2020). GeneMates: an R package for Detecting Horizontal Gene Co-transfer between Bacteria Using Gene-gene Associations Controlled for Population Structure. *BioRxiv*, 2020.02.29.970970. https://doi.org/10.1101/2020.02.29.970970.


