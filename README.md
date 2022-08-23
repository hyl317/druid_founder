# DRUID with extension for founder populations
## Introduction
This is an extension to [DRUID](https://www.sciencedirect.com/science/article/pii/S0002929718301691) (deep relatedness utilizing identity by descent) that additionally takes as input a file describing the population size history in the recent past so that excessive IBD sharing due to small population size can be corrected when inferring pedigree relationships.

## Prerequisite 
Numpy, networkx

## Usage
To view a list of all available command options, run

    python DRUID.py -h

In the most basic use case, three files -i, -s, -m must be provided. 

The -i file gives the pairwise IBD1 and IBD2 proportion . It should have 4 columns, the first two columns being the sample ID, the third column being the IBD1 proportion (IBD1 length/total genome length), and the fourth column being the IBD2 proportion. For example,

    sample1 sample2 0.5 0.25

The -s file gives the pairwise segment file. It should have 6 column, the first two columns being the sample ID, the third column chromosome index, the fourth column being a string describing the type of IBD segment (one of "IBD1" or "IBD2"), and the 5th and 6th columns being the starting and end point of this IBD segment (in cM). For example,

    sample1 sample2 1   IBD2    15.758  25.896

The -m file gives a genetic map. It has 4 columns. The first column is the chromosome index (must match the chromosome index in -s file). The second column is dbSNP ID (can be simply replaced by "."). The third and fourth columns are the genetic (in cM) and physical positions, respectively. For example,

    1       .       2.014081305147059       753498
    1       .       2.0195751834862388      766698
    1       .       2.0202208955007257      767969
    1       .       2.0219964441805227      771799
    1       .       2.0278411895573214      784624

The default for DRUID is to use the simple kinship coefficient to determine close relatives and isolated distant relatives. However, with --useERSA option, it can use the more complicated model described in [ERSA](https://genome.cshlp.org/content/21/5/768.full.html), which has higher power than the simple approach (at the expense of runtime). With --useERSA option set to true, it is recommend to set --FDR as well. It will allow DRUID to use the Benjamini-hochberg procedure to control false discovery rate (by default set to 5%, adjustable by the option --alpha). Otherwise it will use the Bonferroni correction, which is far less powerful. In addition, with --useERSA option, a file in hapIBD output format must be provided (--hapibd), even though this is somewhat redundaunt given that a .seg file is already provided.

To use the correction for foudner population, one must additionally provide a file describing the population size history within the last 100-200 generations. Each line should have two columns, the first column being the generation number and the second column being the diploid population size. Lines that start with # will be ignored. Examples of these Ne files can be found in the directory ./sampleNe. In addition, one needs to specify the IBD length threshold in .seg file with the command line option --minIBD. This parameters gives the minimum length of IBD segments you use here.

