# DRUID with extension for founder populations
## Introduction
This is an extension to [DRUID](https://www.sciencedirect.com/science/article/pii/S0002929718301691) (deep relatedness utilizing identity by descent) that additionally takes as input a file describing the population size history in the recent past so that excessive IBD sharing due to small population size can be corrected when inferring pedigree relationships.

## Prerequisite 
Numpy, networkx

## Usage
To view a list of all available command options, run

    python DRUID.py -h

In the most basic use case, three files -i, -s, -m must be provided. Descriptions of these three files can be found [here](https://github.com/williamslab/druid).

The default for DRUID is to use the simple kinship coefficient to determine close relatives and isolated distant relatives. However, with --useERSA option, it can use the more complicated model described in [ERSA](https://genome.cshlp.org/content/21/5/768.full.html), which has higher power than the simple approach (at the expense of runtime). With --useERSA option set to true, it is recommend to set --FDR as well. It will allow DRUID to use the Benjamini-hochberg procedure to control false discovery rate (by default set to 5%, adjustable by the option --alpha). Otherwise it will use the Bonferroni correction, which is far less powerful. In addition, with --useERSA option, a file in hapIBD output format must be provided (--hapibd), even though this is somewhat redundaunt given that a .seg file is already provided.

To use the correction for foudner population, one must additionally provide a file describing the population size history within the last 100-200 generations. Each line should have two columns, the first column being the generation number and the second column being the diploid population size. Examples of these Ne files can be found in the directory ./sampleNe. In addition, one needs to specify the IBD length threshold in .seg file with the command line option --minIBD. This parameters gives the minimum length of IBD segments you use here.

