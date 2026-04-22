The repository mainly contains a Snakemake pipeline to process HiC data from .bam format after processing with HiCUP.
The scripts of this pipeline are present in workflow/scripts.

The data is processed by the pipeline to a .cool file first. This file is separated into chromosome-wise contact matrices of the desired resolution and ICE normalized in python. The ICE normalized matrix is then fed to an implementation of the Kuramoto model in Julia present in the ICE_processor.jl script. Finally the order parameter R is seen while the parameter K of Kuramoto model is varied, and the K value where R = 0.5 is observed and saved. The final step is done using R scripts.

The exact scripts used can be identified by looking at the Snakemake rules in workflow/rules

Other scripts also used in processing HiC data and related biological data are present in Scripts.
