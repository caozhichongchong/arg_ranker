# arg_ranker

## Install
`pip install arg_ranker`

`conda install -c caozhichongchong arg_ranker`

## Test (download examples and use any of these commands)
`arg_ranker -i example/ARGprofile_example_1.txt -m example/metadata.txt`\
`arg_ranker -i example/ARGprofile_example_2.txt -m example/metadata.txt`\
`arg_ranker -i test`

## How to use it
### Prepare your ARG profile

arg_ranker is suitable for the units of ppm, gene copy per 16S or gene copy per cell

#### Option 1: Use our pipeline

1. Use my traits_finder to search ARGs in genomes and metagenomes (in preparation)\
Now we have both nucleotides and amino acids databases!\ 
https://github.com/caozhichongchong/traits_finder

2. Run\
`arg_ranker -i ARG.profile.txt -m metadata.txt`\
`arg_ranker -i ARG.profile.txt`

#### Option 2: Run your own pipeline using our database

1. Search ARGs-OAP v1.0 database (amino acids) in your data using diamond or blast\
https://github.com/caozhichongchong/arg_ranker/tree/master/arg_ranker/data/SARG.db.fasta*

2. Format your results into example/ARGprofile_example_1.txt or example/ARGprofile_example_2.txt

3. Run\
`arg_ranker -i ARG.profile.txt -m metadata.txt`\
`arg_ranker -i ARG.profile.txt`\
If you see a lot of errors saying: "ARGs in mothertable do not match with the ARGs in ARG_rank.txt.\
Please check something something in ARG.summary.cell.txt!"\
It means that the samples are placed as row names instead of colomn names (which arg_ranker expects).\
Don't worry, please try: `arg_ranker -i ARG.profile.txt.t`\
As we automatically transpose your table to make it work.

#### Option 3: Use results from ARGs-OAP v1.0 (not recommended)

1. If you have already run the ARGs-OAP v1.0 pipeline\
    https://github.com/biofuture/Ublastx_stageone/archive/Ublastx_stageone.tar.gz\
    https://github.com/biofuture/Ublastx_stageone/archive/Ublastx_stageone.zip

2. Check the "extracted.fa.blast6out.txt" and "meta_data_online.txt" in the output_dir

3. Run\
`arg_ranker -f True -fo output_dir`\
`arg_ranker -i formated_table.normalize_cellnumber.gene.tab -m metadata.txt`

### Prepare your metadata for your samples (optional)

Format your metadata of metagenomic samples into example/metadata.txt (not necessarily the same)\
First column matches the sample ID in your ARG profile;\
Other columns contain the metadata of your samples (such as habitat/eco-type, accession number, group...)

## Introduction
ARG_ranker evaluates the risk of antibiotic resistance in metagenomes.\
We designed a framework to rank the risk of ARGs based on three factors: “anthropogenic enrichment”, “mobility”, and “host pathogenicity”, informed by all available bacterial genomes, plasmids, integrons, and 850 metagenomes covering diverse global eco-habitats. The framework prioritizes 3% of ARGs in Rank I (the most at risk of dissemination among pathogens) and 0.3% of ARGs in Rank II (high potential emergence of new resistance in pathogens). 

Requirement: python packages (pandas, argparse)

Requirement: a mothertable of the ARG abundance in all your samples
annotated by ARGs-OAP v1.0 \
(see example/All_sample_cellnumber.txt).

Optimal: a table of the metadata of your samples \
(see example/All_sample_metadata.txt).

## Copyright
Dr. An-Ni Zhang (MIT), Prof. Eric Alm (MIT), Prof. Tong Zhang* (University of Hong Kong)

## Citation
1. Zhang AN, ..., Alm EJ, Zhang T: Choosing Your Battles: Which Resistance Genes Warrant Global Action? (bioRxiv coming soon)
2. Yang Y, ..., Tiedje JM, Zhang T: ARGs-OAP: online analysis pipeline for antibiotic resistance genes detection from metagenomic data using an integrated structured ARG-database. Bioinformatics 2016.

## Contact
anniz44@mit.edu or caozhichongchong@gmail.com
