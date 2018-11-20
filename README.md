# arg_ranker
##Install
pip install arg_ranker
conda install -c caozhichongchong arg_ranker 

##Test (any of these two commands)
arg_ranker -i example/ARGprofile_example_1.txt -m example/metadata.txt
arg_ranker -i example/ARGprofile_example_2.txt -m example/metadata.txt

##Availability
https://anaconda.org/caozhichongchong/arg_ranker
https://pypi.org/project/arg-ranker/

##How to use it
1. Prepare your ARG profile
arg_ranker is suitable for the units of ppm, gene copy per 16S or gene copy per cell

Option 1: Run your own pipeline against our database ()
Format your results into example/ARGprofile_example_1.txt or example/ARGprofile_example_2.txt

Option 2: Run ARGs-OAP v1.0 () and one more step:
'perl dir_to_ARGs_OAP/database/DB/gene_subtype_type_ppm_16s_cellnumber_differentversion.pl dir_to_your_result/extracted.fa.blast6out.txt dir_to_your_result/meta_data_online.txt 25 1e-7 80 your_output_name dir_to_ARGs_OAP/database/DB/other/structure_20160422.list dir_to_ARGs_OAP/database/DB/other/SARG-20160422.mod.fasta'
The ARG profile to input is your_output_name.normalize_cellnumber.gene.tab

2. Prepare your metadata (optional)
Format your metadata of metagenomic samples into example/metadata.txt (not necessarily the same)
First column matches the sample ID in your ARG profile;
Other columns contain the metadata of your samples (such as habitat/eco-type, accession number, group...)

##Introduction
Sample_ranking.py evaluates and assigns the risk and priority levels to environmental samples
based on their profile of antibiotic resistant genes (ARGs).
Requirement: python packages (os, csv, argparse)
Requirement: a mothertable of the ARG abundance in all your samples
annotated by ARGs-OAP v1.0 (see example/All_sample_cellnumber.txt).
Optimal: a table of the metadata of your samples (see example/All_sample_metadata.txt).
Recommend unit of ARG abundance is copy per cell.

##Copyright
Copyright:An-Ni Zhang, Prof. Tong Zhang, University of Hong Kong
Citation:
1. This study
2. Yang Y, Jiang X, Chai B, Ma L, Li B, Zhang A, Cole JR, Tiedje JM, Zhang T: ARGs-OAP: online analysis pipeline for antibiotic resistance genes detection from metagenomic data using an integrated structured ARG-database. Bioinformatics 2016. (optional: antibiotic resistance database)
Contact caozhichongchong@gmail.com
