# ARG_Ranker
Install
pip install arg_ranker

Sample_ranking.py evaluates and assigns the risk and priority levels to environmental samples
based on their profile of antibiotic resistant genes (ARGs).
Requirement: python packages (os, csv, argparse)
Requirement: a mothertable of the ARG abundance in all your samples
annotated by ARGs-OAP v1.0 (see example/All_sample_cellnumber.txt).
Optimal: a table of the metadata of your samples (see example/All_sample_metadata.txt).
Recommend unit of ARG abundance is copy per cell.

The mother table can be obtained by ARGs-OAP v1.0.
'perl dir_to_ARGs_OAP/database/DB/gene_subtype_type_ppm_16s_cellnumber_differentversion.pl dir_to_your_result/extracted.fa.blast6out.txt dir_to_your_result/meta_data_online.txt 25 1e-7 80 output_mothertable_name dir_to_ARGs_OAP/database/DB/other/structure_20160422.list dir_to_ARGs_OAP/database/DB/other/SARG-20160422.mod.fasta'
The mother table to input is output_mothertable_name.normalize_cellnumber.gene.tab

Copyright:An-Ni Zhang, Prof. Tong Zhang, University of Hong Kong
Citation:
1. This study
2. Yang Y, Jiang X, Chai B, Ma L, Li B, Zhang A, Cole JR, Tiedje JM, Zhang T: ARGs-OAP: online analysispipeline for antibiotic resistance genes detection from metagenomic data using an integrated structured ARG-database. Bioinformatics 2016. (optional: antibiotic resistance database)
Contact caozhichongchong@gmail.com
