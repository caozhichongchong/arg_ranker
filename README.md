# arg_ranker

## Install
pip install arg_ranker

conda install -c caozhichongchong arg_ranker

## Test (any of these two commands)
`arg_ranker -i example/ARGprofile_example_1.txt -m example/metadata.txt`\
`arg_ranker -i example/ARGprofile_example_2.txt -m example/metadata.txt`

## How to use it
### Prepare your ARG profile

arg_ranker is suitable for the units of ppm, gene copy per 16S or gene copy per cell

#### Option 1: Run your own pipeline against our database

1. Download the ARGs-OAP v1.0 database\
https://github.com/caozhichongchong/arg_ranker/tree/master/arg_ranker/data/SARG.db.fasta*

2. Format your results into example/ARGprofile_example_1.txt or example/ARGprofile_example_2.txt

#### Option 2: Run ARGs-OAP v1.0 and format the results by ARG_Ranker

1. Download ARGs-OAP v1.0 pipeline and run the pipeline\
    https://github.com/biofuture/Ublastx_stageone/archive/Ublastx_stageone.tar.gz
    https://github.com/biofuture/Ublastx_stageone/archive/Ublastx_stageone.zip

    A brief introduction on how to use ARGs-OAP v1.0\
    Please refer to the README.md of ARGs-OAP v1.0 for more details

    Prepare your metadata for your samples into example/metadata.txt (separated by tab)\
    SampleID (a number for the sample) | Name (metagenomic samples name) | Category (metadata of habitat, or group)\
    `./ublastx_stage_one  -i inputfqs -o testoutdir -m meta-data.txt -n 2`

        Usage: ./ublastx_stage_one -i <Fq input dir> -m <Metadata_map.txt> -o <output dir>
        -n [number of threads] -f [fa|fq] -z -h  -c   
            -i Input files directory, required\
            -m meta data file, required
            -o Output files directory, default current directory
            -n number of threads used for usearch, default 1
            -f the format of processed files, default fq
            -z whether the fq files were .gz format, if -z, then firstly gzip -d, default(none)
            -c This option fulfill copy number correction by Copywriter database to transfrom 16S information into cell number [ direct searching hyper variable region database by usearch; default 1]
            -h print this help information

2. Check the "extracted.fa.blast6out.txt" and "meta_data_online.txt" in the output_dir

3. Run\
`arg_ranker -f True -fo output_dir`\
`arg_ranker -i formated_table.normalize_cellnumber.gene.tab -m metadata.txt`

### Prepare your metadata for your samples (optional)

Format your metadata of metagenomic samples into example/metadata.txt (not necessarily the same)\
First column matches the sample ID in your ARG profile;\
Other columns contain the metadata of your samples (such as habitat/eco-type, accession number, group...)

## Introduction
Sample_ranking.py evaluates and assigns the risk and priority levels to environmental samples
based on their profile of antibiotic resistant genes (ARGs).

Requirement: python packages (pandas, argparse)

Requirement: a mothertable of the ARG abundance in all your samples
annotated by ARGs-OAP v1.0 (see example/All_sample_cellnumber.txt).

Optimal: a table of the metadata of your samples (see example/All_sample_metadata.txt).

## Copyright
Dr. An-Ni Zhang (MIT), Prof. Tong Zhang (University of Hong Kong)

## Citation
1. Zhang AN, ..., Alm EJ, Zhang T: Whom to Fight: Top Risk Antibiotic Resistances for Global Action (Under Review)
2. (Optional: antibiotic resistance database)\
Yang Y, ..., Tiedje JM, Zhang T: ARGs-OAP: online analysis pipeline for antibiotic resistance genes detection from metagenomic data using an integrated structured ARG-database. Bioinformatics 2016.

## Contact
anniz44@mit.edu or caozhichongchong@gmail.com
