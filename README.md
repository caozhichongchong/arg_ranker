# arg_ranker
arg_ranker evaluates the risk of ARGs in genomes and metagenomes

## Install
`pip install arg_ranker`

## Requirement
* python 3
* kraken2: `conda install -c bioconda kraken2`\
download kraken2 database: `kraken2-build --standard --db $KRAKENDB` \
where $krakenDB is your preferred database name/location\
* diamond: `conda install -c bioconda diamond`\
* blast+: `conda install -c bioconda blast`

## How to use it
* put all your genomes (.fa or .fasta) and metagenomes (.fq or .fastq) into one folder ($INPUT)
* run `arg_ranker -i $INPUT --kkdb $KRAKENDB`
* run `sh arg_ranking/script_output/arg_ranker.sh`

## Output
* Sample_ranking_results.txt (Table 1)

    |Sample|Rank_I_per|Rank_II_per|Rank_III_per|Rank_IV_per|Unassessed_per|Total_abu|Rank_code|Rank_I_risk|Rank_II_risk|Rank_III_risk|Rank_IV_risk|ARGs_unassessed_risk|note1|
    | :--------: | :--------: | :--------: | :--------: | :--------: | :--------: | :--------: | :--------: | :--------: | :--------: | :--------: | :--------: | :--------: | :--------: |
    |WEE300_all-trimmed-decont_1.fastq|4.6E-02|0.0E+00|6.8E-02|7.5E-01|1.3E-01|5.4E-04|1.0-0.0-0.5-1.7-0.3|1.0|0.0|0.5|1.7|0.3|hospital_metagenome|
    |EsCo_genome.fasta|0.0E+00|0.0E+00|0.0E+00|1.0E+00|0.0E+00|2.0E+00|0.0-0.0-0.0-2.2-0.0|0.0|0.0|0.0|2.2|0.0|E.coli_genome|

1. We compute the abundance of ARGs as the copy number of ARGs divided by the 16S copy number in a sample\
For metagenomes, the copy number of ARGs and 16S were computed as the number of reads mapped to them divided by their gene length.\
Rank_I_per - Unassessed_per: percentage of ARGs of a risk rank\
Total_abu: total abundance of all ARGs
2. We compute the risk of ARGs as the average abundance of ARGs of a risk rank divided the average abundance of all ARGs\
Rank_I_risk - Unassessed_risk: the risk of ARGs of a risk rank\
Rank_code: a code of ARG risk from Rank I to Unassessed

* Sample_ARGpresence.txt:\
The abundance, the gene family, and the antibiotic of resistance of ARGs detected in the input samples

## Test
run `arg_ranker -i example --kkdb $KRAKENDB`\
run `sh arg_ranking/script_output/arg_ranker.sh`\
The arg_ranking/Sample_ranking_results.txt should look like Table 1

## Metadata for your samples (optional)
arg_ranker can merge your sample metadata into the results of ARG ranking (i.e. note1 in Table 1).\
Simply put all information you would like to include into a tab-delimited table\
Make sure that your sample names are listed as the first column (check example/metadata.txt).

## Copyright
Dr. An-Ni Zhang (MIT), Prof. Eric Alm (MIT), Prof. Tong Zhang* (University of Hong Kong)

## Citation
1. Zhang AN, ..., Alm EJ, Zhang T: Choosing Your Battles: Which Resistance Genes Warrant Global Action? (bioRxiv coming soon)
2. Yang Y, ..., Tiedje JM, Zhang T: ARGs-OAP: online analysis pipeline for antibiotic resistance genes detection from metagenomic data using an integrated structured ARG-database. Bioinformatics 2016.

## Contact
anniz44@mit.edu or caozhichongchong@gmail.com
