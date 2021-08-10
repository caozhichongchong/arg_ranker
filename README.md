# arg_ranker
arg_ranker evaluates the risk of ARGs in genomes and metagenomes

## Install
`pip install arg_ranker`

## Requirement
* python 3
* diamond: `conda install -c bioconda diamond` (https://github.com/bbuchfink/diamond)
* blast+: `conda install -c bioconda blast` (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* For metagenomes:
    * kraken2: `conda install -c bioconda kraken2`(https://github.com/DerrickWood/kraken2/wiki)
        * to compute the abundance of ARGs as copy number of ARGs per bacterial cell (recommended)
            * download the kraken2 standard database (50 GB of disk space): `kraken2-build --standard --db $KRAKENDB` \
            where $KRAKENDB is your preferred database name/location
            * MicrobeCensu: `git clone https://github.com/snayfach/MicrobeCensus && cd MicrobeCensus && python setup.py install` to estimate the average genome size for metagenomes.
            (https://github.com/snayfach/MicrobeCensus)
        * to compute the abundance of ARGs as copy number of ARGs per 16S
            * download the kraken2 16S database (73.2 MB of disk space): `kraken2-build --db $DBNAME --special greengenes`

## How to use it
* put all your genomes (.fa or .fasta) and metagenomes (.fq or .fastq) into one folder ($INPUT)
* run `arg_ranker -i $INPUT` (genomes only)
* run `arg_ranker -i $INPUT -kkdb $KRAKENDB` (genomes/metagenomes + kraken2 standard database)
    * or run `arg_ranker -i $INPUT -kkdb $KRAKENDB -kkdbtype 16S` (kraken2 16S database)
* run `sh arg_ranking/script_output/arg_ranker.sh`

## Output
* Sample_ranking_results.txt (Table 1)

    |Sample|Rank_I_per|Rank_II_per|Rank_III_per|Rank_IV_per|Unassessed_per|Total_abu|Rank_code|Rank_I_risk|Rank_II_risk|Rank_III_risk|Rank_IV_risk|ARGs_unassessed_risk|note1|
    | :--------: | :--------: | :--------: | :--------: | :--------: | :--------: | :--------: | :--------: | :--------: | :--------: | :--------: | :--------: | :--------: | :--------: |
    |WEE300_all-trimmed-decont_1.fastq|4.6E-02|0.0E+00|6.8E-02|7.5E-01|1.3E-01|5.4E-04|1.0-0.0-0.5-1.7-0.3|1.0|0.0|0.5|1.7|0.3|hospital_metagenome|
    |EsCo_genome.fasta|0.0E+00|0.0E+00|0.0E+00|1.0E+00|0.0E+00|2.0E+00|0.0-0.0-0.0-2.2-0.0|0.0|0.0|0.0|2.2|0.0|E.coli_genome|

1. Rank_I_per - Unassessed_per: percentage of ARGs of a risk Rank\
Total_abu: total abundance of all ARGs
2. For genomes, we output the copy number of ARGs detected in each genome.
3. For metagenomes, we compute the abundance of ARGs as the copy number of ARGs divided by the bacterial cell number or 16S copy number in the same metagenome.\
If you downloaded the kraken2 standard database, we compute the copy number of ARGs divided by the bacterial cell number.\
If you downloaded the kraken2 16S database, we compute the copy number of ARGs divided by the 16S copy number.\
The copy number of ARGs, 16S, and bacterial cells were computed as the number of reads mapped to them divided by their gene/genome length.
4. We compute the contribution of each ARG risk Rank as the average abundance of ARGs of a risk Rank divided by the average abundance of all ARGs\
Rank_I_risk - Unassessed_risk: the contribution of ARGs of a risk Rank\
Rank_code: a code of contribution from Rank I to Unassessed

* Sample_ARGpresence.txt:\
The abundance, the gene family, and the antibiotic of resistance of ARGs detected in the input samples

## Test
run `arg_ranker -i example -kkdb $KRAKENDB`\
run `sh arg_ranking/script_output/arg_ranker.sh`\
The arg_ranking/Sample_ranking_results.txt should look like Table 1 (using kraken2 standard database)

## Metadata for your samples (optional)
arg_ranker can merge your sample metadata into the results of ARG ranking (i.e. note1 in Table 1).\
Simply put all information you would like to include into a tab-delimited table\
Make sure that your sample names are listed as the first column (check example/metadata.txt).

## Copyright
Dr. An-Ni Zhang (MIT), Prof. Eric Alm (MIT), Prof. Tong Zhang* (University of Hong Kong)

## Citation
Zhang, AN., Gaston, J.M., Dai, C.L. et al. An omics-based framework for assessing the health risk of antimicrobial resistance genes. Nat Commun 12, 4765 (2021). https://doi.org/10.1038/s41467-021-25096-3

## Contact
anniz44@mit.edu or caozhichongchong@gmail.com
