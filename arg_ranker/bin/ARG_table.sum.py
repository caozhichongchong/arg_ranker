import argparse
import os
import glob

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                    help="input folder of ARG searching result",
                     type=str, default='arg_ranking/search_output',metavar='search_output')
parser.add_argument("-d",
                    help="database mapping file",
                     type=str, default='SARG.structure.txt ',metavar='SARG.structure.txt')
optional.add_argument('-kkdbtype',
                          help="Optional: type of kraken2 database (default = standard)",
                          choice = ['standard','16S'],
                          metavar="standard or 16S",
                          action='store', default='standard', type=str)

################################################## Definition ########################################################
args = parser.parse_args()
allsearchoutput = glob.glob(os.path.join(args.i,'*.blast.txt.filter'))
################################################### Function #######################################################
def load_ARG_mapping():
    ARG_mapping = dict()
    for lines in open(args.d,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        ARG,genotype,phenotype = lines_set[0:3]
        ARG_mapping.setdefault(ARG,'%s\t%s\t%s'%(ARG,genotype,phenotype))
    return ARG_mapping

def load_ARG_length():
    ARG_length = dict()
    for lines in open(args.d.replace('.structure.txt','.db.fasta.length.txt'),'r'):
        lines_set = lines.split('\n')[0].split('\t')
        ARG,ARG2,genelength = lines_set[0:3]
        ARG_length.setdefault(ARG,int(genelength)*3) # AA to DNA
    return ARG_length

def sum_ARG(allsearchoutput):
    ARG_mapping = load_ARG_mapping()
    ARG_length = load_ARG_length()
    allsampleARG = []
    allsamplename = []
    for searchoutput in allsearchoutput:
        sampleARG = dict()
        # load kraken results
        kraken = glob.glob(searchoutput.replace('.blast.txt.filter','.kraken.kreport'))
        copy_16S = 1
        gene_length = 1550 # 16S
        if kraken!= []:
            # metagenomes
            for lines in open(kraken[0],'r'):
                lines_set = lines.split('\n')[0].split('\t')
                if 'Bacteria' in lines:
                    copy_16S = float(lines_set[1])*2 # pair end
                    break
        if args.kkdbtype != '16S':
            # load average genome size
            AGS_result = glob.glob(searchoutput.replace('.blast.txt.filter', '.AGS.txt'))[0]
            for lines in open(AGS_result,'r'):
                if lines.startswith('average_genome_size'):
                    lines_set = lines.split('\n')[0].split('\t')
                    gene_length = float(lines_set[1])
        # load ARG blast results
        samplename = os.path.split(searchoutput)[-1].split('.blast.txt.filter')[0]
        allsamplename.append(samplename)
        for lines in open(searchoutput,'r'):
            lines_set = lines.split('\n')[0].split('\t')
            ARG = lines_set[1]
            sampleARG.setdefault(ARG,0)
            if copy_16S == 1: # genomes
                sampleARG[ARG] += 1.0 / copy_16S
            else: # metagenomes
                sampleARG[ARG] += 1.0 / copy_16S / ARG_length[ARG] * gene_length  # normalize against ARG length and 16S/genome length
        allsampleARG.append(sampleARG)
    # sum ARG blast results
    alloutput = []
    alloutput.append('\tARG_gene_family\tARG_antibiotics\t%s'%('\t'.join(allsamplename)))
    for ARG in ARG_mapping:
        temp_line = ARG_mapping[ARG]
        for sampleARG in allsampleARG:
            if ARG in sampleARG:
                temp_line += '\t%s'%(sampleARG[ARG])
            else:
                temp_line += '\t0'
        alloutput.append(temp_line)
    f1 = open('%s/../Sample_ARGpresence.txt'%(args.i),'w')
    f1.write('\n'.join(alloutput)+'\n')
    f1.close()

################################################### Programme #######################################################
sum_ARG(allsearchoutput)
