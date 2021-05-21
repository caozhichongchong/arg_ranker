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

def sum_ARG(allsearchoutput):
    ARG_mapping = load_ARG_mapping()
    allsampleARG = []
    allsamplename = []
    for searchoutput in allsearchoutput:
        sampleARG = dict()
        # load kraken results
        kraken = glob.glob(searchoutput.replace('.blast.txt.filter','.kraken.kreport'))
        copy_16S = 1
        if kraken!= []:
            # metagenomes
            for lines in open(kraken[0],'r'):
                lines_set = lines.split('\n')[0].split('\t')
                if 'Bacteria' in lines:
                    copy_16S = float(lines_set[1])*2 # pair end
                    break
        # load ARG blast results
        samplename = os.path.split(searchoutput)[-1].split('.blast.txt.filter')[0]
        allsamplename.append(samplename)
        for lines in open(searchoutput,'r'):
            lines_set = lines.split('\n')[0].split('\t')
            ARG = lines_set[1]
            sampleARG.setdefault(ARG,0)
            sampleARG[ARG] += 1.0/copy_16S
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
