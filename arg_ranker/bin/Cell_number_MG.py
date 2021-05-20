import os
from Bio import SeqIO
import argparse
import glob


############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-m",
                    help="file name of your cell number mapping file", type=str,
                    default='all_KO30_name.list',
                    metavar='all_KO30_name.list')
parser.add_argument("-r",
                    help="output directory or folder of your results",
                    type=str, default='Result_traits',metavar='Result_traits')

################################################## Definition ########################################################
args = parser.parse_args()
try:
    os.mkdir(os.path.join(args.r))
except OSError:
    pass


################################################### Function ########################################################
def cell_calculate(inputfile,mapping_list):
    output = dict()
    for lines in open(inputfile):
        gene = lines.split('\t')[1]
        maplength = float(lines.split('\t')[3])
        if gene not in mapping_list:
            print(gene,' not in mapping list!')
        else:
            geneID = mapping_list[gene][0]
            genelength = mapping_list[gene][1]
            if geneID not in output:
                output.setdefault(geneID,maplength/genelength)
            else:
                output[geneID] += maplength/genelength
    average_output = 0
    for geneID in output:
        average_output += output[geneID]
    return average_output/float(len(output))


################################################### Programme #######################################################
# checkfiles
try:
    ftry = open(os.path.join(os.path.join(args.r),'cell.copynum.all.txt'), 'r')
except (IOError,FileNotFoundError):
    # load all files
    # should load 16S
    inputfiles=glob.glob(os.path.join(args.r16,'*/*.ESSG.diamond.txt'))
    # load cell number mapping file
    # gene ID: genotype ID, length
    Mapping = dict()
    for lines in open(args.m, 'r'):
        Mapping.setdefault(lines.split('\t')[0],
                           [lines.split('\t')[1],float(lines.split('\t')[2].split('\r')[0].split('\n')[0])])
    # process cell num
    Cellnum=dict()
    for filenames in inputfiles:
        filename = os.path.split(filenames)[-1].split('_1.')[0].split('_2.')[0].split('.fastq')[0].split('.fq')[0].split('.ESSG.diamond.txt')[0]
        if Cell_out == 0:
            if filename not in Cellnum:
                Cellnum.setdefault(filename,[])
            Cellnum[filename].append(cell_calculate(filenames,Mapping))

    # output results
    if Cell_out == 0:
        f1=open(os.path.join(os.path.join(args.r),'cell.copynum.all.txt'),'w')
        for filename in Cellnum:
            f1.write(filename)
            for number in Cellnum[filename]:
                f1.write('\t'+str(number))
            f1.write('\n')
        f1.close()