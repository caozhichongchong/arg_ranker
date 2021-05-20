import argparse
import os
from Bio import SeqIO


############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                    help="file name of your genome", type=str, default='.fa',metavar='.fa')


################################################## Definition ########################################################
args = parser.parse_args()


################################################### Function ########################################################
def check_genome(filename):
    for record in SeqIO.parse(open(filename, 'r'), 'fasta'):
        if 'GCA_' in record.id:
            try:
                str(record.id).split('_')[2]
                return 0
            except IndexError:
                return 1
        else:
            return 0


def formate_genome(filename,Necessary):
    if Necessary == 1:
        os.system('mv '+filename +' '+filename+'.unformat')
        f1=open(filename,'w')
        for record in SeqIO.parse(open(filename+'.unformat', 'r'), 'fasta'):
            newid=str(record.id)+'_'+str(record.description).split(' ')[0]
            f1.write('>' + str(newid) + '\n' + str(record.seq) + '\n')
        f1.close()


################################################### Programme #######################################################
formate_genome(args.i,check_genome(args.i))
