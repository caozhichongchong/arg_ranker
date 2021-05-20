import argparse
import os
from Bio import SeqIO


############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-db",
                    help="file name of your input database",
                     type=str, default='Butyrate.pro.aa',metavar='database.aa')
parser.add_argument("-dbf",
                    help="sequence format of your input database\
                    (1: nucleotide; 2: protein), \
                    (default \'1\' for nucleotide)",
                    metavar="1 or 2",
                    choices=[1, 2],
                    action='store', default=1, type=int)
parser.add_argument("-i",
                    help="input dir", type=str, default='.',metavar='current dir (.)')
parser.add_argument("-f",
                    help="input filename", type=str, default='.blast.txt',metavar='.blast.txt')
parser.add_argument('-s',
                    help="set the method to search the your database \
                    (1: blast; 2: hmm), \
                    (default \'1\' for blast search)",
                    metavar="1 or 2",
                    choices=[1, 2],
                    action='store', default=1, type=int)
# optional parameters
parser.add_argument('--g',
                    help="input genomes (T) or metagenomes (F) \
                    (default \'F\' for metagenomes)",
                    metavar="1 or 2",
                    choices=['T', 'F'],
                    action='store', default='F', type=str)
# optional search parameters
parser.add_argument('--id',
                    default=80.0, action='store', type=float, metavar='60.0',
                    help='Optional: set the amno acid based identity cutoff for blast (default is 60.0)\n'
                         'Leave it alone if hmm is used')
parser.add_argument('--ht',
                    default=80.0, action='store', type=float, metavar='60.0',
                    help='Optional: set the amno acid based hit-length cutoff for blast (default is 60.0)\n'
                         'Leave it alone if hmm is used')
parser.add_argument('--e',
                    default=1e-5, action='store', type=float, metavar='1e-5',
                    help='Optional: set the evalue cutoff for blast or hmm (default is 1e-5)')

################################################## Definition ########################################################
args = parser.parse_args()


################################################### Function #######################################################
def Calculate_length(file_name):
    DB_length=dict()
    try:
        for lines in open(file_name + '.length', 'r'):
            DB_length.setdefault(str(lines.split('\t')[0]), len(str(lines.split('\t')[-1]).replace('\n','')))
    except (IOError,FileNotFoundError):
        Fasta_name = open(file_name, 'r')
        f = open(file_name + '.length', 'w')
        for record in SeqIO.parse(Fasta_name, 'fasta'):
            f.write(str(record.id) + '\t'  + str(
                len(record.seq)) + '\n')
            DB_length.setdefault(str(record.id), len(str(record.seq)))
        f.close()
    return DB_length


def compare_blast(line1,line2):
    ID1=float(str(line1).split('\t')[2])
    ID2=float(str(line2).split('\t')[2])
    Loci1=(float(str(line1).split('\t')[6])+float(str(line1).split('\t')[7]))/2.0
    Loci2=(float(str(line2).split('\t')[6])+float(str(line2).split('\t')[7]))/2.0
    Length1 = abs(float(str(line1).split('\t')[6]) - float(str(line1).split('\t')[7]))
    Length2 = abs(float(str(line2).split('\t')[6]) - float(str(line2).split('\t')[7]))
    if abs(Loci1 - Loci2) > 500:
        return [line1,line2]
    elif ID1 > ID2:
        return [line1]
    elif ID2 > ID1:
        return [line2]
    elif ID1 == ID2:
        # same identity, choose the longer one
        if Length1 >= Length2:
            return [line1]
        else:
            return [line2]


def blast_list(file, Cutoff_identity,Cutoff_hitlength):
    f1=open(file+'.filter','w')
    All_hit = dict()
    if args.dbf == 1:
        ratio = 1.0
    else:
        ratio = 3.0
    for line in open(file, 'r'):
        if float(str(line).split('\t')[2]) >= Cutoff_identity:
            try:
                if args.g == 'T':
                    length_cut = float(str(line).split('\t')[3]) >= Cutoff_hitlength * float(
                        DB_length.get(str(line).split('\t')[1],0.0)) / 100
                else:
                    length_cut = float(str(line).split('\t')[3]) >= Cutoff_hitlength * 100.0/ratio/100.0
                if length_cut:
                    if str(line).split('\t')[0] not in All_hit:
                        All_hit.setdefault(str(line).split('\t')[0],
                        [line])
                    else:
                        All_hit[str(line).split('\t')[0]] = compare_blast(All_hit[str(line).split('\t')[0]][0],
                        line)
            except TypeError:
                print (str(line).split('\t')[1])
    for genes in All_hit:
        for hits in All_hit[genes]:
            f1.write(hits)
    f1.close()


################################################### Programme #######################################################
if args.s == 1:
    DB_length = dict()
    if args.ht > 0.0:
        DB_length = Calculate_length(args.db)
    blast_list(os.path.join(args.i, args.f), float(args.id), float(args.ht))
else:
    pass
