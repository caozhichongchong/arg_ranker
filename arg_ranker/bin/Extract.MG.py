import argparse
from Bio import SeqIO
import os

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                    help="input dir", type=str, default='.',metavar='current dir (.)')
parser.add_argument("-f",
                    help="input filename", type=str, default='input.faa',metavar='input.faa')
parser.add_argument("-ni",
                    help="prefix name for input file", type=str, default='',metavar='.diamond.txt.aa')
parser.add_argument("-n",
                    help="prefix name for usearch result", type=str, default='.diamond.txt',metavar='.diamond.txt')
parser.add_argument("-r",
                    help="output dir", type=str, default='.',metavar='current dir (.)')
parser.add_argument("-d",
                    help="extra distance outside the target gene (default: 0 for 0bp)",
                    type=int, default=0,metavar='500 for 0.5kbp')
parser.add_argument("-p",
                    help="extract whole sequences or only the hit part \
                    (1: whole; 2: hit), \
                    (default \'1\' for blast search)",
                    metavar="1 or 2",
                    choices=[1, 2],
                    action='store', default=2, type=int)

################################################## Definition ########################################################
args = parser.parse_args()


################################################### Function #######################################################
def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N': 'N'}
    try:
        return ''.join([complement[base] for base in dna[::-1]])
    except KeyError:
        return dna


def Extractaa(root, searchfile, orffile, resultdir):
    # extract the query aa sequences according to a usearch or diamond alignment output
    # generate a smaller data of potential intI1 or sul1 for blastp search
    # input the query ORF sequences
    AA_seq = dict()
    try:
        f1 = open(os.path.join(resultdir, searchfile + '.aa'), 'w')
        if args.d > 0:
            f2 = open(os.path.join(resultdir, searchfile +
            '.extra' + str(args.d) + '.aa'), 'w')
        try:
            for line in open(os.path.join(resultdir, searchfile), 'r'):
                    AA = str(line).split('\t')[0].split(' ')[0]
                    loci1=int(str(line).split('\t')[6])
                    loci2=int(str(line).split('\t')[7])
                    if args.p == 2:
                        if AA not in AA_seq:
                            AA_seq.setdefault(AA,[[loci1,loci2]])
                        elif [loci1,loci2] not in AA_seq[AA]:
                            AA_seq[AA].append([loci1,loci2])
                    else:
                        AA_seq.setdefault(AA, ['whole'])
            tag = 'fasta'
            if orffile.split('.')[-1] in ['fastq','fq']:
                tag = 'fastq'
            for record in SeqIO.parse(open(os.path.join(root, orffile), 'r'), tag):
                AA = str(record.id)
                total_length = len(str(record.seq))
                if AA in AA_seq:
                    for locus in AA_seq[AA]:
                        if args.p == 2:
                            # extract hit sequences
                            loci1=int(locus[0])
                            loci2=int(locus[1])
                        else:
                            # extract whole sequences
                            loci1 = 1
                            loci2 = len(str(record.seq))
                        if loci1 < loci2:
                            # avoid duplicate ORF
                            f1.write('>%s_%s_%s\n' %(AA,str(loci1-1),str(loci2)) +
                                     str(record.seq)[(loci1-1):loci2] + '\n')
                            if args.d > 0:
                                f2.write('>%s_%s_%s\n' %(AA,str(max(loci1-1-args.d,0)),
                                                         str(min(loci2+args.d,total_length)))+
                                         str(record.seq)[max(loci1-1-args.d,0):
                                         min(loci2+args.d,total_length)] + '\n')
                        else:
                            f1.write('>%s_%s_%s\n' %(AA,str(loci2-1),str(loci1)) +
                                     str(reverse_complement(str(record.seq)[(loci2-1):loci1])) + '\n')
                            if args.d > 0:
                                f2.write('>%s_%s_%s\n' %(AA,str(max(loci2-1-args.d,0)),
                                                         str(min(loci1+args.d,total_length))) +
                                         str(reverse_complement(str(record.seq)[max(loci2-1-args.d,0):
                                         min(loci1+args.d,total_length)])) + '\n')
                    # finish extracting AA
                    AA_seq[AA]=''
        except (IOError,FileNotFoundError):
            pass
        f1.close()
    except (IOError,FileNotFoundError):
        print ('Files were missing: ' + orffile)


################################################### Programme #######################################################
if args.ni != '':
    input_file = args.f.split(args.ni)[0]+args.n
else:
    input_file = args.f+ args.n
Extractaa(args.i, input_file,
          args.f, args.r)
