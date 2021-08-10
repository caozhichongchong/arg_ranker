import argparse
import os
import arg_ranker
import sys
import pandas as pd
import glob

################################################### Decalration #######################################################
print ("\
------------------------------------------------------------------------\n\
arg_ranker evaluates the risk of ARGs in genomes/metagenomes\n\
Requirement: blast, diamond, kraken, python 3\n\
Optimal: a table of the metadata of your samples (see example/All_sample_metadata.txt).\n\
Copyright:An-Ni Zhang, MIT; Prof. Tong Zhang, University of Hong Kong\n\
Contact caozhichongchong@gmail.com\n\
------------------------------------------------------------------------\n\
")

def main():
    usage = ("usage: arg_ranker -i metagenomes/ -dm diamond")
    version_string = 'arg_ranker {v}, on Python {pyv[0]}.{pyv[1]}.{pyv[2]}'.format(
        v=arg_ranker.__version__,
        pyv=sys.version_info,
    )
    ############################################ Arguments and declarations ##############################################
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i',
                        default="example",
                        action='store', type=str,
                        metavar='metagenomes/',
                        help="input directory or metagenomes/genomes (.fa, .fastq, .fq, or .fastq)")
    optional.add_argument('-m',
                        default="None", action='store', type=str, metavar='metadata.txt',
                        help="input table containing the metadata of your samples\n"+
                             "\'None\' for no metadata input")
    optional.add_argument('-o',
                        default="arg_ranking", action='store', type=str, metavar='arg_ranking',
                        help="output directory to store results")
    optional.add_argument('-t',
                        help="Optional: set the thread number assigned for running arg_rankers (default 1)",
                        metavar="1 or more", action='store', default=1, type=int)
    optional.add_argument('-dm',
                          help="Optional: complete path to diamond if not in PATH",
                          metavar="/usr/local/bin/diamond",
                          action='store', default='diamond', type=str)
    optional.add_argument('-bl',
                          help="Optional: complete path to blast if not in PATH",
                          metavar="/usr/local/bin/blast",
                          action='store', default='blast', type=str)
    optional.add_argument('-mc',
                          help="Optional: complete path to MicrobeCensus if not in PATH",
                          metavar="/usr/local/bin/run_microbe_census.py",
                          action='store', default='run_microbe_census.py', type=str)
    optional.add_argument('-kk',
                          help="Optional: complete path to kraken2 if not in PATH",
                          metavar="/usr/local/bin/kraken2",
                          action='store', default='kraken2', type=str)
    optional.add_argument('-kkdb',
                          help="Optional: complete path to kraken2 database if not in PATH",
                          metavar="/usr/local/bin/kraken2/krakendb",
                          action='store', default='krakendb', type=str)
    optional.add_argument('-kkdbtype',
                          help="Optional: type of kraken2 database (default = standard)",
                          choice = ['standard','16S'],
                          metavar="standard or 16S",
                          action='store', default='standard', type=str)

    ################################################## Definition ########################################################
    args = parser.parse_args()
    workingdir = os.path.abspath(os.path.dirname(__file__))
    Metadata=args.m
    inputfasta = os.path.abspath(args.i)
    output_dir = os.path.abspath(args.o)
    if args.i == 'example':
        inputfasta = '%s/example/'%(workingdir)
        Metadata = '%s/example/metadata.txt'%(workingdir)
        try:
            f1 = open('%s/example/WEE300_all-trimmed-decont_1.fastq'%(workingdir),'r')
        except IOError:
            print('loading example files')
            os.system('unzip %s/example/WEE300_all-trimmed-decont_1.fastq.zip -d %s/example/'%(workingdir,workingdir))
    ARGranks=os.path.join(workingdir,
    'data', 'ARG_rank.txt')
    ARGdatabase = os.path.join(workingdir,
                            'data', 'SARG.db.fasta')
    ARGdatabase_mapping = os.path.join(workingdir,
                               'data', 'SARG.structure.txt')
    ESSGdatabase = os.path.join(workingdir,
                               'data', 'all_KO30.pro.fasta')
    try:
        os.mkdir(output_dir)
    except OSError:
        pass
    ################################################## Function ########################################################
    def Rank_num(rank, RKN, RK_profile):
        # calculate the nuber of genes in each rank
        try:
            RKN[RK_profile[rank]] += 1.0
        except (KeyError, IndexError, TypeError, ValueError):
            pass

    def Level_ranking(row, ARGlist, RK, RKN, fout, RK_profile, MD, inputfile):
        Abu = [0.0, 0.0, 0.0, 0.0, 0.0]
        ARGposition = 1
        # calculate the overall rank-based contribution by each rank
        for ARG in ARGlist:
            rank = RK.get(ARG, 'None')
            if rank == 'None':
                if ARG != 'Reference':
                    print ('ARGs in mothertable do not match with the ARGs in ARG_rank.txt.\nPlease check ' \
                           + ARG + ' in ' + inputfile + '!\n')
            else:
                try:
                    Abu[RK_profile[rank]] += float(row[ARGposition])
                except (KeyError, IndexError, TypeError, ValueError):
                    pass
            ARGposition += 1
        total_abu_all = sum(Abu)
        total_abu_I_IV = sum(Abu[0:4])
        samplename = row[0]
        if total_abu_all > 0:
            Abu1 = []
            Abu2 = []
            for abu in range(0, 5):
                Abu1.append(float(Abu[abu]) / total_abu_all)
            Num_all = sum(RKN[0:5])  # considering un-assessed ARGs
            if total_abu_I_IV > 0:
                for abu in range(0, 5):  # considering un-assessed ARGs
                    Abu2.append((Abu[abu] / RKN[abu]) / (total_abu_all / Num_all))  # considering un-assessed ARGs
            else:
                Abu2 = [0.0, 0.0, 0.0, 0.0, 0.0]  # considering un-assessed ARGs
            # output results
            if MD != {}:
                metadata_sample = MD.get(samplename, 'None')
                fout.write(str(samplename) + '\t' +
                           '\t'.join(str('%.1E' % abu) for abu in Abu1) +
                           '\t' + str('%.1E' % total_abu_all) + '\t' +
                           str('-'.join(str("%.1f" % abu) for abu in Abu2)) +
                           '\t' + '\t'.join(str('%.1f' % abu) for abu in Abu2) +
                           '\t' + str(metadata_sample) + '\n')
            else:
                fout.write(str(samplename) + '\t' +
                           '\t'.join(str('%.1E' % abu) for abu in Abu1) +
                           '\t' + str('%.1E' % total_abu_all) + '\t' +
                           str('-'.join(str("%.1f" % abu) for abu in Abu2)) +
                           '\t' + '\t'.join(str('%.1f' % abu) for abu in Abu2) + '\n')
        else:
            fout.write('%s\t0\t0\t0\t0\t0\t0\t0.0-0.0-0.0-0.0-0.0\t0\t0\t0\t0\t0\t0\n'%(samplename))

    def split_string_last(input_string, substring):
        return input_string[0: input_string.rfind(substring)]

    def makedatabase(search_method, db_file):
        # aa database
        if 'blast' in search_method:
            try:
                f1 = open("%s.phr" % (db_file), 'r')
            except IOError:
                os.system('%s -in %s -dbtype prot' %
                          (os.path.join(os.path.split(search_method)[0], 'makeblastdb'), db_file))
        if 'diamond' in search_method:
            try:
                f1 = open("%s.dmnd" % (db_file), 'r')
            except IOError:
                os.system('%sdiamond makedb --in %s -d %s.dmnd' %
                          (split_string_last(args.dm, 'diamond'), db_file, db_file))

    def setcutoff(sample):
        if '.fq' in sample or '.fastq' in sample:
            # metagenomes fastq
            return [80, 75, 1e-7, 'F']
        else:
            # genomes fasta
            return [90, 80, 1e-5, 'T']

    def searchARG(allsamples):
        cmds = '#!/bin/bash\n'
        #cmds += 'source ~/.bashrc\npy37\nexport LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n'
        for sample in allsamples:
            id_ARG, hit_ARG, evalue_ARG, sampletype = setcutoff(sample)
            samplename = os.path.split(sample)[-1]
            if sampletype == 'F':
                # run diamond for metagneomes
                sampleoutput = os.path.join(search_output, samplename + '.diamond.txt')
                try:
                    f1 = open(sampleoutput, 'r')
                except IOError:
                    cmds += "%sdiamond blastx --query %s --db %s.dmnd --out %s --outfmt 6 --max-target-seqs 1 --evalue %s --id %s --query-cover %s --threads %s \n" % (
                        split_string_last(args.dm, 'diamond'), sample, ARGdatabase, sampleoutput, evalue_ARG, id_ARG,
                        hit_ARG, args.t)
                    # extract candidate seqs
                    cmds += 'python %s/bin/Extract.MG.py -p 1 -i %s -f %s -n .diamond.txt -r %s\n' % (
                        workingdir, inputfasta, samplename, search_output)
                # compute taxonomy
                try:
                    sampleoutput = os.path.join(search_output, samplename + '.kraken.kreport')
                    f1 = open(sampleoutput, 'r')
                except IOError:
                    cmds += '%skraken2 --db %s %s --output %s --report %s.kreport --threads %s\n' % (
                        split_string_last(args.kk, 'kraken'), args.kkdb, sample, sampleoutput, sampleoutput, args.t)
                # compute average genome size
                if args.kkdbtype != '16S':
                    try:
                        sampleoutput = os.path.join(search_output, samplename + '.AGS.txt')
                        f1 = open(sampleoutput, 'r')
                    except IOError:
                        cmds += '%run_microbe_census.py -t %s %s %s \n' % (
                            split_string_last(args.mc, 'run_microbe_census'),args.t, sample, sampleoutput)
                sample = os.path.join(search_output, samplename + '.diamond.txt.aa')
            # run blast
            sampleoutput = os.path.join(search_output, samplename + '.blast.txt')
            try:
                f1 = open(sampleoutput, 'r')
            except IOError:
                cmds += "%sblastx -query %s -db %s -out %s -outfmt 6 -evalue %s -num_threads %s\n" % (
                    split_string_last(args.bl, 'blast'), sample, ARGdatabase, sampleoutput, evalue_ARG, args.t)
                # filter searching result
                cmds += 'python %s/bin/Filter.MG.py --g %s -i %s -f %s -db %s -dbf 2 -s 1 --ht %s --id %s --e %s \n' % (
                    workingdir, sampletype, search_output, sampleoutput, ARGdatabase, hit_ARG, id_ARG, evalue_ARG)
        return cmds

    def ranking_arg(inputfasta,Metadata):
        # set output file
        Mothertable = inputfasta
        fout = open(os.path.join(output_dir, 'Sample_ranking_results.txt'), 'w')
        # input metadata
        MD = dict()
        i = 0
        if Metadata != 'None':
            for lines in open(Metadata, 'r'):
                lines = str(lines).split('\r')[0].split('\n')[0]
                # output the lable in output file
                if i == 0:
                    fout.write('Sample\tRank_I_per\tRank_II_per\tRank_III_per\tRank_IV_per' +
                               '\tUnassessed_per\tTotal_abu\tRank_code\t' +
                               'Rank_I_risk\tRank_II_risk\tRank_III_risk\tRank_IV_risk' +
                               '\tUnassessed_risk\t%s\n' % '\t'.join(str(lines).split('\t')[1:]))
                else:
                    try:
                        # valid metadata input
                        MD.setdefault(str(lines).split('\t')[0], '\t'.join(str(lines).split('\t')[1:]))
                    except KeyError:
                        pass
                i += 1
        else:
            fout.write('Sample\tRank_I_per\tRank_II_per\tRank_III_per\tRank_IV_per' +
                       '\tUnassessed_per\tTotal_abu\tRank_code\t' +
                       'Rank_I_risk\tRank_II_risk\tRank_III_risk\tRank_IV_risk' +
                       '\tUnassessed_risk\n')
        # input ARG ranks
        RK = dict()
        RK_profile = {'I': 0, 'II': 1, 'III': 2, 'IV': 3,
                      'Unassessed': 4}
        RKN = [0.0, 0.0, 0.0, 0.0, 0.0]
        for lines in open(ARGranks, 'r'):
            lines = str(lines).split('\r')[0].split('\n')[0]
            RK.setdefault(lines.split('\t')[0], lines.split('\t')[-1])
            Rank_num(lines.split('\t')[-1], RKN, RK_profile)
        # input ARG mothertable
        # transpose the mothertable
        df = pd.read_csv(Mothertable, index_col=None, header=None, skipinitialspace=True, sep='\t')
        df.dropna(axis=0, thresh=2, subset=None, inplace=True)
        os.system('rm %s/Sample_ARGpresence.txt.t'%(output_dir))
        df = df.T
        df.to_csv(str(Mothertable) + '.t', header=None, index=None, sep='\t', mode='a')
        i = 0
        for row in df.itertuples(index=True, name='Pandas'):
            if i == 0:
                ARGlist = row[2:]
            elif i > 2:
                Level_ranking(row[1:], ARGlist, RK, RKN, fout, RK_profile, MD, Mothertable)
            i += 1
        fout.close()
        os.system('rm %s/Sample_ARGpresence.txt.t' % (output_dir))
        print('Finished ranking ARGs\nPlease check your results in ' +
              str(os.path.join(output_dir, 'Sample_ranking_results.txt')))

    ################################################### Programme #######################################################
    # run ARG ranking
    try:
        os.mkdir(output_dir)
    except OSError:
        pass
    if 'Sample_ARGpresence.txt' in inputfasta:
        ranking_arg(inputfasta,Metadata)
    # search ARGs in samples
    else:
        try:
            search_output = output_dir + '/search_output/'
            os.mkdir(search_output)
        except OSError:
            pass
        try:
            scripts_output = output_dir + '/script_output/'
            os.mkdir(scripts_output)
        except OSError:
            pass
        allsamples = glob.glob('%s/*.fa' % (inputfasta)) +\
                     glob.glob('%s/*.fasta' % (inputfasta)) + \
                     glob.glob('%s/*.fq' % (inputfasta)) + \
                     glob.glob('%s/*.fastq' % (inputfasta))
        print(allsamples,inputfasta,workingdir)
        # make database
        makedatabase(args.bl, ARGdatabase)
        makedatabase(args.dm, ARGdatabase)
        makedatabase(args.dm, ESSGdatabase)
        # search ARGs in all samples
        cmds = searchARG(allsamples)
        # output summary table
        cmds += 'python %s/bin/ARG_table.sum.py -i %s -d %s -kkdbtype %s\n' % (workingdir,search_output,ARGdatabase_mapping,args.kkdbtype)
        # risk ranking of ARGs
        cmds += 'arg_ranker -i %s/Sample_ARGpresence.txt -m %s -o %s\n'%(output_dir,Metadata,output_dir)
        f1 = open('%s/arg_ranker.sh'%(scripts_output),'w')
        f1.write(cmds)
        f1.close()
        print('Please run \nsh %s/arg_ranker.sh'%(scripts_output))

if __name__ == '__main__':
    main()
