import os
from Bio import SeqIO
import argparse
import glob
import statistics
import traits_finder
import sys

################################################### Decalration #######################################################
print ("\
------------------------------------------------------------------------\n\
traits_finder searches and summarizes traits in genomes and metagenomes\n\
input: reference database and folder of genomes/metagenomes\n\
requirement: blast \n\n\
optional: diamond, bwa, hs-blastn, usearch \n\n\
Copyright:An Ni Zhang, Prof. Eric Alm, MIT\n\n\
Citation:\n\
Contact anniz44@mit.edu\n\
------------------------------------------------------------------------\n\
")

def main():
    usage = ("usage: traits_finder -t your.otu.table -s your.otu.seqs")
    version_string = 'traits_finder {v}, on Python {pyv[0]}.{pyv[1]}.{pyv[2]}'.format(
        v=traits_finder.__version__,
        pyv=sys.version_info,
    )
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
    parser.add_argument("-m",
                        help="mapping file of traits to function", type=str,
                        default='Butyrate.pro.mapping.txt',
                        metavar='Butyrate.pro.mapping.txt')
    parser.add_argument("-i",
                        help="input folder of metagenomes or genomes", type=str,
                        default='.',metavar='current dir (.)')
    parser.add_argument("-inf",
                        help="input format\
                        (1: genomes; 2: metagenomes), \
                        (default \'1\' for genomes)",
                        metavar="1 or 2",
                        choices=[1, 2],
                        action='store', default=1, type=int)
    parser.add_argument("-fa",
                        help="input format of genome/metagenome sequence",
                        type=str, default='.fa', metavar='.fasta, .fna, .fastq or .fa')
    parser.add_argument('-s',
                        help="set the method to search the your database \
                        (1: blast; 2: hmm), \
                        (default \'1\' for blast search)",
                        metavar="1 or 2",
                        choices=[1, 2],
                        action='store', default=1, type=int)
    # optional parameters
    parser.add_argument("--l",
                        help="input list of metagenomes", type=str,
                        default='None',metavar='rep_metagenomes.txt')
    parser.add_argument("--meta",
                        help="metadata  of metagenomes", type=str,
                        default='None',
                        metavar='metadata.metagenomes.txt')
    parser.add_argument("--orf",
                        help="input format of genomes orfs", type=str, default='.genes.faa',metavar='.faa')
    # optional output setup
    parser.add_argument("--r",
                        help="output directory or folder of your results",
                        type=str, default='Result',metavar='Result')
    parser.add_argument("--r16",
                        help="output directory or folder of your 16S sequences",
                        type=str, default='Result',metavar='Result')
    # optional search parameters
    parser.add_argument('--t',
                        help="Optional: set the thread number assigned for running XXX (default 1)",
                        metavar="1 or more", action='store', default=1, type=int)
    parser.add_argument('--id',
                        default=75.0, action='store', type=float, metavar='60.0',
                        help='Optional: set the amno acid based identity cutoff for blast (default is 80.0)\n'
                             'Leave it alone if hmm is used')
    parser.add_argument('--ht',
                        default=75.0, action='store', type=float, metavar='60.0',
                        help='Optional: set the amno acid based hit-length cutoff for blast (default is 80.0)\n'
                             'Leave it alone if hmm is used')
    parser.add_argument('--e',
                        default=1e-2, action='store', type=float, metavar='1e-5',
                        help='Optional: set the evalue cutoff for blast or hmm (default is 1e-5)')
    # requirement for software calling
    parser.add_argument('--u',
                        help="Optional: use two-step method for blast search,"+
                             " \'None\' for using one step, \'usearch\' or \'diamond\' for using two-step \
                             (complete path to usearch or diamond if not in PATH, \
                             please make sure the search tools can be directly called), (default: \'None\')",
                        metavar="None or usearch",
                        action='store', default='None', type=str)
    parser.add_argument('--hmm',
                        help="Optional: complete path to hmmscan if not in PATH,",
                        metavar="/usr/local/bin/hmmscan",
                        action='store', default='hmmscan', type=str)
    parser.add_argument('--bp',
                        help="Optional: complete path to blastp or blastn if not in PATH, \'None\' for no blast search",
                        metavar="/usr/local/bin/blastp",
                        action='store', default='blastp', type=str)
    parser.add_argument('--bwa',
                        help="Optional: complete path to bwa if not in PATH,",
                        metavar="/usr/local/bin/bwa",
                        action='store', default='None', type=str)
    ################################################## Definition ########################################################
    args = parser.parse_args()
    workingdir=os.path.abspath(os.path.dirname(__file__))
    ################################################### Programme #######################################################
    f1 = open ('traits_finder.log','w')
    if args.inf == 1:
        cmd = ('python '+workingdir+'/Traits_WG.py -db %s -dbf %s -i %s -s %s --fa %s --orf %s --r %s --r16 %s --t %s --id %s --ht %s --e %s --u %s --hmm %s --bp %s --bwa %s\n'
        % (str(args.db),str(args.dbf),str(args.i),str(args.s),str(args.fa),str(args.orf),str(args.r),str(args.r16),str(args.t),str(args.id),str(args.ht),str(args.e),str(args.u),str(args.hmm),str(args.bp),str(args.bwa)))
        #cmd += ('python '+workingdir+'/scripts/Traits_summary_WG.py -t %s -db %s --fa %s --orf %s -i %s -m %s --r %s --r16 %s --s %s -c %s\n'
        #%(str(os.path.split(args.db)[1]),str(args.db),str(args.fa),str(args.orf),str(args.i),str(args.m),str(args.r),str(args.r16),str(os.path.join(args.r,'summary')),str(args.id)))
        f1.write(cmd)
        os.system(cmd)
    else:
        cmd = ('python '+workingdir+'/Traits_MG.py -db %s -dbf %s -i %s -s %s --fa %s -l %s --r %s --r16 %s --t %s --id %s --ht %s --e %s --u %s --hmm %s --bp %s --bwa %s\n'
        % (str(args.db),str(args.dbf),str(args.i),str(args.s),str(args.fa),str(args.l),str(args.r),str(args.r16),str(args.t),str(args.id),str(args.ht),str(args.e),str(args.u),str(args.hmm),str(args.bp),str(args.bwa)))
        #if args.s == 1:
        #   cmd += ('python '+workingdir+'/Copynumber_traits.py -i %s -i16 %s -o %s -f %s.blast.txt -f16 %s.16S.txt -l %s.length -tf %s -m %s -c %s'
        #    % (str(os.path.join(args.r,'search_output/0'))),str(args.r16),str(os.path.join(args.r,'summary'))),str(args.fa),str(args.fa),str(args.db),str(args.m),str(args.meta),str(args.id))
        #else:
            #cmd += ('python '+workingdir+'/Copynumber_traits_hmm.py -i %s -i16 %s -o %s -f %s.hmm2.txt -f16 %s.16S.txt -l %s.length -tf %s -m %s -c %s'
            #% (str(os.path.join(args.r,'search_output/0'))),str(args.r16),str(os.path.join(args.r,'summary'))),str(args.fa),str(args.fa),str(args.db),str(args.m),str(args.meta),str(args.e))
        #cmd += ('python '+workingdir+'/scripts/Traits_summary_MG.py -t %s -db %s --fa %s --orf %s -i %s -m %s --r %s --r16 %s --s %s -c %s\n'
        #%(str(os.path.split(args.db)[1]),str(args.db),str(args.fa),str(args.orf),str(args.i),str(args.m),str(args.r),str(args.r16),str(os.path.join(args.r,'summary')),str(args.id)))
        f1.write(cmd)
        os.system(cmd)

################################################## Function ########################################################

if __name__ == '__main__':
    main()
