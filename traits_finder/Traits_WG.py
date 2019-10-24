import os
from Bio import SeqIO
import argparse
import glob


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
                    help="input dir of WGD", type=str, default='.',metavar='current dir (.)')
parser.add_argument('-m',
                    help="set the model to predict the results \
                    (1: simple; 2: phylogenetic; 3: advanced), \
                    (default \'1\' for simple)",
                    metavar="1 or 2 or 3",
                    choices=[1, 2, 3],
                    action='store', default=1, type=int)
parser.add_argument('-s',
                    help="set the method to search the your database \
                    (1: blast; 2: hmm), \
                    (default \'1\' for blast search)",
                    metavar="1 or 2",
                    choices=[1, 2],
                    action='store', default=1, type=int)
# optional parameters
parser.add_argument("--fa",
                    help="input format of genome sequence", type=str,
                    default='.fna.add',metavar='.fasta, .fna or .fa')
parser.add_argument("--orf",
                    help="input format of genome orfs", type=str,
                    default='.genes.faa',metavar='.faa')
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
try:
    os.mkdir(args.r)
except OSError:
    pass
try:
    os.mkdir('subscripts')
except OSError:
    pass
workingdir=os.path.abspath(os.path.dirname(__file__))

################################################### Programme #######################################################
# set the search for database type
diamond_set = args.u
if args.dbf == 1:
    blast_set=args.bp.replace('blastp','blastn')
    if 'diamond' in args.u:
        print('Diamond cannot use dna database\nDirect use blast!')
        diamond_set = 'None'
else:
    blast_set=args.bp.replace('blastn','blastp')
    diamond_set=args.u.replace('diamond','diamond-blastp')
# search the database in all genomes
cmds = 'python '+ workingdir +'/scripts/Search.WG.py -i ' + args.i  +\
                ' -db ' + args.db + ' -dbf ' + str(args.dbf) + ' -s ' + str(args.s) + ' --r ' + str(args.r) + ' --t ' + str(args.t) + \
                ' --u ' + str(diamond_set) + ' --hmm ' + str(args.hmm) + ' --bp ' + str(blast_set) + \
                ' --ht ' + str(args.ht) + ' --id ' + str(args.id) + ' --fa ' + str(args.fa) + \
                ' --e ' + str(args.e) + ' --orf ' + str(args.orf) + ' --r16 ' + str(args.r16)+\
        ' --bwa ' + str(args.bwa) + ' \n'
#print(cmds)
os.system(cmds)

# run all bash
list_of_files = glob.glob('subscripts/*.sh')
f1 = open("all.sh", 'w')
f1.write("#!/bin/bash \nmodule add c3ddb/blast+/2.7.1 \n")
# make all database
if args.s == 1:
    if args.bwa != 'None':
        try:
            ftest=open(args.db+'.bwt','r')
        except IOError:
            f1.write('bwa index %s \n' % (args.db))
    if args.dbf == 1:
        try:
            ftest = open(args.db + '.nin', 'r')
        except IOError:
            f1.write(os.path.join(os.path.split(args.bp)[0],'makeblastdb -in %s -dbtype nucl\n' % (args.db)))
            if args.u != 'None':
                if 'usearch' in args.u:
                    f1.write(os.path.join(os.path.split(args.u)[0], 'usearch -makeudb_usearch %s -output %s.udb\n' % (args.db,args.db)))
                elif 'hs-blastn' in args.u:
                    f1.write(os.path.join(os.path.split(args.u)[0],
                                          'windowmasker -in %s -infmt blastdb -mk_counts -out %s.counts\n' % (args.db, args.db)))
                    f1.write(os.path.join(os.path.split(args.u)[0],
                                          'windowmasker -in %s.counts -sformat obinary -out %s.counts.obinary -convert\n' % (
                                          args.db, args.db)))
                    f1.write(os.path.join(os.path.split(args.u)[0],
                                          'hs-blastn index %s\n' % (
                                          args.db)))
    else:
        try:
            ftest = open(args.db + '.pin', 'r')
        except IOError:
            f1.write(os.path.join(os.path.split(args.bp)[0],'makeblastdb -in %s -dbtype prot\n' % (args.db)))
            if args.u != 'None':
                if 'usearch' in args.u:
                    f1.write(os.path.join(os.path.split(args.u)[0], 'usearch -makeudb_usearch %s -output %s.udb\n' % (args.db,args.db)))
                elif 'diamond' in args.u:
                    f1.write(os.path.join(os.path.split(args.u)[0],
                                          'diamond makedb --in %s -d %s.dmnd\n' % (args.db, args.db)))
# makedb for 16S
try:
    ftest = open(workingdir +'/database/85_otus.fasta.udb', 'r')
except IOError:
    if args.u != 'None':
        if 'usearch' in args.u:
            f1.write(os.path.join(os.path.split(args.u)[0],
                                  'usearch -makeudb_usearch %s -output %s.udb\n' % (workingdir +'/database/85_otus.fasta', workingdir +'/database/85_otus.fasta')))
try:
    ftest = open(workingdir +'/database/85_otus.fasta.counts', 'r')
except IOError:
    if args.u != 'None':
            if 'hs-blastn' in args.u:
                f1.write(os.path.join(os.path.split(args.u)[0],
                                      'windowmasker -in %s -infmt blastdb -mk_counts -out %s.counts\n' % (
                                          workingdir +'/database/85_otus.fasta', workingdir +'/database/85_otus.fasta')))
                f1.write(os.path.join(os.path.split(args.u)[0],
                                      'windowmasker -in %s.counts -sformat obinary -out %s.counts.obinary -convert\n' % (
                                          workingdir +'/database/85_otus.fasta', workingdir +'/database/85_otus.fasta')))
                f1.write(os.path.join(os.path.split(args.u)[0],
                                      'hs-blastn index %s\n' % (
                                          workingdir +'/database/85_otus.fasta')))
# bash for all subscripts
for file_name in list_of_files:
    if all(keys not in file_name for keys in ['SearchWG.sh','all.sh']):
        f1.write("sbatch -p defq -c 40 -t 5-00:00:00 --mem=100000 -J " + str(file_name) + "_traits -o " + str(file_name) + ".out -e " + str(file_name) + ".err " + str(file_name) +" \n")
    #f1.write("nohup sh " + str(file_name) + '  & \n')
f1.close()
#os.system("nohup sh all.sh & \n")
