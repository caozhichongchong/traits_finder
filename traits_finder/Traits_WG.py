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
                    (1: blast; 2: hmm; 3: alignment), \
                    (default \'1\' for blast search)",
                    metavar="1 or 2",
                    choices=[1, 2, 3],
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
parser.add_argument('--u','--usearch',
                        help="Optional: use two-step method for blast search,"+
                             " \'None\' for using one step, \'usearch\' for using two-step \
                             (complete path to usearch if not in PATH), (default: \'None\')",
                        metavar="None or usearch",
                        action='store', default='None', type=str)
parser.add_argument('--dm', '--diamond',
                      help="Optional: use two-step method for blast search," +
                           " \'None\' for using one step, \'diamond\' for using two-step \
                           (complete path to diamond if not in PATH), (default: \'None\')",
                      metavar="None or diamond",
                      action='store', default='None', type=str)
parser.add_argument('--hs',
                      help="Optional: use two-step method for blast search," +
                           " \'None\' for using one step, \'hs-blastn\' for using two-step \
                           (complete path to hs-blastn if not in PATH), (default: \'None\')",
                      metavar="None or hs-blastn",
                      action='store', default='None', type=str)
parser.add_argument('--hmm',
                    help="Optional: complete path to hmmscan if not in PATH,",
                    metavar="/usr/local/bin/hmmscan",
                    action='store', default='hmmscan', type=str)
parser.add_argument('--bp',
                    help="Optional: complete path to blast if not in PATH, \'None\' for no blast search",
                    metavar="/usr/local/bin/blast",
                    action='store', default='blast', type=str)
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
# search the database in all genomes
cmds = 'python '+ workingdir +'/scripts/Search.WG.py -i ' + args.i  +\
                ' -db ' + args.db + ' -dbf ' + str(args.dbf) + ' -s ' + str(args.s) + ' --r ' + str(args.r) + ' --t ' + str(args.t) + \
                ' --dm ' + str(args.dm) + ' --u ' + str(args.u) + ' --hs ' + str(args.hs) + ' --hmm ' + str(args.hmm) + ' --bp ' + str(args.bp) + \
                ' --ht ' + str(args.ht) + ' --id ' + str(args.id) + ' --fa ' + str(args.fa) + \
                ' --e ' + str(args.e) + ' --orf ' + str(args.orf) + ' --r16 ' + str(args.r16)+\
        ' --bwa ' + str(args.bwa) + ' \n'
#print(cmds)
os.system(cmds)

# run all bash
list_of_files = glob.glob('subscripts/*.sh')
f1 = open("all.sh", 'w')
f1.write("#!/bin/bash \nmodule add c3ddb/blast+/2.7.1 \n")

# bash for all subscripts
for file_name in list_of_files:
    if all(keys not in file_name for keys in ['SearchWG.sh','all.sh']):
        f1.write("sbatch -p defq,sched_mem1TB -c 40 -t 3-00:00:00 --mem=100000 -J " + os.path.split(str(file_name))[-1] + "_"+ \
                 str(os.path.split(args.db)[1]) + " -o " + str(file_name) + ".out -e " + str(file_name) + ".err " + str(file_name) +" \n")
    #f1.write("nohup sh " + str(file_name) + '  & \n')
f1.close()
#os.system("nohup sh all.sh & \n")
