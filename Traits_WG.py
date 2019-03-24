import os
from Bio import SeqIO
import argparse
import glob


############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-db",
                    help="file name of your input database", type=str, default='Butyrate.pro.aa',metavar='database.aa')
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
                    help="input format of genome sequence", type=str, default='.fna.add',metavar='.fasta, .fna or .fa')
parser.add_argument("--orf",
                    help="input format of genome orfs", type=str, default='.genes.faa',metavar='.faa')
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
                    default=1e-5, action='store', type=float, metavar='1e-5',
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
                    help="Optional: complete path to blastp if not in PATH,",
                    metavar="/usr/local/bin/blastp",
                    action='store', default='blastp', type=str)

################################################## Definition ########################################################
args = parser.parse_args()
try:
    os.mkdir(args.r)
except OSError:
    pass


################################################### Programme #######################################################
# search the database in all genomes
cmds = 'python scripts/Search.py -i ' + args.i  +\
                ' -db ' + args.db + ' -s ' + str(args.s) + ' --r ' + str(args.r) + ' --t ' + str(args.t) + \
                ' --u ' + str(args.u) + ' --hmm ' + str(args.hmm) + ' --bp ' + str(args.bp) + \
                ' --ht ' + str(args.ht) + ' --id ' + str(args.id) + ' --fa ' + str(args.fa) + \
                ' --e ' + str(args.e) + ' --orf ' + str(args.orf) + ' --r16 ' + str(args.r16)+ ' \n'
os.system(cmds)

# run all bash
list_of_files = glob.glob('*.sh')
f1 = open("all.sh", 'a')
f1.write("#!/bin/bash \nmodule add c3ddb/blast+/2.7.1 \n")
#f1.write('#!/bin/bash\nexport PATH=/scratch/users/anniz44/bin/miniconda3/bin:$PATH\n'+\
#         'export PATH=/scratch/users/anniz44/bin/miniconda3/bin/bin:$PATH\n'+\
#         'cd /scratch/users/anniz44/bin/miniconda3/bin/traits_search/\n')
for file_name in list_of_files:
    f1.write("sbatch -p defq -c 30 -t 1-00:00:00 --mem=100000 -J traits -o " + str(file_name) + ".out -e " + str(file_name) + ".err /scratch/users/anniz44/scripts/traits_search/" + str(file_name) +" \n")
    #f1.write("nohup sh " + str(file_name) + '  & \n')
f1.close()
#os.system("nohup sh all.sh & \n")
