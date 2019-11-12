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
                    help="prefix name for input file", type=str, default='',metavar='.usearch.txt.aa')

parser.add_argument("-n",
                    help="prefix name for usearch result", type=str, default='.usearch.txt',metavar='.usearch.txt')
parser.add_argument("-r",
                    help="output dir", type=str, default='.',metavar='current dir (.)')


################################################## Definition ########################################################
args = parser.parse_args()


################################################### Function #######################################################
def Extractaa(root, searchfile, orffile, resultdir):
    # extract the query aa sequences according to a usearch or diamond alignment output
    # generate a smaller data of potential intI1 or sul1 for blastp search
    # input the query ORF sequences
    AA_seq = dict()
    try:
        for record in SeqIO.parse(open(os.path.join(root, orffile), 'r'), 'fasta'):
            AA_seq.setdefault(str(record.id), str(record.seq))
        f1 = open(os.path.join(resultdir, searchfile + '.aa'), 'w')
        try:
            for line in open(os.path.join(resultdir, searchfile), 'r'):
                try:
                    AA = str(line).split('\t')[0].split(' ')[0]
                    if AA_seq[AA] != '':
                        # avoid duplicate ORF
                        f1.write('>' + AA + '\n' +
                                 str(AA_seq[AA]) + '\n')
                        AA_seq[AA]=''
                except KeyError:
                    print ('AA not found for ' + AA)
                    flog.write('AA not found for ' + AA + '\n')
        except IOError:
            pass
        f1.close()
    except (IOError):
        print ('Files were missing: ' + orffile)


################################################### Programme #######################################################
if args.ni != '':
    input_file = args.f.split(args.ni)[0]+args.n
else:
    input_file = args.f+ args.n
Extractaa(args.i, input_file,
          args.f, args.r)