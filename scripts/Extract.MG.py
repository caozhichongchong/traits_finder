import argparse
from Bio import SeqIO
import os

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                    help="input dir", type=str, default='.',metavar='current dir (.)')
parser.add_argument("-f",
                    help="input filename", type=str, default='input.faa',metavar='input.faa')
parser.add_argument("-n",
                    help="prefix name for usearch result", type=str, default='.usearch.txt',metavar='.usearch.txt')
parser.add_argument("-r",
                    help="output dir", type=str, default='.',metavar='current dir (.)')
parser.add_argument("-d",
                    help="extra distance outside the target gene (default: 0 for 0bp)",
                    type=int, default=0,metavar='2000 for 2kbp')


################################################## Definition ########################################################
args = parser.parse_args()


################################################### Function #######################################################
def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N': 'N'}
    return ''.join([complement[base] for base in dna[::-1]])


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
                    AA_seq.setdefault(AA,[loci1,loci2])
            tag = 'fasta'
            if 'fastq' in orffile:
                tag = 'fastq'
            for record in SeqIO.parse(open(os.path.join(root, orffile), 'r'), tag):
                AA = str(record.id)
                total_length = len(str(record.seq))
                if AA in AA_seq:
                    loci1=int(AA_seq[AA][0])
                    loci2=int(AA_seq[AA][1])
                    if loci1 < loci2:
                        # avoid duplicate ORF
                        f1.write('>' + AA + '\n' +
                                 str(record.seq)[(loci1-1):loci2] + '\n')
                        if args.d > 0:
                            f2.write('>' + AA + '\n' +
                                     str(record.seq)[max(loci1-1-args.d,0):
                                     min(loci2+args.d,total_length)] + '\n')
                    else:
                        f1.write('>' + AA + '\n' +
                                 str(reverse_complement(str(record.seq)[(loci2-1):loci1])) + '\n')
                        if args.d > 0:
                            f2.write('>' + AA + '\n' +
                                     str(reverse_complement(str(record.seq)[max(loci2-1-args.d,0):
                                     min(loci1+args.d,total_length)])) + '\n')
                    AA_seq[AA]=''
        except IOError:
            pass
        f1.close()
    except (IOError):
        print ('Files were missing: ' + orffile)


################################################### Programme #######################################################
Extractaa( args.i, args.f+args.n, args.f, args.r)
