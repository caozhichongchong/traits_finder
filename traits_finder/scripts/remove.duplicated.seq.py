from Bio import SeqIO
import argparse
import os


############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                    help="input sequence name", type=str,
                    default='input.fasta',metavar='input.fasta')

################################################## Definition ########################################################
args = parser.parse_args()
input_fasta=args.i
################################################### Function #######################################################


################################################### Programme #######################################################
if 'dereplicated' not in input_fasta:
    f1 = open(input_fasta + '.dereplicated.id.fasta', 'w')
    input_id = []
    for record in SeqIO.parse(open(input_fasta, 'r'), 'fasta'):
        if str(record.id) not in input_id:
            input_id.append(str(record.id))
            f1.write('>%s\n%s\n'%(str(record.id),str(record.seq)))
    f1.close()

