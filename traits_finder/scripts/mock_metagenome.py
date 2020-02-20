import os
from Bio import SeqIO
import argparse
import glob
import random


############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                    help="file name of your input fasta",
                     type=str, metavar='database.aa',
                    default='/scratch/users/anniz44/genomes/GHM/new_markers/summary/31_marker.fas.all.traits.dna.fasta')

################################################## Definition ########################################################
args = parser.parse_args()


################################################### Function ########################################################
def mocking(inputfile):
    fout = open(inputfile + '.mc.fasta','w')
    fout_dep = open(inputfile + '.mc.fasta.depth', 'w')
    fout_dep_list = []
    fout_list = []
    Depth = [30,20,10,10,10,5,5,5,5,5,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,
             1,1,1]
    window = 100
    num_seq = 0
    for record in SeqIO.parse(inputfile, 'fasta'):
        ID = str(record.id)
        Seq = str(record.seq)
        length = len(Seq)
        new_depth = random.choice(Depth)
        fout_dep_list.append('%s\t%s\t\n'%(ID,new_depth))
        if length<= window:
            for j in range(0,new_depth):
                fout_list.append('>%s_%s\n%s\n'%(ID,j,Seq))
        else:
            for i in range(0,length-window):
                for j in range(0, new_depth):
                    fout_list.append('>%s_%s_%s\n%s\n' % (ID, i,j,Seq[i:(i+window)]))
        num_seq += 1
        if num_seq%10 == 0:
            fout.write(''.join(fout_list))
            fout_dep.write(''.join(fout_dep_list))
            fout_list = []
            fout_dep_list = []
    fout.write(''.join(fout_list))
    fout_dep.write(''.join(fout_dep_list))
    fout.close()
    fout_dep.close()


################################################### Programme #######################################################
mocking(args.i)
