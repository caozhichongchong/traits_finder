import os
from Bio import SeqIO
import argparse
import glob
import statistics

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-t",
                    help="trait name", type=str, default='ARG',metavar='trait name')
parser.add_argument("-i",
                    help="input dir of search results", type=str,
                    default='/scratch/users/anniz44/Metagenomes/ARG/search_output/0',metavar='.')
parser.add_argument("-f",
                    help="input format of searching file", type=str,
                    default='.blast.txt.filter',metavar='.blast.txt.filter')
parser.add_argument("-c",
                    help="input cell number summary", type=str,
                    default='/scratch/users/anniz44/Metagenomes/summary/cell.num.all.txt',metavar='cell.num.all.txt')
parser.add_argument("-16S",
                    help="input cell number summary", type=str,
                    default='/scratch/users/anniz44/Metagenomes/summary/16S.num.all.txt',metavar='16S.num.all.txt')
parser.add_argument("-m",
                    help="file name of trait length mapping file", type=str,
                    default='/scratch/users/anniz44/scripts/database/SARG.db.fasta.length',
                    metavar='db.fasta.length')

# optional  setup
parser.add_argument("--i16",
                    help="input dir of 16S search results", type=str,
                    default='/scratch/users/anniz44/Metagenomes/16S/0',metavar='.')
parser.add_argument("--f16",
                    help="input format of 16S searching file", type=str,
                    default='.16S.txt',metavar='.16S.txt')

parser.add_argument("--r",
                    help="output directory or folder of your results",
                    type=str, default='/scratch/users/anniz44/Metagenomes/summary',
                    metavar='Result')


################################################## Definition ########################################################
args = parser.parse_args()
try:
    os.mkdir(args.r)
except OSError:
    pass


################################################### Function ########################################################
def check_16S(inputfile):
    try:
        coverage=0
        for lines in open(os.path.join(args.i16,inputfile.replace(args.f,args.f16))):
            coverage += float(lines.split('\t')[3])/1500.0
        return coverage
    except IOError:
        return 0


def check_traits(inputfile):
    output = dict()
    for lines in open(os.path.join(args.i, inputfile)):
        gene = lines.split('\t')[1]
        maplength = float(lines.split('\t')[3])
        genelength = Mapping[gene]
        if gene not in output:
            output.setdefault(gene, maplength / genelength)
        else:
            output[gene] += maplength / genelength
    return output


################################################### Programme #######################################################
# load all files
searchfiles=glob.glob(os.path.join(args.i,'*'+args.f))

# input cell number
Cellnum=dict()
for lines in open(args.c,'r'):
    cellnum = []
    for numbers in lines.split('\t')[1:]:
        cellnum.append(float(numbers))
    averagecell = statistics.mean(cellnum)
    Cellnum.setdefault(lines.split('\t')[0].replace('.fasta','').replace('.blast.txt.filter',''),averagecell)

# input 16S number
Cell16S=dict()
for lines in open(args.c,'r'):
    cellnum = []
    averagecell = 0
    for numbers in lines.split('\t')[1:]:
        if float(numbers) > 0:
            cellnum.append(float(numbers))
    averagecell = statistics.mean(cellnum)
    if averagecell > 0:
        Cell16S.setdefault(lines.split('\t')[0].replace('.fasta','').replace('.blast.txt.filter',''),averagecell)

# gene ID: length
Mapping = dict()
for lines in open(args.m, 'r'):
    Mapping.setdefault(lines.split('\t')[0],
                       float(lines.split('\t')[1].split('\r')[0].split('\n')[0]))

# process all traits
#Cell16S=dict()
Trait=dict()
for searchfile in searchfiles:
    filedir, filename = os.path.split(searchfile.split('\r')[0].split('\n')[0])
    # check 16S
    #if filename.split('_1')[0].split('_2')[0] not in Cell16S:
    #    Cell16S.setdefault(filename.split('_1')[0].split('_2')[0],[])
    #Cell16S[filename.split('_1')[0].split('_2')[0]].append(check_16S(filename))
    # calculate trait
    filename_short = filename.split('_1')[0].split('_2')[0].replace('.fasta','').replace('.blast.txt.filter','')
    if filename_short not in Trait and '_1_2' not in filename:
        Trait.setdefault(filename_short,[])
    Tempset=[]
    Tempoutput = check_traits(filename)
    for traits in Mapping:
        try:
            Tempset.append(Tempoutput[traits])
        except KeyError:
            Tempset.append(0)
    Trait[filename_short].append(Tempset)


# output results
#f1=open(os.path.join(args.r,'16S.num.all.txt'),'w')
#for filename in Cell16S:
#    f1.write(filename)
#    for number in Cell16S[filename]:
#        f1.write('\t'+str(number))
#    f1.write('\n')
#f1.close()


# summarize result
fout1 = open(os.path.join(args.r,args.t+'.summary.16S.txt'),'w')
fout2 = open(os.path.join(args.r,args.t+'.summary.cell.txt'),'w')
fout3 = open(os.path.join(args.r,args.t+'.summary.copy.txt'),'w')
fout1.write('SampleID')
fout2.write('SampleID')
fout3.write('SampleID')
for traits in Mapping:
    fout1.write('\t'+traits)
    fout2.write('\t' + traits)
    fout3.write('\t' + traits)
fout1.write('\n')
fout2.write('\n')
fout3.write('\n')
for samples in Trait:
    if len(Trait[samples]) > 1:
        Tempset = []
        i = 0
        while i < len(Trait[samples][0]):
            Tempset.append((Trait[samples][0][i]+Trait[samples][1][i])/2.0)
            i+=1
        Trait[samples]=Tempset
    else:
        Trait[samples]=Trait[samples][0]
    if samples in Cell16S:
        fout1.write(samples)
        for numbers in Trait[samples]:
            fout1.write('\t'+str(numbers/Cell16S[samples]))
        fout1.write('\n')
    if samples in Cellnum:
        fout2.write(samples)
        for numbers in Trait[samples]:
            fout2.write('\t'+str(numbers/Cellnum[samples]))
        fout2.write('\n')
    fout3.write(samples)
    for numbers in Trait[samples]:
        fout3.write('\t' + str(numbers))
    fout3.write('\n')

fout1.close()
fout2.close()
fout3.close()
