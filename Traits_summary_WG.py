import os
from Bio import SeqIO
import argparse
import glob


############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-t",
                    help="trait name", type=str, default='ARG',metavar='trait name')
parser.add_argument("-db",
                    help="file name of your input database", type=str,
                    default='Butyrate.pro.aa',metavar='database.aa')
parser.add_argument("-i",
                    help="input dir of Filelist", type=str, default='.',metavar='Filelist.txt')
parser.add_argument("-m",
                    help="mapping file of traits to function", type=str,
                    default='Butyrate.pro.mapping.txt',
                    metavar='Butyrate.pro.mapping.txt')
# optional parameters
parser.add_argument("--fa",
                    help="input format of genome sequence", type=str, default='.fna.add',metavar='.fasta, .fna or .fa')
parser.add_argument("--orf",
                    help="input format of genome orfs", type=str, default='.genes.faa',metavar='.faa')
# optional output setup
parser.add_argument("--r",
                    help="output directory or folder of your results of Traits_WGD.py",
                    type=str, default='Result',metavar='Result')
parser.add_argument("--r16",
                    help="output directory or folder of your 16S sequences of Traits_WGD.py",
                    type=str, default='Result',metavar='Result')
parser.add_argument("--s",
                    help="output directory or folder of your results of traits summary",
                    type=str, default='summary',metavar='summary')
parser.add_argument("-c",
                    help="cutoff for traits",
                    type=float, default=60.0,
                    metavar='60.0')


################################################## Definition ########################################################
args = parser.parse_args()
fasta_format = args.fa
orfs_format = args.orf
if '.add' not in orfs_format:
    orfs_format = orfs_format + '.add'
try:
    os.mkdir(args.s)
except OSError:
    pass


################################################### Function ########################################################
def check_16S(inputfile):
    file16S=os.path.join(args.r16 + '/0',
                    inputfile.replace(orfs_format, fasta_format) + '.16S.txt.fasta')
    Has16S = 0
    for lines in open(file16S,'r'):
        if str(lines)!='':
            # 16S file not empty
            Has16S = 1
    if Has16S == 1:
        # merge 16S fasta
        for record in SeqIO.parse(file16S, 'fasta'):
            f16s.write('>'+str(record.id).split('_final')[0]+'\n'+str(record.seq)+'\n')
    return Has16S


def check_traits(inputfile,outputfile_aa,outputfile_aa_2000,outputfile_blast,outputfile_summary,file_subfix,i):
    blastout = args.r + '/search_output/0/' + inputfile + '.blast.txt.filter'
    aaout = blastout + '.aa'
    aaout2000 = blastout + '.extra2000.aa'
    Hastraits = 0
    for lines in open(blastout,'r'):
        if str(lines)!='':
            # blastout file not empty
            Hastraits = 1
    if Hastraits == 1:
        # merge traits fasta
        for record in SeqIO.parse(aaout, 'fasta'):
            if filename.split(file_subfix)[0] in str(record.id):
                outputfile_aa.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
            else:
                outputfile_aa.write('>'+filename.split(file_subfix)[0] + '_'+ str(record.id)+
                '\n'+str(record.seq)+'\n')
        # merge traits extend 2000 fasta
        if outputfile_aa_2000 != 'None':
            for record in SeqIO.parse(aaout2000, 'fasta'):
                if filename.split(file_subfix)[0] in str(record.id):
                    outputfile_aa_2000.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
                else:
                    outputfile_aa_2000.write('>'+filename.split(file_subfix)[0] + '_'+ str(record.id)+
                    '\n'+str(record.seq)+'\n')
        # format traits output
        totaltraits=[]
        temp = []
        while i > 0:
            totaltraits.append(0)
            i = i-1
        for lines in open(blastout, 'r'):
            if filename.split(file_subfix)[0] not in lines:
                lines = filename.split(file_subfix)[0] + '_' + lines
            traits_gene=str(lines).split('\t')[1]
            # output blast results
            if Function != dict():
                outputfile_blast.write(Function.get(traits_gene,'None')+'\t'+str(lines))
            else:
                # direct traits output
                outputfile_blast.write(str(lines))
            # calculate copy number of each gene
            if float(lines.split('\t')[2]) >= args.c:
                totaltraits[Functionlist[Function[traits_gene]]] += 1
            for copy_number in totaltraits:
                temp.append(str(copy_number))
        outputfile_summary.write(inputfile.split(file_subfix)[0]+'\t'+'\t'.join(temp)+'\n')
    return Hastraits


################################################### Programme #######################################################
# load all files
in_dir=args.i
Targetroot=dict()
for root,dirs,files in os.walk(in_dir):
    list_fasta1 = glob.glob(os.path.join(root, '*'+orfs_format))
    if list_fasta1!=[]:
        for files in list_fasta1:
            Targetroot.setdefault(files, orfs_format)


# load traits mapping file
Function = dict()
Functionlist=dict()
genenum=0
for lines in open(args.m,'r'):
    Function.setdefault(str(lines).split('\t')[0],
    str(lines).split('\t')[1].split('\r')[0].split('\n')[0])
    if (lines.split('\t')[-1].split('\r')[0].split('\n')[0]) not in Functionlist:
        Functionlist.setdefault((lines.split('\t')[-1].split('\r')[0].split('\n')[0]),genenum)
        genenum+=1

# set output files
f16s=open(os.path.join(args.s,args.t+'.all.16S.fasta'),'w')
faa=open(os.path.join(args.s,args.t+'.all.traits.aa.fasta'),'w')
ftraits=open(os.path.join(args.s,args.t+'.all.traits.aa.txt'),'w')
fdna=open(os.path.join(args.s,args.t+'.all.traits.dna.fasta'),'w')
fdna_2000=open(os.path.join(args.s,args.t+'.all.traits.dna.extra2000.fasta'),'w')
ftraits_dna=open(os.path.join(args.s,args.t+'.all.traits.dna.txt'),'w')
fsum_aa = open(os.path.join(args.s,args.t+'.all.traits.aa.summarize.'+str(args.c)+'.txt'),'w')
fsum_aa.write('SampleID\t')
for functions in Functionlist:
    fsum_aa.write('\t'+str(functions)+'_copy')
fsum_aa.write('\n')
fsum_dna = open(os.path.join(args.s,args.t+'.all.traits.dna.summarize.'+str(args.c)+'.txt'),'w')
fsum_dna.write('SampleID')
for functions in Functionlist:
    fsum_dna.write('\t'+str(functions)+'_copy')
fsum_dna.write('\n')

# process all traits
for filenames in Targetroot:
    filedir, filename = os.path.split(filenames.split('\r')[0].split('\n')[0])
    # check 16S
    if check_16S(filename) == 1:
        # check and summarize traits
        check_traits(filename,faa,'None',ftraits,fsum_aa,orfs_format,genenum)
        check_traits(filename.replace(orfs_format,fasta_format),fdna,fdna_2000,ftraits_dna,fsum_dna,fasta_format,genenum)


# end of processing all traits
f16s.close()
faa.close()
ftraits.close()
fdna.close()
fdna_2000.close()
ftraits_dna.close()
fsum_aa.close()
fsum_dna.close()
