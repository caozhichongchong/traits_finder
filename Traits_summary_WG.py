import os
from Bio import SeqIO
import argparse
import glob


############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-db",
                    help="file name of your input database", type=str, default='Butyrate.pro.aa',metavar='database.aa')
parser.add_argument("-i",
                    help="input dir of Filelist", type=str, default='.',metavar='Filelist.txt')
# optional parameters
parser.add_argument("--fa",
                    help="input format of genome sequence", type=str, default='.fna.add',metavar='.fasta, .fna or .fa')
parser.add_argument("--orf",
                    help="input format of genome orfs", type=str, default='.genes.faa',metavar='.faa')
parser.add_argument("--m",
                    help="mapping file of traits to function", type=str, default='None',metavar='Butyrate.pro.mapping.txt')
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


################################################## Definition ########################################################
args = parser.parse_args()
fasta_format = args.fa
orfs_format = args.orf
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


def check_traits(inputfile):
    blastout = args.r + '/search_output/0/' + inputfile + '.add.blast.txt.filter'
    aaout = blastout + '.aa'
    Hastraits = 0
    for lines in open(blastout,'r'):
        if str(lines)!='':
            # blastout file not empty
            Hastraits = 1
    if Hastraits == 1:
        traits_inputfile = []
        # merge traits fasta
        for record in SeqIO.parse(aaout, 'fasta'):
            faa.write('>'+str(record.id).split(orfs_format)[0]+'\t'
                          +str(record.id).split(orfs_format+'_')[1]+'\n'+str(record.seq)+'\n')
        # format traits output
        for lines in open(blastout, 'r'):
            traits_gene=str(lines).split('\t')[1]
            # logic of functional genes
            if Traits != dict():
                if traits_gene in Traits:
                    traits_inputfile.append(Traits.get(traits_gene))
                else:
                    print([traits_gene,blastout])
            else:
                # direct traits output
                ftraits.write(str(lines).split('\t')[0] + '\t1\t' + str(lines).split('\t')[1]+'\n')
        if traits_inputfile != []:
            # output functional genes with logic
            #if any(str(logic.replace('\'','')) in traits_inputfile for logic in mapping_logic):
            if 'But' in traits_inputfile or ('Ptb' in traits_inputfile and 'Buk' in traits_inputfile):
                ftraits.write(str(inputfile).split(orfs_format)[0] + '\t1\n')
            else:
                ftraits.write(str(inputfile).split(orfs_format)[0] + '\t0\n')
    else:
        # files with no traits
        ftraits.write(str(inputfile).split(orfs_format)[0] + '\t0\n' )
    return Hastraits


################################################### Programme #######################################################
# load all files
flist=open(args.i,'r')
Targetroot=dict()
f16s=open(os.path.join(args.s,'all.16S.fasta'),'w')
faa=open(os.path.join(args.s,'all.traits.aa.fasta'),'w')
ftraits=open(os.path.join(args.s,'all.traits.txt'),'w')


# load traits mapping file: trait_gene_ID, gene_function
Traits = dict()
mapping_logic = []
if args.m != 'None':
    for lines in open(args.m,'r'):
        if '#' in str(lines):
            mapping_logic.append(str(lines).split('#')[1].split('\r')[0].split('\n')[0].replace(' ','').split('and'))
        else:
            Traits.setdefault(str(lines).split('\t')[0],str(lines).split('\t')[1].split('\r')[0].split('\n')[0])


# process all traits
for filenames in flist:
    filedir, filename = os.path.split(filenames.split('\r')[0].split('\n')[0])
    # check 16S
    if check_16S(filename) == 1:
        # check and summarize traits
        check_traits(filename)


# end of processing all traits
f16s.close()
faa.close()
ftraits.close()