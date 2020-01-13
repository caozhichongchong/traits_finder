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
                    help="input dir of genomes", type=str,
                    default='.',metavar='.')
parser.add_argument("-m",
                    help="mapping file of traits to function", type=str,
                    default='Butyrate.pro.mapping.txt',
                    metavar='Butyrate.pro.mapping.txt')
parser.add_argument("-dbf",
                        help="sequence format of your input database\
                        (1: nucleotide; 2: protein), \
                        (default \'1\' for nucleotide)",
                        metavar="1 or 2",
                        choices=[1, 2],
                        action='store', default=1, type=int)
parser.add_argument("--mge",
                        help="whether your input sequences are genomes or mobile genetic elements (MGEs)\
                        (1: genomes; 2: mge), \
                        (default \'1\' for genomes)",
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
                    help="input directory or folder of your previous results by Traits_WGD.py",
                    type=str, default='Result',metavar='Result')
parser.add_argument("--r16",
                    help="input directory or folder of your previous 16S sequences extracted by Traits_WGD.py",
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
try:
    os.mkdir(os.path.join(args.s,'sub_sequences/'))
except OSError:
    pass

################################################### Function ########################################################
def check_file(listoffiles):
    for files in listoffiles:
        try:
            f1 = open(files,'r')
            return files
        except FileNotFoundError:
            pass


def check_16S(inputfile):
    try:
        file16S= check_file(glob.glob(os.path.join(args.r16 + '/*',
                        inputfile.replace(orfs_format, fasta_format) + '.16S.txt.fasta')))
        Has16S = 0
        if file16S != None:
            for lines in open(file16S,'r'):
                if str(lines)!='':
                    # 16S file not empty
                    Has16S = 1
                break
        if Has16S == 1:
            # merge 16S fasta
            for record in SeqIO.parse(file16S, 'fasta'):
                f16s.write('>'+str(record.id).split('_final')[0]+'\n'+str(record.seq)+'\n')
    except FileNotFoundError:
        pass


def check_traits(inputfile,outputfile_aa,outputfile_aa_500,outputfile_blast,outputfile_summary,
outputfile_summary_fun,file_subfix,i):
    blastout = check_file(glob.glob(args.r + '/search_output/*/'+ inputfile + '.blast.txt.filter'))
    Hastraits = 0
    if blastout != None:
        for lines in open(blastout,'r'):
            if str(lines)!='':
                # blastout file not empty
                Hastraits = 1
            break
    if Hastraits == 1:
        Functionset = dict()
        # format traits output for genes
        totaltraits=[]
        temp = []
        while i > 0:
            totaltraits.append(0)
            i = i-1
        # format traits output for functions
        for functions in allfunction:
            allfunction[functions]=0
        for lines in open(blastout, 'r'):
            if inputfile.split(file_subfix)[0] not in lines:
                lines = inputfile.split(file_subfix)[0] + '_' + lines
            lines_set = str(lines).split('\t')
            traits_gene=lines_set[1]
            # output blast results
            if Function != dict():
                if args.mge == 2:
                    # MGEs
                    outputfile_blast.write(Function.get(traits_gene,'None')+'\tmge_'+str(lines))
                else:
                    # genomes
                    outputfile_blast.write(Function.get(traits_gene,'None')+'\t'+str(lines))
                if outputfile_aa_500 == 'None':
                    # amino acid search results
                    Functionset.setdefault(lines_set[0],
                                           Function.get(traits_gene,'None'))
                else:
                    # dna search results
                    loci1=min(int(lines_set[6]),int(lines_set[7].split('\r')[0].split('\n')[0]))
                    loci2=max(int(lines_set[6]),int(lines_set[7].split('\r')[0].split('\n')[0]))
                    Functionset.setdefault('%s_%s_%s' %(lines_set[0],
                    str(loci1-1),str(loci2)), Function.get(traits_gene,'None'))
            else:
                # direct traits output
                if args.mge == 2:
                    # MGEs
                    outputfile_blast.write('mge_'+str(lines))
                else:
                    # genomes
                    outputfile_blast.write(str(lines))
            # calculate copy number of each gene and each function
            if float(lines.split('\t')[2]) >= args.c:
                try:
                    totaltraits[Functionlist[traits_gene]] += 1
                    allfunction[Function.get(traits_gene,'None')] += 1
                except KeyError: #reference genes not in database
                    pass
        for copy_number in totaltraits:
            temp.append(str(copy_number))
        outputfile_summary.write(inputfile.split(file_subfix)[0]+'\t'+'\t'.join(temp)+'\n')
        outputfile_summary_fun.write(inputfile.split(file_subfix)[0])
        for functions in allfunction:
            outputfile_summary_fun.write('\t'+str(allfunction[functions]))
        outputfile_summary_fun.write('\n')
        # extract sequences
        aaout = blastout + '.aa'
        aaout500 = blastout + '.extra*.aa'
        # merge traits fasta into functions
        for record in SeqIO.parse(aaout, 'fasta'):
                if inputfile.split(file_subfix)[0] in str(record.id):
                    if str(record.id) in Functionset:
                        # output according to its function
                        outputfile_aa_file = open(outputfile_aa.replace('fasta',Functionset[str(record.id)]+'.fasta'),'a')
                        if args.mge == 2:
                            # MGEs
                            outputfile_aa_file.write('>mge_%s\n%s\n'%(
                                str(record.id),
                                str(record.seq)))
                        else:
                            # genomes
                            outputfile_aa_file.write('>%s\n%s\n'%(
                                str(record.id),
                                str(record.seq)))
                        outputfile_aa_file.close()
                    else:
                        print('%s not found in blast output'%(str(record.id)))
                else:
                    if inputfile.split(file_subfix)[0] + '_'+ str(record.id) in Functionset:
                        # output according to its function
                        function_name = Functionset[inputfile.split(file_subfix)[0] +
                                                                                    '_' + str(record.id)]
                        outputfile_aa_file = open(outputfile_aa.replace('fasta',function_name
                                                                        +'.fasta'),'a')
                        if args.mge == 2:
                            # MGEs
                            outputfile_aa_file.write('>mge_%s_%s\n%s\n'%(
                                inputfile.split(file_subfix)[0],
                                str(record.id),
                                str(record.seq)))
                        else:
                            # genomes
                            outputfile_aa_file.write('>%s_%s\n%s\n'%(
                                inputfile.split(file_subfix)[0],
                                str(record.id),
                                str(record.seq)))
                        outputfile_aa_file.close()
                    else:
                        print('%s not found in blast output'%(inputfile.split(file_subfix)[0] + '_'+ str(record.id)))
        # merge traits extend 500 fasta
        if outputfile_aa_500 != 'None':
            for record in SeqIO.parse(glob.glob(aaout500)[0], 'fasta'):
                if args.mge == 2:
                    # MGEs
                    if inputfile.split(file_subfix)[0] in str(record.id):
                        outputfile_aa_500.write('>mge_'+str(record.id)+'\n'+str(record.seq)+'\n')
                    else:
                        outputfile_aa_500.write('>mge_'+inputfile.split(file_subfix)[0] + '_'+ str(record.id)+
                        '\n'+str(record.seq)+'\n')
                else:
                    # genomes
                    if inputfile.split(file_subfix)[0] in str(record.id):
                        outputfile_aa_500.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
                    else:
                        outputfile_aa_500.write('>'+inputfile.split(file_subfix)[0] + '_'+ str(record.id)+
                        '\n'+str(record.seq)+'\n')
    else:
        outputfile_summary.write(inputfile.split(file_subfix)[0] + '\tNo_hit\n')
        outputfile_summary_fun.write(inputfile.split(file_subfix)[0] + '\tNo_hit\n')


################################################### Programme #######################################################
# load all files
in_dir=args.i
Targetroot=dict()
for root,dirs,files in os.walk(in_dir):
    list_fasta1 = glob.glob(os.path.join(root, '*'+orfs_format))
    if list_fasta1!=[]:
        for files in list_fasta1:
            Targetroot.setdefault(files, orfs_format)


# set output files
# 16S sequences
f16s=open(os.path.join(args.s,args.t+'.all.16S.fasta'),'w')
# blast summary
ftraits=open(os.path.join(args.s,args.t+'.all.traits.aa.txt'),'w')
ftraits_dna=open(os.path.join(args.s,args.t+'.all.traits.dna.txt'),'w')
# traits summary
fsum_aa = open(os.path.join(args.s,args.t+'.all.traits.aa.summarize.gene.'+str(args.c)+'.txt'),'w')
fsum_aa_fun = open(os.path.join(args.s,args.t+'.all.traits.aa.summarize.function.'+str(args.c)+'.txt'),'w')
fsum_dna = open(os.path.join(args.s,args.t+'.all.traits.dna.summarize.gene.'+str(args.c)+'.txt'),'w')
fsum_dna_fun = open(os.path.join(args.s,args.t+'.all.traits.dna.summarize.function.'+str(args.c)+'.txt'),'w')
# sequences
faa=os.path.join(args.s,'sub_sequences/'+args.t+'.all.traits.aa.fasta')
fdna=os.path.join(args.s,'sub_sequences/'+args.t+'.all.traits.dna.fasta')
fdna_500=open(os.path.join(args.s,args.t+'.all.traits.dna.extra500.fasta'),'w')
# reset output sequence files for all sequences
outputfile_aa_file = open(faa,'w')
outputfile_aa_file.close()
outputfile_aa_file = open(fdna,'w')
outputfile_aa_file.close()


# load traits mapping file
Function = dict()
Functionlist=dict()
genenum=0
allfunction=dict()
for lines in open(args.m,'r'):
    lines_set = str(lines).split('\t')
    gene = lines_set[0]
    gene_fun = lines_set[1].split('\r')[0].split('\n')[0]
    Function.setdefault(gene,
    gene_fun)
    if gene_fun not in allfunction:
        allfunction.setdefault(gene_fun,0)
    # reset output sequence files
    outputfile_aa_file = open(faa.replace('fasta',gene_fun+'.fasta'),'w')
    outputfile_aa_file.close()
    outputfile_aa_file = open(fdna.replace('fasta',gene_fun+'.fasta'),'w')
    outputfile_aa_file.close()
    if gene not in Functionlist:
        Functionlist.setdefault(gene,genenum)
        genenum+=1


# merge reference sequences and output sequences (amino acid only)
if args.dbf == 2:
    for record in SeqIO.parse(args.db, 'fasta'):
            if str(record.id) in Function:
                outputfile_aa_file = open(faa.replace('fasta',Function[str(record.id)]+'.fasta'),'a')
                outputfile_aa_file.write('>reference_'+str(record.id)+'\n'+str(record.seq)+'\n')
                outputfile_aa_file.close()
else:
    for record in SeqIO.parse(args.db, 'fasta'):
            if str(record.id) in Function:
                outputfile_aa_file = open(fdna.replace('fasta',Function[str(record.id)]+'.fasta'),'a')
                outputfile_aa_file.write('>reference_'+str(record.id)+'\n'+str(record.seq)+'\n')
                outputfile_aa_file.close()


# output
fsum_aa.write('SampleID')
for functions in Functionlist:
    fsum_aa.write('\t'+str(functions))
fsum_aa.write('\n')
fsum_aa_fun.write('SampleID')
for functions in allfunction:
    fsum_aa_fun.write('\t'+str(functions))
fsum_aa_fun.write('\n')
fsum_dna.write('SampleID')
for functions in Functionlist:
    fsum_dna.write('\t'+str(functions))
fsum_dna.write('\n')
fsum_dna_fun.write('SampleID')
for functions in allfunction:
    fsum_dna_fun.write('\t'+str(functions))
fsum_dna_fun.write('\n')

# process all traits
for filenames in Targetroot:
    filedir, filename = os.path.split(filenames)
    # check 16S
    check_16S(filename)
    # check and summarize traits
    check_traits(filename,faa,'None',
    ftraits,fsum_aa,fsum_aa_fun,orfs_format,genenum)
    check_traits(filename.replace(orfs_format,fasta_format),fdna,fdna_500,
    ftraits_dna,fsum_dna,fsum_dna_fun,fasta_format,genenum)


# merge all sequences into one file
os.system('cat %s > %s' %(faa.replace('fasta','*.fasta'),
                          os.path.join(args.s, os.path.split(faa)[-1])))
os.system('cat %s > %s' %(fdna.replace('fasta','*.fasta'),
                          os.path.join(args.s, os.path.split(fdna)[-1])))


# end of processing all traits
f16s.close()
ftraits.close()
fdna_500.close()
ftraits_dna.close()
fsum_aa.close()
fsum_dna.close()
fsum_aa_fun.close()
fsum_dna_fun.close()
