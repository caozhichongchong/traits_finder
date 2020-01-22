import os
from Bio import SeqIO
import argparse
import glob
import statistics


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
parser.add_argument("-t",
                    help="trait name", type=str, default='ARG',metavar='trait name')
parser.add_argument("-fa",
                    help="input format of metagenomes sequence", type=str,
                    default='.fasta',metavar='.fasta or .fastq')
parser.add_argument('-s',
                    help="set the method to search the your database \
                    (1: blast; 2: hmm), \
                    (default \'1\' for blast search)",
                    metavar="1 or 2",
                    choices=[1, 2],
                    action='store', default=1, type=int)
parser.add_argument("-m",
                        help="mapping file of traits to function", type=str,
                        default='Butyrate.pro.mapping.txt',
                        metavar='database.mapping.txt')
# optional  setup
parser.add_argument("--meta",
                        help="metadata  of metagenomes", type=str,
                        default='None',
                        metavar='metadata.metagenomes.txt')
parser.add_argument("--r",
                    help="output directory or folder of your results",
                    type=str, default='Result_traits',metavar='Result_traits')
parser.add_argument("--r16",
                    help="output directory or folder of your 16S sequences",
                    type=str, default='Result_16S',metavar='Result_16S')


################################################## Definition ########################################################
args = parser.parse_args()
workingdir=os.path.abspath(os.path.dirname(__file__))
try:
    os.mkdir(os.path.join(args.r,'summary'))
except OSError:
    pass


################################################### Function ########################################################
def check_traits(inputfile):
    output = dict()
    fun_output = dict()
    for lines in open(inputfile,'r'):
        if args.s == 1:
            gene = lines.split('\t')[1]
            maplength = float(lines.split('\t')[3])
            genelength = Mapping[gene]
        else:
            # for hmm, we treat all genes as 1000bp
            gene = lines.split('\t')[2]
            maplength = float(lines.split('\t')[-6])
            genelength = 1000
        function_gene = Function.get(gene)
        gene_copy = maplength / genelength
        if function_gene not in fun_output:
            fun_output.setdefault(function_gene, gene_copy)
        else:
            fun_output[function_gene] += gene_copy
        if gene not in output:
            output.setdefault(gene,gene_copy)
        else:
            output[gene] += gene_copy
    return [output,fun_output]


def norm_pairend(Dataset, sample):
    if len(Dataset[sample]) > 1:
        Tempset = []
        i = 0
        while i < len(Dataset[sample][0]):
            Tempset.append((Dataset[sample][0][i] + Dataset[sample][1][i]) / 2.0)
            i += 1
        Dataset[sample] = Tempset
    else:
        Dataset[sample] = Dataset[sample][0]


def output_copy(Dataset,sample,metadata, fout1,fout2,fout3):
    # calculate copy per 16S
    if sample in Cell16S:
        fout1.write(sample)
        for numbers in Dataset[sample]:
            fout1.write('\t' + str(numbers / Cell16S[sample]))
        fout1.write(metadata + '\n')
    # calculate copy per cell
    if sample in Cellnum:
        fout2.write(sample)
        for numbers in Dataset[sample]:
            fout2.write('\t' + str(numbers / Cellnum[sample]))
        fout2.write(metadata + '\n')
    # output copy number
    fout3.write(sample)
    for numbers in Dataset[sample]:
        fout3.write('\t' + str(numbers))
    fout3.write(metadata + '\n')


def sum_output():
    fout1 = open(os.path.join(os.path.join(args.r, 'summary'), args.t + '.gene.summary.16S.txt'), 'w')
    fout2 = open(os.path.join(os.path.join(args.r, 'summary'), args.t + '.gene.summary.cell.txt'), 'w')
    fout3 = open(os.path.join(os.path.join(args.r, 'summary'), args.t + '.gene.summary.copy.txt'), 'w')
    fout4 = open(os.path.join(os.path.join(args.r, 'summary'), args.t + '.function.summary.16S.txt'), 'w')
    fout5 = open(os.path.join(os.path.join(args.r, 'summary'), args.t + '.function.summary.cell.txt'), 'w')
    fout6 = open(os.path.join(os.path.join(args.r, 'summary'), args.t + '.function.summary.copy.txt'), 'w')
    fout1.write('SampleID')
    fout2.write('SampleID')
    fout3.write('SampleID')
    fout4.write('SampleID')
    fout5.write('SampleID')
    fout6.write('SampleID')
    for traits in Mapping:
        fout1.write('\t' + traits)
        fout2.write('\t' + traits)
        fout3.write('\t' + traits)
    for functions in Function_list:
        fout1.write('\t' + functions)
        fout2.write('\t' + functions)
        fout3.write('\t' + functions)
    fout1.write('\tMetadata\n')
    fout2.write('\tMetadata\n')
    fout3.write('\tMetadata\n')
    fout4.write('\tMetadata\n')
    fout5.write('\tMetadata\n')
    fout6.write('\tMetadata\n')
    for samples in Trait:
        norm_pairend(Trait, samples)
        norm_pairend(Trait_fun, samples)
        # load metadata
        metadata = '\t'
        if samples in Meta:
            metadata += Meta[samples]
        output_copy(Trait, samples, metadata, fout1, fout2, fout3)
        output_copy(Trait_fun, samples, metadata, fout4, fout5, fout6)
    fout1.close()
    fout2.close()
    fout3.close()
    fout4.close()
    fout5.close()
    fout6.close()


################################################### Programme #######################################################
# load all files
if args.s == 1:
    searchfiles=glob.glob(os.path.join(os.path.join(args.r,'search_output'),'*/*.blast.txt.filter'))
    # gene ID: length
    Mapping = dict()
    for lines in open(args.db + '.length', 'r'):
        Mapping.setdefault(lines.split('\t')[0],
                           float(lines.split('\t')[1].split('\r')[0].split('\n')[0]))
elif args.s == 2:
    searchfiles = glob.glob(os.path.join(os.path.join(args.r, 'search_output'), '*/*.hmm2.txt'))
else:
    searchfiles = glob.glob(os.path.join(os.path.join(args.r, 'bwa'), '*/*..sorted.bam.avgcov'))
    Mapping = dict()
    try:
        for lines in open(args.db + '.length', 'r'):
            Mapping.setdefault(lines.split('\t')[0],
                               float(lines.split('\t')[1].split('\r')[0].split('\n')[0]))
    except IOError:
        Fasta_name = open(args.db, 'r')
        f = open(args.db + '.length', 'w')
        for record in SeqIO.parse(Fasta_name, 'fasta'):
            f.write(str(record.id) + '\t' + str(
                len(record.seq)) + '\n')
            Mapping.setdefault(str(record.id), len(str(record.seq)))
        f.close()


# load traits function
Function = dict()
Function_list = []
for lines in open(args.m,'r'):
    gene = str(lines).split('\t')[0]
    gene_fun = str(lines).split('\t')[1].split('\r')[0].split('\n')[0]
    Function.setdefault(gene,
    gene_fun)
    if gene_fun not in Function_list:
        Function_list.append(gene_fun)


Meta = dict()
if args.meta != 'None':
    # load metadata
    for lines in open(args.meta, 'r'):
        Meta.setdefault(lines.split('\t')[0].split(args.fa)[0].split('.blast.txt.filter')[0].split('_1')[0].split('_2')[0],
                        '\t'.join(lines.split('\r')[0].split('\n')[0].split('\t')[1:]))


# calculate cell number and 16S number
os.system('python %s/Cell_number_MG.py -m %s -fa %s -r16 %s -r %s'
          %(workingdir,workingdir + "/../database/all_KO30_name.list",
            args.fa, args.r16, args.r))
print('python %s/Cell_number_MG.py -m %s -fa %s -r16 %s -r %s'
          %(workingdir,workingdir + "/../database/all_KO30_name.list",
            args.fa, args.r16, args.r))

# input cell number
Cellnum=dict()
for lines in open(os.path.join(os.path.join(args.r,'summary'),'cell.copynum.all.txt'),'r'):
    cellnum = []
    for numbers in lines.split('\t')[1:]:
        cellnum.append(float(numbers))
    averagecell = statistics.mean(cellnum)
    Cellnum.setdefault(lines.split('\t')[0],averagecell)


# input 16S number
Cell16S=dict()
for lines in open(os.path.join(os.path.join(args.r,'summary'),'16S.copynum.all.txt'),'r'):
    cellnum = []
    averagecell = 0
    for numbers in lines.split('\t')[1:]:
        if float(numbers) > 0:
            cellnum.append(float(numbers))
    averagecell = statistics.mean(cellnum)
    if averagecell > 0:
        Cell16S.setdefault(lines.split('\t')[0],averagecell)


# process all traits
Trait = dict()
Trait_fun = dict()
for searchfile in searchfiles:
    filedir, filename = os.path.split(searchfile.split('\r')[0].split('\n')[0])
    # calculate trait
    filename_short = filename.split(args.fa)[0].split('.blast.txt.filter')[0].split('_1')[0].split('_2')[0]
    if filename_short not in Trait and '_1_2' not in filename:
        Trait.setdefault(filename_short,[])
        Trait_fun.setdefault(filename_short, [])
    Tempoutput = check_traits(searchfile)
    Tempset = []
    for traits in Mapping:
        try:
            Tempset.append(Tempoutput[0][traits])
        except KeyError:
            Tempset.append(0)
    Trait[filename_short].append(Tempset)
    Tempset = []
    for functions in Function_list:
        try:
            Tempset.append(Tempoutput[1][functions])
        except KeyError:
            Tempset.append(0)
    Trait_fun[filename_short].append(Tempset)


# summarize result
sum_output()
