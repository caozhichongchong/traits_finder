import os
from Bio import SeqIO
import argparse
import glob
import statistics
from datetime import datetime
import numpy as np
import random


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
                    (1: blast; 2: hmm; 3: alignment), \
                    (default \'1\' for blast search)",
                    metavar="1 or 2 or 3",
                    choices=[1, 2, 3],
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
parser.add_argument('--bcf',
                    help="Optional: complete path to bcftools if not in PATH,",
                    metavar="/usr/local/bin/bcftools",
                    action='store', default='bcftools', type=str)
parser.add_argument('--sam',
                    help="Optional: complete path to bwa if not in PATH,",
                    metavar="/usr/local/bin/samtools",
                    action='store', default='samtools', type=str)
parser.add_argument('--vcf',
                    help="Optional: complete path to bwa if not in PATH,",
                    metavar="/usr/local/bin/vcftools",
                    action='store', default='vcftools', type=str)
parser.add_argument('--vcfstats',
                    help="Optional: complete path to vcfstats if not in PATH,",
                    metavar="/usr/local/bin/vcfstats",
                    action='store', default='vcfstats', type=str)
parser.add_argument('--strainfinder',
                    help="Optional: complete path to strainfinder",
                    metavar="/scratch/users/anniz44/bin/miniconda3/bin/strainfinder",
                    action='store',
                    default='/scratch/users/anniz44/bin/miniconda3/bin/strainfinder',
                    type=str)


################################################## Definition ########################################################
args = parser.parse_args()
workingdir=os.path.abspath(os.path.dirname(__file__))
try:
    os.mkdir(os.path.join(args.r,'summary'))
except OSError:
    pass
# import strainfinder
import os, sys
if args.strainfinder!='None':
    sys.path.append(args.strainfinder)
    from strainFinder import fitStrains, genomes_given_fracs
# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3


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


def strain_finder(SNP_file):
    try:
        foutput = open(SNP_file + '.abu', 'r')
    except IOError:
        try:
            counts = np.genfromtxt(SNP_file, int, skip_header=1)
            foutput = open(SNP_file + '.abu','w')
            ## fit strains and relative abundances
            fracs, e, ll, dof = fitStrains(counts)
            best_fracs = fracs.values()[-2]
            best_fracs_new = []
            for abu in best_fracs:
                best_fracs_new.append(str('%.3f'%(abu)))
            foutput.write('%s\n' % ('\t'.join(best_fracs_new)))
            ## get ML genomes
            n = len(best_fracs)  # fitted number of strains
            k = counts.shape[1]  # number of alleles in alphabet
            perm_file = os.path.join(args.strainfinder, 'presence_matrices/strains_%d.alleles_%d.npy' % (n, k))
            allele_perm = np.load(perm_file)
            genomes = genomes_given_fracs(counts, best_fracs, e[n], alleles=['A', 'C', 'G', 'T'])
            ## output genomes
            genome_file = SNP_file + '.genome'
            try:
                np.savetxt(genome_file, genomes, '%s', '\t')
            except IOError:
                pass
        except IndexError:
            pass


def freq_call_sub(vcf_file_list):
    SNP = dict()
    SNP_count = dict()
    Position_diff = []
    Chr = ''
    Position = 0
    for vcf_file in vcf_file_list:
        for lines in open(vcf_file, 'r'):
            if not lines.startswith("#"):
                lines_set = lines.split('\t')
                Chr_new = lines_set[0]
                Position_new = int(lines_set[1])
                Chr_position = '%s_%s' % (Chr_new, Position_new)
                Allels_frq = [0, 0, 0, 0]
                SNP.setdefault(Chr_position, Allels_frq)
                SNP_count.setdefault(Chr,0)
                Allels_frq = SNP[Chr_position]
                if "INDEL" not in lines_set[7] and lines_set[6] != 'LowQual':
                    # new SNPs on Chr
                    SNP_count[Chr] += 1
                    # calculate position apart
                    if Chr_new == Chr:
                        if Position_new != Position:
                            Position_diff.append(abs(Position_new - Position) % 3)
                    Chr = Chr_new
                    Position = Position_new
                    Poly = lines_set[7].split('DP4=')[1].split(';')[0].split(',')
                    # Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
                    Sum_ref = int(Poly[0]) + int(Poly[1])
                    Sum_alt = int(Poly[2]) + int(Poly[3])
                    Allels_frq[Allels[lines_set[3]]] += Sum_ref
                    try:
                        Allels_frq[Allels[lines_set[4]]] += int(Sum_alt)
                    except KeyError:
                        allels_set = lines_set[4].split(',')
                        for allels in allels_set:
                            if allels in Allels:
                                print('multiallelic sites not separated %s %s'%(vcf_file, lines))
                                Allels_frq[Allels[allels]] += int(Sum_alt * 0.5)
                            else:
                                pass
                    SNP[Chr_position] = Allels_frq
    return [SNP,Position_diff,SNP_count]


def freq_call(vcf_file,cov_file):
    Output = []
    Output2 = []
    Output3 = []
    try:
        SNP, Position_diff, SNP_count = freq_call_sub([vcf_file])
        sample_name = cov_file.split('.fasta')[0].split('.fastq')[0].split(args.fa)[0]
        os.system('cat %s | cut -f 1,2,3 | sort -n > %s.sort' % (cov_file, cov_file))
        Depth = dict()
        for lines in open('%s.sort' % (cov_file), 'r'):
            lines = lines.replace('\r', '').replace('\n', '')
            line_set = lines.split('\t')
            Chr = line_set[0]
            Chr_position = '%s_%s' % (line_set[0], line_set[1])
            Depth.setdefault(Chr, [[], 0])  # depth coverage
            new_cov = int(line_set[2])
            # if new_cov >=2:
            Depth[Chr][0].append(new_cov)
            Depth[Chr][1] += 1
            if Chr_position not in SNP:
                Output.append('%s\t%s\t0\t0\t0\n' % (lines, line_set[2]))
            else:
                Depth_position = sum(SNP[Chr_position])
                if Depth_position > 0 and not any(Depth_position == allels for allels in SNP[Chr_position]):
                    new_line = '%s\t%s\t%s\t%s\t%s\t%s\n' % (
                    '\t'.join(line_set[0:2]), Depth_position, SNP[Chr_position][0],
                    SNP[Chr_position][1],
                    SNP[Chr_position][2],
                    SNP[Chr_position][3])
                    Output.append(new_line)
                    Output2.append(new_line)
                    Output3.append('%s\t%s\t%s\t%s\n' % (SNP[Chr_position][0],
                                                         SNP[Chr_position][1],
                                                         SNP[Chr_position][2],
                                                         SNP[Chr_position][3]))
                else:
                    Output.append('%s\t%s\t0\t0\t0\n' % (lines, line_set[2]))
        for Chr in Depth:
            Gene_length = 1000
            if Chr in Mapping:
                Gene_length = Mapping[Chr]
            temp_cov = Depth[Chr][0]
            temp_cov.sort()
            small_90 = int(len(temp_cov) * 0.9)
            Depth[Chr][0] = float(sum(temp_cov[0:small_90])) / float(small_90)
            Depth[Chr][1] = float(Depth[Chr][1]) / float(Gene_length)
            all_output.write('%s\t%s\t%s\t%s\t%s\n' % (sample_name, Chr,
                                                       '%.3f' % (Depth[Chr][0]),
                                                       '%.3f' % (Depth[Chr][1]), Gene_length))
        os.system('rm -rf %s.sort' % cov_file)
        Position_diff_new = []
        for diff in Position_diff:
            Position_diff_new.append(str(diff))
        foutput = open(cov_file + '.frq', 'w')
        foutput.write('#CHR\tPOS\tDEP\tA\tT\tG\tC\n')
        foutput.write(''.join(Output))
        foutput.close()
        foutput = open(cov_file + '.frq.snp', 'w')
        foutput.write('#CHR\tPOS\tDEP\tA\tT\tG\tC\n')
        foutput.write(''.join(Output2))
        foutput.close()
        foutput = open(cov_file + '.frq.snp.clean', 'w')
        foutput.write('#A\tT\tG\tC\n')
        foutput.write(''.join(Output3))
        foutput.close()
        foutput = open(cov_file + '.frq.snp.position.diff', 'w')
        foutput.write('\n'.join(Position_diff_new))
        foutput.close()
    except IOError:
        pass
    if args.strainfinder!= 'None':
        strain_finder(cov_file + '.frq.snp')


def SNP_dynamics(output_files, SNP_outputfile):
    SNP_dynamics_pair = dict()
    SNP_dynamics_pair_output = []
    Total = len(output_files)
    for i in range(1, Total):
        print('%s calculate SNP dynamics for %s samples by 10 iterations'
              % (datetime.now(), i))
        for repeat_time in range(0,10):
            temp_list = set()
            while len(temp_list) < i:
                temp_list.add(random.choice(output_files))
            SNP, Position_diff,SNP_count = freq_call_sub(temp_list)
            for Chr in SNP_count:
                Chr_i = '%s\t%s'%(Chr,i)
                SNP_dynamics_pair.setdefault(Chr_i,[])
                SNP_dynamics_pair[Chr_i].append(SNP_count[Chr])
    print('%s output calculate SNP dynamics'
          % (datetime.now()))
    for i in SNP_dynamics_pair:
        for Chr_i in SNP_dynamics_pair[i]:
            SNP_dynamics_pair_output.append('%s\t%s\n' %(Chr_i,numbers))
    SNP_outputfile.write(''.join(SNP_dynamics_pair_output))


################################################### Programme #######################################################
# load database length file
Mapping = dict()
try:
    for lines in open(args.db + '.length', 'r'):
            Mapping.setdefault(lines.split('\t')[0],
                               float(lines.split('\t')[1].split('\r')[0].split('\n')[0]))
except (IOError,FileNotFoundError):
    if args.s != 2:
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


# load all files
if args.s == 1:
    searchfiles=glob.glob(os.path.join(os.path.join(args.r,'search_output'),'*/*.blast.txt.filter'))
    # calculate cell number and 16S number
    os.system('python %s/Cell_number_MG.py -m %s -fa %s -r16 %s -r %s'
              % (workingdir, workingdir + "/../database/all_KO30_name.list",
                 args.fa, args.r16, args.r))
    print('python %s/Cell_number_MG.py -m %s -fa %s -r16 %s -r %s'
          % (workingdir, workingdir + "/../database/all_KO30_name.list",
             args.fa, args.r16, args.r))

    # input cell number
    Cellnum = dict()
    for lines in open(os.path.join(os.path.join(args.r, 'summary'), 'cell.copynum.all.txt'), 'r'):
        cellnum = []
        for numbers in lines.split('\t')[1:]:
            cellnum.append(float(numbers))
        averagecell = statistics.mean(cellnum)
        Cellnum.setdefault(lines.split('\t')[0], averagecell)

    # input 16S number
    Cell16S = dict()
    for lines in open(os.path.join(os.path.join(args.r, 'summary'), '16S.copynum.all.txt'), 'r'):
        cellnum = []
        averagecell = 0
        for numbers in lines.split('\t')[1:]:
            if float(numbers) > 0:
                cellnum.append(float(numbers))
        averagecell = statistics.mean(cellnum)
        if averagecell > 0:
            Cell16S.setdefault(lines.split('\t')[0], averagecell)

    # process all traits
    Trait = dict()
    Trait_fun = dict()
    for searchfile in searchfiles:
        filedir, filename = os.path.split(searchfile.split('\r')[0].split('\n')[0])
        # calculate trait
        filename_short = filename.split(args.fa)[0].split('.blast.txt.filter')[0].split('_1')[0].split('_2')[0]
        if filename_short not in Trait and '_1_2' not in filename:
            Trait.setdefault(filename_short, [])
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
elif args.s == 2:
    searchfiles = glob.glob(os.path.join(os.path.join(args.r, 'search_output'), '*/*.hmm2.txt'))
else:
    # calculate average coverage
    try:
        all_output = open((os.path.join(args.r, 'summary/all.bam.cov.sum')), 'r')
    except IOError:
        all_output = open((os.path.join(args.r, 'summary/all.bam.cov.sum')), 'w')
        all_output.write('metagenome\tgenome\tdepth\tcoverage\tgene_length\n')
        # calculate allel frequency
        output_files = glob.glob(os.path.join(args.r, 'bwa/0/*.flt.vcf'))
        i = 0
        print('%s summarizing allele frequency and SNPs for %s' % (datetime.now(), os.path.join(args.r, 'bwa/0')))
        for vcf_file in output_files:
            cov_file = vcf_file.replace('.flt.vcf', '.sorted.bam.cov')
            i += 1
            if i%100 == 0:
                print('%s summarizing allele frequency and SNPs for %s files' % (datetime.now(), i))
            freq_call(vcf_file, cov_file)
        all_output.close()
    print('%s start calculate SNP dynamics into %s' % (datetime.now(),
                                                       os.path.join(args.r, 'summary/all.snp.dynamics.txt')))
    foutput = open(os.path.join(args.r, 'summary/all.snp.dynamics.txt'), 'w')
    foutput.write('reference\tsample_number\ttotal_snps\n')
    SNP_dynamics(output_files, foutput)
    foutput.close()

