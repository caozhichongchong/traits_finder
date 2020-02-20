import os
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import glob
import statistics
from datetime import datetime
import numpy as np
import random
import operator
#from Bio.codonalign.codonseq import _get_codon_list, CodonSeq, cal_dn_ds


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
parser.add_argument('--th',
                    help="Optional: set the thread number assigned for running XXX (default 40)",
                    metavar="1 or more", action='store', default=40, type=int)
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
try:
    os.mkdir(os.path.join(args.r,'summary/vcf'))
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
Allels_order = ['A','T','G','C']
# Set up N or S
N_S_set = dict()
N_S_set['N']=0
N_S_set['S']=1
purines=['A','G']
pyrimidines=['C','T']
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


################################################### new class #########################################################
__metaclass__ = type


class SNP_gene:
    # create a class to store SNP_gene
    'a class to store SNP_gene'

    def init(self, gene):
        self.gene = gene
        self.position = []
        # [[N,S],freq]
        self.SNP_pair = {'A-T': [0,0],
                         'A-C': [0,0],
                         'G-C': [0,0],
                         'G-T': [0,0],
                         'A-G': [0,0],
                         'G-A': [0,0]}
        self.SNP_pair_freq = {'A-T': 0,
                         'A-C': 0,
                         'G-C': 0,
                         'G-T': 0,
                         'A-G': 0,
                         'G-A': 0}
        self.depth = []
        self.cov = 0
        self.mutposition = []
        self.NSratio = [0,0]

    def addmutposition(self, position):
        self.mutposition.append(position)

    def addposition(self, position,depth):
        self.position.append(position)
        self.cov += 1
        self.depth.append(depth)

    def addSNP_pair(self, pair,position,count,depth):
        self.SNP_pair[pair][position] += count
        self.SNP_pair_freq[pair] += count/depth
        self.NSratio[position] += count

    def sum_SNP_pair(self):
        self.SNP_pair_sum = {'A-T': [0, 0],
                         'A-C': [0, 0],
                         'G-C': [0, 0],
                         'G-T': [0, 0],
                         'A-G': [0, 0],
                         'G-A': [0, 0]}
        for pair in self.SNP_pair:
            self.SNP_pair_sum[pair][0] += self.SNP_pair[pair][0] * self.SNP_pair_freq[pair]
            self.SNP_pair_sum[pair][1] += self.SNP_pair[pair][1] * self.SNP_pair_freq[pair]

    def dN_dS(self):
        self.expectNSratio = 'No_expect'
        expectNSratio = [0, 0]
        refSNP_pair_sum = Ref_NSratio.get(self.gene, 'None')
        if refSNP_pair_sum != 'None':
            for pair in refSNP_pair_sum:
                expectNSratio[0] += refSNP_pair_sum[pair][0] * self.SNP_pair_freq[pair]
                expectNSratio[1] += refSNP_pair_sum[pair][1] * self.SNP_pair_freq[pair]
            try:
                self.expectNSratio = expectNSratio[0] / expectNSratio[1]
            except ZeroDivisionError:
                if self.NSratio[0] == 0:
                    self.expectNSratio = 'expect_None'
                else:
                    self.expectNSratio = 'expect_N_only'
        try:
            self.NSratiosum = self.NSratio[0] / self.NSratio[1]
        except (ZeroDivisionError, TypeError):
            if self.NSratio[0] == 0:
                self.NSratiosum = 'observe_None'
            else:
                self.NSratiosum = 'observe_N_only'
        self.dNdS = self.expectNSratio
        if type(self.expectNSratio) == float and type(self.NSratiosum) == float:
            self.dNdS = self.NSratiosum / self.expectNSratio
            self.dNdS = '%.3f' % (self.dNdS)
            self.NSratiosum = '%.3f' % (self.NSratiosum)
            self.expectNSratio = '%.3f' % (self.expectNSratio)
        elif type(self.expectNSratio) == float :
            self.expectNSratio = '%.3f' % (self.expectNSratio)
        elif type(self.NSratiosum) == float:
            self.NSratiosum = '%.3f' % (self.NSratiosum)


################################################### Function ########################################################
# Set DNA translate to protein
def translate(seq):
    seq = Seq(seq)
    try:
        return seq.translate()
    except ValueError:
        return seq.translate(seq.complement())


def dnORds(amino1, amino2):
    if amino1 == amino2:
        return 'S'
    else:
        return 'N'


def causeSNP(seq,position,ALT):
    seq = list(seq)
    seq[position - 1]=ALT
    return ''.join(seq)


def ALT_freq(Allels_count):
    Major_ALT = []
    Minor_ALT = []
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 4):
        ALT_frq = int(Allels_count[alleles])
        if ALT_frq > 0:
            ALT_set.setdefault(ALT_frq, set())
            ALT_set[ALT_frq].add(alleles)
            ALT_frq_set.add(ALT_frq)
    ALT_frq_set = sorted(ALT_frq_set,reverse=True)
    for ALT_frq in ALT_frq_set:
        for alleles in ALT_set[ALT_frq]:
            if Major_ALT == []:
                Major_ALT = [Allels_order[alleles],ALT_frq]
            else:
                Minor_ALT.append([Allels_order[alleles],ALT_frq])
    return [Major_ALT,Minor_ALT]


def transitions(REF,ALT):
    if REF in pyrimidines:
        REF = complement[REF]
        ALT = complement[ALT]
    return '%s-%s'%(REF,ALT)
    #if (REF in purines and ALT in purines):
    #    return [0,[REF,ALT]]
    #elif (REF in pyrimidines and ALT in pyrimidines):
    #    return [0,[REF,ALT]]
    #return [1,[REF,ALT]]


def expectNSsub(record_name,record_seq,position=0):
    Total = len(record_seq)
    temp_SNP_gene = SNP_gene()
    temp_SNP_gene.init(record_name)
    for i in range(0, (int(Total / 3) - position)):
        codon = record_seq[(i * 3 + position):((i + 1) * 3 + position)]
        codon_NSratio = codontable_NSratio[codon]
        for pair in codon_NSratio.SNP_pair:
            temp_SNP_gene.addSNP_pair(pair, 0, codon_NSratio.SNP_pair[pair][0], 1)
            temp_SNP_gene.addSNP_pair(pair, 1, codon_NSratio.SNP_pair[pair][1], 1)
    temp_SNP_gene.sum_SNP_pair()
    return temp_SNP_gene.SNP_pair_sum

def expectNS(record_name,record_seq):
    try:
        return expectNSsub(record_name, record_seq)
    except KeyError:
        try:
            return expectNSsub(record_name, record_seq,1)
        except KeyError:
            try:
                return expectNSsub(record_name, record_seq,2)
            except KeyError:
                return 'None'


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
    except (IOError, FileNotFoundError):
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
            except (IOError, FileNotFoundError):
                pass
        except IndexError:
            pass


def freq_call_sub(vcf_file_list, alloutput = True):
        SNP = dict()
        SNP_count = dict()
        Position_diff = []
        Chr = ''
        Position = 0
        for vcf_file in vcf_file_list:
            for lines in open(vcf_file, 'r'):
                if not lines.startswith("#"):
                    lines_set = lines.split('\t')
                    Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
                    Chr_new = lines_set[0]
                    SNP_count.setdefault(Chr_new, 0)
                    if alloutput:
                        Position_new = int(lines_set[1])
                        Chr_position = '%s--%s' % (Chr_new, Position_new)
                        Allels_frq = [0, 0, 0, 0]
                        SNP.setdefault(Chr_position, Allels_frq)
                        Allels_frq = SNP[Chr_position]
                    if "INDEL" not in lines_set[7] and (lines_set[6] != 'LowQual' or Depth >= 10):
                        #Depth >= 10 for genome or Depth >= 100 for metagenomes
                        # new SNPs on Chr
                        SNP_count[Chr_new] += 1
                        if alloutput:
                            # calculate position apart
                            if Chr_new == Chr:
                                if Position_new != Position:
                                    Position_diff.append(abs(Position_new - Position) % 3)
                            Chr = Chr_new
                            Position = Position_new
                            allels_set = [lines_set[3]]
                            if '.' not in lines_set[4]:
                                allels_set += lines_set[4].split(',')
                            Subdepth = lines_set[9].split(':')[-1].replace('\n', '').split(',')
                            for num_allels in range(0, len(allels_set)):
                                allels = allels_set[num_allels]
                                Subdepth_alleles = Subdepth[num_allels]
                                if allels in Allels:
                                    Allels_frq[Allels[allels]] += int(Subdepth_alleles)
                                else:
                                    pass
                            SNP[Chr_position] = Allels_frq
        return [SNP,Position_diff,SNP_count]


def dN_dS_ratio(seq1,seq2,method = 'ML'):
    if method=='ML':
        return cal_dn_ds(CodonSeq(seq1), CodonSeq(seq2), method='ML')
    elif method == 'NG86':
        return cal_dn_ds(CodonSeq(seq1), CodonSeq(seq2), method='NG86')
    elif method == 'LWL85':
        return cal_dn_ds(CodonSeq(seq1), CodonSeq(seq2), method='LWL85')
    elif method == 'YN00':
        return cal_dn_ds(CodonSeq(seq1), CodonSeq(seq2), method='YN00')


def freq_call(vcf_file,cov_file):
    Output =  []
    Output2 = []
    Output3 = []
    all_SNP_gene_temp = dict()
    try:
        SNP, Position_diff, SNP_count = freq_call_sub([vcf_file])
        SNP_all, Position_diff_all, SNP_count_all = freq_call_sub([cov_file])
        sample_name = os.path.split(cov_file)[1].split('.fasta')[0].split('.fastq')[0].split(args.fa)[0]
        #Depth = dict()
        for Chr_position in SNP_all:
            Depth_position = sum(SNP_all[Chr_position])
            Chr = Chr_position.split('--')[0]
            if Chr not in all_SNP_gene_temp:
                SNP_gene_temp = SNP_gene()
                SNP_gene_temp.init(Chr)
                all_SNP_gene_temp.setdefault(Chr,SNP_gene_temp)
            SNP_gene_temp = all_SNP_gene_temp[Chr]
            #N_S.setdefault(Chr, [0, 0])
            position = Chr_position.split('--')[1]
            SNP_gene_temp.addposition(position, Depth_position)
            lines = '%s\t%s\t%s\t%s\t%s\t%s\t%s' % (Chr,position,
                    Depth_position, SNP_all[Chr_position][0],
                    SNP_all[Chr_position][1],
                    SNP_all[Chr_position][2],
                    SNP_all[Chr_position][3])
            Output.append(lines + '\n')
            if Chr_position in SNP:
                position = int(position)
                Depth_position_snp = sum(SNP[Chr_position])
                if Depth_position_snp > 0 and not any(Depth_position_snp == allels for allels in SNP[Chr_position]):
                    # a SNP
                    # calculate N or S
                    codon_start = position - 1 - int(position%3)
                    Ref_seq_chr = Ref_seq[Chr]
                    SNP_seq_chr = Ref_seq_chr
                    Major_ALT, Minor_ALT = ALT_freq(SNP[Chr_position])
                    REF = Major_ALT[0]
                    Ref_frq = Major_ALT[1]
                    Ref_seq_chr = causeSNP(Ref_seq_chr, position, REF)
                    Ref_seq_codon = translate(Ref_seq_chr[codon_start:(codon_start + 3)])[0]
                    for minor in Minor_ALT:
                        ALT =  minor[0]
                        ALT_frq =  minor[1]
                        SNP_seq_chr = causeSNP(SNP_seq_chr, position,ALT)
                        SNP_seq_codon = translate(SNP_seq_chr[codon_start:(codon_start + 3)])[0]
                        temp_NorS = dnORds(Ref_seq_codon, SNP_seq_codon)
                        SNP_pair = transitions(REF, ALT)
                        SNP_gene_temp.addSNP_pair(SNP_pair,N_S_set[temp_NorS],
                                                  ALT_frq,
                                                  Depth_position_snp)
                        SNP_gene_temp.addmutposition(position)
                        #result = dN_dS_ratio(Ref_seq_chr, SNP_seq_chr)
                        #try:
                        #    N_S[Chr][-1].append('%.3f' % (result[0] / result[1]))
                        #except ZeroDivisionError:
                        #    N_S[Chr][-1].append('N_only')
                        #lines += '\t%s:%s:%s\t%s:%s:%s\t%s\t%s:%s' % (REF, Ref_frq, Ref_seq_codon,
                        #                                      ALT,ALT_frq, SNP_seq_codon, temp_NorS,result[0],result[1])
                        lines += '\t%s:%s:%s\t%s:%s:%s\t%s' % (REF, Ref_frq, Ref_seq_codon,
                                                                     ALT, ALT_frq, SNP_seq_codon, temp_NorS)
                    Output2.append(lines + '\n')
                    Output3.append('%s\t%s\t%s\t%s\n' % (SNP[Chr_position][0],
                                                         SNP[Chr_position][1],
                                                         SNP[Chr_position][2],
                                                         SNP[Chr_position][3]))
        all_output_list = []
        for Chr in all_SNP_gene_temp:
            SNP_gene_temp = all_SNP_gene_temp[Chr]
            Gene_length = 1000
            if Chr in Mapping:
                Gene_length = Mapping[Chr]
            temp_depth = SNP_gene_temp.depth
            temp_depth.sort()
            temp_std = statistics.stdev(temp_depth)
            temp_mean = statistics.mean(temp_depth)
            temp_range = [temp_mean-1.5*temp_std,temp_mean+1.5*temp_std]
            temp_depth_filter = []
            for position in range(0,len(temp_depth)):
                a = temp_depth[position]
                if (temp_range[0] <= a and a <= temp_range[1]):
                    temp_depth_filter.append(a)
                else:
                    pass
                    #Output.pop(Depth[Chr][2][position], None)
                    # Output2.pop(Depth[Chr][2][position], None)
                    # Output3.pop(Depth[Chr][2][position], None)
            SNP_gene_temp.depth = statistics.mean(temp_depth_filter)
            SNP_gene_temp.cov = float(SNP_gene_temp.cov) / float(Gene_length)
            new_line = '%s\t%s\t%s\t%s\t%s' % (sample_name, Chr,
                                                       '%.3f' % (SNP_gene_temp.depth),
                                                       '%.3f' % (SNP_gene_temp.cov), Gene_length)
            N_temp = SNP_gene_temp.NSratio[0]
            S_temp = SNP_gene_temp.NSratio[1]
            if N_temp + S_temp > 0:
                SNP_gene_temp.dN_dS()
                new_line += ('\t%s:%s\t%s' % (N_temp, S_temp, SNP_gene_temp.NSratiosum))
                new_line += '\t%s\t%s' % (SNP_gene_temp.expectNSratio,SNP_gene_temp.dNdS)
                for pair in SNP_gene_temp.SNP_pair:
                    pair_freq = SNP_gene_temp.SNP_pair_freq[pair]
                    pair_N = SNP_gene_temp.SNP_pair[pair][0]
                    pair_S = SNP_gene_temp.SNP_pair[pair][1]
                    try:
                        pair_NSratio = pair_N/pair_S
                    except ZeroDivisionError:
                        if pair_N == 0:
                            pair_NSratio = 'None'
                        else:
                            pair_NSratio = 'N_only'
                    new_line += ('\t%s\t%s\t%s:%s\t%s' % (pair,pair_freq,pair_N,pair_S,pair_NSratio))
            new_line += '\n'
            all_output_list.append(new_line)
        all_output.write(''.join(all_output_list))
        #os.system('rm -rf %s.sort' % cov_file)
        Position_diff_new = []
        for diff in Position_diff:
            Position_diff_new.append(str(diff))
        Output = sorted(Output, key=operator.itemgetter(1))
        Output = sorted(Output, key=operator.itemgetter(0))
        Output2 = sorted(Output2, key=operator.itemgetter(1))
        Output2 = sorted(Output2, key=operator.itemgetter(0))
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
    except (IOError, FileNotFoundError):
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
            SNP, Position_diff,SNP_count = freq_call_sub(temp_list,False)
            for Chr in SNP_count:
                Chr_i = '%s\t%s'%(Chr,i)
                SNP_dynamics_pair.setdefault(Chr_i,[])
                SNP_dynamics_pair[Chr_i].append(SNP_count[Chr])
    i = Total
    SNP, Position_diff, SNP_count = freq_call_sub(output_files,False)
    for Chr in SNP_count:
        Chr_i = '%s\t%s' % (Chr, i)
        SNP_dynamics_pair.setdefault(Chr_i, [])
        SNP_dynamics_pair[Chr_i].append(SNP_count[Chr])
    print('%s output calculate SNP dynamics'
          % (datetime.now()))
    for Chr_i in SNP_dynamics_pair:
        for numbers in SNP_dynamics_pair[Chr_i]:
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
    # load codon freq
    # Set up codon table
    codontable = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    codontable_NSratio = dict()
    for codon in codontable:
        SNP_gene_temp = SNP_gene()
        SNP_gene_temp.init(codon)
        codontable_NSratio.setdefault(codon,SNP_gene_temp)
        for position in range(0, 3):
            REF = codon[position]
            for ALT in Allels_order:
                if ALT != REF:
                    new_codon = causeSNP(codon, position, ALT)
                    temp_NorS = dnORds(translate(codon)[0], translate(new_codon)[0])
                    SNP_pair = transitions(REF, ALT)
                    SNP_gene_temp.addSNP_pair(SNP_pair, N_S_set[temp_NorS], 1,1)
    # load database seq
    Ref_seq = dict()
    Ref_NSratio = dict()
    for record in SeqIO.parse(args.db, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        Ref_seq.setdefault(record_id,record_seq)
        Ref_NSratio.setdefault(record_id,
                               expectNS(record_id,record_seq))
    print(Ref_NSratio)
    # all bam file call vcf
    tempbamoutput = glob.glob(os.path.join(args.r, 'summary/vcf/all.flt.vcf*'))
    if tempbamoutput == []:
        tempbamoutput = os.path.join(args.r, 'summary/vcf/all.flt.vcf')
        output_files = glob.glob(os.path.join(args.r, 'bwa/0/*.sorted.bam'))
        total_sample = len(output_files)
        cmds = ''
        if total_sample <= 100:
            cmds += ('%s mpileup --threads %s -a FMT/AD -q30 -B -Ou -d3000 -f %s %s  | %s call --ploidy 1 -A --threads %s -m > %s' % (
                args.bcf, min(args.th, 40), args.db, ' '.join(output_files), args.bcf, min(args.th, 40), tempbamoutput.replace('.flt.vcf','.raw.vcf')))
            cmds += ('\n%s filter --threads %s -s LowQual -e \'DP>100\' %s > %s \n' % (
                args.bcf,min(args.th,40), tempbamoutput.replace('.flt.vcf','.raw.vcf'), tempbamoutput))
            cmds += ('\n%s view -v snps %s > %s \n' % (
                args.bcf, tempbamoutput, tempbamoutput.replace('.flt.vcf','.flt.snp.vcf')))
        else:
            for i in range(1,int(total_sample/100)):
                cmds += ('%s mpileup --threads %s -a FMT/AD -q30 -B -Ou -d3000 -f %s %s | %s call --ploidy 1 -A --threads %s -m > %s.%s' % (
                        args.bcf, min(args.th, 40), args.db, ' '.join(output_files[(i-1)*100:(i*100)]), args.bcf, min(args.th, 40),
                         tempbamoutput.replace('.flt.vcf', '.raw.vcf'),i))
                cmds += ('\n%s filter --threads %s -s LowQual -e \'DP>=100\' %s.%s > %s.%s \n' % (
                    args.bcf, min(args.th, 40), tempbamoutput.replace('.flt.vcf', '.raw.vcf'), i, tempbamoutput,i))
                cmds += ('\n%s view -v snps %s.%s > %s.%s \n' % (
                    args.bcf, tempbamoutput, i, tempbamoutput.replace('.flt.vcf', '.flt.snp.vcf'),i))
            cmds += (
                        '%s mpileup --threads %s -a FMT/AD -q30 -B -Ou -d3000 -f %s %s | %s call --ploidy 1 -A --threads %s -m > %s.%s' % (
                    args.bcf, min(args.th, 40), args.db, ' '.join(output_files[(i * 100):total_sample]), args.bcf,
                    min(args.th, 40),
                    tempbamoutput.replace('.flt.vcf', '.raw.vcf'), i+1))
            cmds += ('\n%s filter --threads %s -s LowQual -e \'DP>100\' %s.%s > %s.%s \n' % (
                args.bcf, min(args.th, 40), tempbamoutput.replace('.flt.vcf', '.raw.vcf'), i+1, tempbamoutput,i+1))
            cmds += ('\n%s view -v snps %s.%s > %s.%s \n' % (
                args.bcf, tempbamoutput, i+1, tempbamoutput.replace('.flt.vcf', '.flt.snp.vcf'), i+1))
        os.system(cmds)
        #fscripts = open('summary.sh','w')
        #fscripts.write('#!/bin/bash\nsource ~/.bashrc\npy37\n%s'%(''.join(cmds)))
        #fscripts.close()
    # calculate average coverage
    try:
        all_output = open((os.path.join(args.r, 'summary/all.bam.cov.sum.snp.NS.ratio')), 'r')
    except (IOError,FileNotFoundError):
        all_output = open((os.path.join(args.r, 'summary/all.bam.cov.sum.snp.NS.ratio')), 'w')
        all_output.write('metagenome\tgenome\tdepth\tcoverage\tgene_length\tN:S\tobserved_ratio\texpected_ratio\tdNdS\tA-T\tfreq\tN:S\tNSratio'+
                         '\tA-C\tfreq\tN:S\tNSratio\tG-C\tfreq\tN:S\tNSratio\tG-T\tfreq\tN:S\tNSratio\tA-G\tfreq\tN:S\tNSratio\tG-A\tfreq\tN:S\tNSratio\n')
        # calculate allel frequency
        output_files = glob.glob(os.path.join(args.r, 'bwa/0/*.flt.snp.vcf'))
        i = 0
        print('%s summarizing allele frequency and SNPs for %s' % (datetime.now(), os.path.join(args.r, 'bwa/0')))
        for vcf_file in output_files:
            cov_file = vcf_file.replace('.flt.snp.vcf', '.raw.vcf')
            i += 1
            if i%100 == 0:
                print('%s summarizing allele frequency and SNPs for %s files' % (datetime.now(), i))
            freq_call(vcf_file, cov_file)
        all_output.close()
    print('%s start calculate SNP dynamics into %s' % (datetime.now(),
                                                       os.path.join(args.r, 'summary/all.snp.dynamics.txt')))
    try:
        foutput = open(os.path.join(args.r, 'summary/all.snp.dynamics.txt'), 'r')
    except (IOError,FileNotFoundError):
        foutput = open(os.path.join(args.r, 'summary/all.snp.dynamics.txt'), 'w')
        foutput.write('reference\tsample_number\ttotal_snps\n')
        output_files = glob.glob(os.path.join(args.r, 'bwa/0/*.flt.snp.vcf'))
        SNP_dynamics(output_files, foutput)
        foutput.close()

