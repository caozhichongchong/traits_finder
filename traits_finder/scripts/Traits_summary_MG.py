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
                        metavar='/scratch/users/anniz44/scripts/1MG/metadata/all_MGD_GMD_metagenome.metadata.txt')
parser.add_argument("--taxa",
                      help="metadata  of genomes", type=str,
                      default='None',
                      metavar='/scratch/users/anniz44/scripts/1MG/metadata/GTDB_taxon_CG_GMC.brief.habitat.species.selected.all')
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
                    default='None',
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
# Set up cutoff
Diff_gene_cutoff = 0.05 #50 SNP per 1kb

################################################### new class #########################################################
__metaclass__ = type


class SNP_gene:
    # create a class to store SNP_gene
    'a class to store SNP_gene'

    def init(self, gene):
        self.gene = gene
        self.position = dict()
        self.position.setdefault(gene,set())
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
        self.mutposition = dict()
        self.mutposition.setdefault(gene, set())
        self.mutpositiontype = []
        self.NSratio = [0,0]
        self.protein = ''
        self.minor_freq = []

    def addmutposition(self, gene, position,mutation):
        self.mutposition.setdefault(gene, set())
        self.mutposition[gene].add(position)
        self.mutpositiontype.append(mutation)

    def deletemutposition(self, gene,position):
        self.mutposition.setdefault(gene, set())
        if position == 0 :
            self.mutposition.pop(gene, None)
        else:
            self.mutposition[gene].discard(position)

    def addposition(self, gene, position,depth):
        self.position.setdefault(gene, set())
        self.position[gene].add(position)
        self.cov += 1
        self.depth.append(depth)

    def deleteposition(self, gene, position):
        if position == 0 :
            self.position.pop(gene, None)
        else:
            self.position[gene].discard(position)


    def addprotein(self, aa):
        self.protein += aa

    def mutpositioncal(self):
        self.mutpositionsum1 = {0: 0,
                         1: 0,
                         2: 0}
        self.mutpositionsum2 = {0: 0,
                                1: 0,
                                2: 0}
        self.mutpositionsum = 0
        for Chr in self.mutposition:
            self.mutposition2 = list(self.mutposition[Chr])
            self.mutposition2.sort()
            Total = len(self.mutposition2)
            self.mutpositionsum += Total
            if Total >= 1:
                self.mutpositionsum1[abs(self.mutposition2[0]) % 3] += 1
                self.mutpositionsum2[(self.mutposition2[-1]) % 3] += 1
                if Total != 1:
                    for i in range(0,len(self.mutposition2)-1):
                        self.mutpositionsum1[abs(self.mutposition2[i+1] - self.mutposition2[i]) % 3] += 1
                        self.mutpositionsum2[(self.mutposition2[i]) % 3] += 1

    def addSNP_pair(self, pair, position, count, unique_snp_count,depth = 1):
        self.SNP_pair[pair][position] += unique_snp_count
        self.SNP_pair_freq[pair] += unique_snp_count
        self.NSratio[position] += unique_snp_count

    def addalt(self,AllALT_frq):
        self.minor_freq.append(AllALT_frq)

    def deleteSNP_pair(self,SNP_gene):
        for position in [0,1]:
            self.NSratio[position] -= SNP_gene.NSratio[position]
            for pair in SNP_gene.SNP_pair:
                self.SNP_pair[pair][position] -= SNP_gene.SNP_pair[pair][position]
                if position == 0:
                    self.SNP_pair_freq[pair] -= SNP_gene.SNP_pair_freq[pair]

    def sum_SNP_pair(self):
        self.SNP_pair_sum = {'A-T': [0, 0],
                         'A-C': [0, 0],
                         'G-C': [0, 0],
                         'G-T': [0, 0],
                         'A-G': [0, 0],
                         'G-A': [0, 0]}
        self.NSsum = self.NSratio[0] + self.NSratio[1]
        for pair in self.SNP_pair:
            self.SNP_pair_sum[pair][0] += self.SNP_pair[pair][0]
            self.SNP_pair_sum[pair][1] += self.SNP_pair[pair][1]

    def dN_dS(self,SNP_gene_all,normalize=0):
        self.expectNSratio = 'No_expect'
        expectNSratio = [0, 0]
        refSNP_pair_sum_all = Ref_NSratio.get(self.gene, 'None')
        if refSNP_pair_sum_all != 'None':
            refSNP_pair_sum = refSNP_pair_sum_all[0]
            for pair in refSNP_pair_sum:
                if normalize == 0:
                    expectNSratio[0] += refSNP_pair_sum[pair][0] * self.SNP_pair_freq[pair]
                    expectNSratio[1] += refSNP_pair_sum[pair][1] * self.SNP_pair_freq[pair]
                else:
                    expectNSratio[0] += refSNP_pair_sum[pair][0] * SNP_gene_all.SNP_pair_freq[pair]
                    expectNSratio[1] += refSNP_pair_sum[pair][1] * SNP_gene_all.SNP_pair_freq[pair]
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
        self.dNdS = self.NSratiosum
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
# load metadata
def load_meta(metadata):
    Meta = dict()
    i = 0
    Keytitles = ['species',
            'habitat',
                 'selected']
    Col_num = []
    for lines in open(metadata,'r'):
        lines = lines.replace('\n', '').replace('\r', '')
        line_set = lines.split('\t')
        if i == 0:
            for j in range(1,len(line_set)):
                title = lines.split('\t')[j].lower()
                if title in Keytitles:
                    Col_num.append(j)
        else:
            sampleID = line_set[0]
            meta_data = []
            for numbers in Col_num:
                meta_data.append(line_set[numbers].replace('.','_'))
            Meta.setdefault(sampleID,'.'.join(meta_data))
        i+=1
    return Meta


def group_samples(output_files,Meta):
    for files in output_files:
        file_name = os.path.split(files)[-1].split(args.fa)[0].split('.fastq')[0].split('.fasta')[0].split('_1')[0].split('_2')[0]
        if file_name not in Meta:
            for samples in Meta:
                if samples in file_name:
                    file_name = samples
                    break
        if file_name in Meta:
            Sample_group.setdefault(files,'')
            Sample_group[files] += Meta[file_name]
        else:
            print('sample %s (%s) not found in metadata' %(files, file_name))


# Set DNA translate to protein
def translate(seq):
    seq = Seq(seq)
    try:
        return seq.translate()
    except ValueError:
        try:
            return seq.translate(seq.complement())
        except ValueError:
            return ['None']


def NorS(amino1, amino2):
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
    Total = int(len(record_seq)/3)
    temp_SNP_gene = SNP_gene()
    temp_SNP_gene.init(record_name)
    for i in range(0, (Total - position)):
        codon = record_seq[(i * 3 + position):((i + 1) * 3 + position)]
        try:
            codon_NSratio = codontable_NSratio[codon]
            temp_SNP_gene.addprotein(codontable[codon])
            for pair in codon_NSratio.SNP_pair:
                temp_SNP_gene.addSNP_pair(pair, 0, codon_NSratio.SNP_pair[pair][0],codon_NSratio.SNP_pair[pair][0])
                temp_SNP_gene.addSNP_pair(pair, 1, codon_NSratio.SNP_pair[pair][1],codon_NSratio.SNP_pair[pair][1])
        except KeyError:
            pass
    temp_SNP_gene.sum_SNP_pair()
    return [temp_SNP_gene.SNP_pair_sum,temp_SNP_gene.protein,position]


def expectNS(record_name,record_seq):
    Total = int(len(record_seq) / 3)
    temp_result = expectNSsub(record_name, record_seq)
    if len(temp_result[1]) < 0.8 * Total:
        temp_result = expectNSsub(record_name, record_seq,1)
        if len(temp_result[1]) < 0.8 * Total:
            temp_result = expectNSsub(record_name, record_seq,2)
            if len(temp_result[1]) < 0.8 * Total:
                return 'None'
    return temp_result


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


def freq_call_sub(vcf_file,all_wrong_gene,vcf_calculate=False,SNP_gene_list=set()):
        SNP = dict()
        SNP_count = dict()
        SNP_count_genome_count = dict()
        SNP_type = dict()
        SNP_noref = dict()
        Map_length = dict()
        Total = 0
        SNP_gene_list_fixed_position = set()
        for lines in open(vcf_file, 'r'):
            if not lines.startswith("#"):
                lines_set = lines.split('\t')
                if Total == 0:
                    Total = len(lines_set)-9
                Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
                if Depth >= 1:
                    Chr = lines_set[0]
                    if (Chr in SNP_gene_list or vcf_calculate) \
                            and "INDEL" not in lines_set[7] \
                            and (lines_set[6] != 'LowQual' or Depth >= 100):
                        # Depth >= 10 for genome or Depth >= 100 for metagenomes
                        # new SNPs on Chr
                        SNP_type.setdefault(Chr, [''] * Total)
                        SNP_count.setdefault(Chr, 0)
                        SNP_noref.setdefault(Chr, 0)
                        Map_length.setdefault(Chr, 0)
                        Position = int(lines_set[1])
                        Chr_position = '%s--%s' % (Chr, Position)
                        Allels_frq = [0, 0, 0, 0]
                        SNP_count[Chr] += 1
                        genome_ID = 0
                        allels_set = [lines_set[3]]
                        if '.' not in lines_set[4]:
                            allels_set += lines_set[4].split(',')
                        Total_alleles = len(allels_set)
                        SNP_count_genome_count.setdefault(Chr_position, [[0] * Total_alleles, ''])
                        for Subdepth_all in lines_set[9:]:
                            Subdepth = Subdepth_all.split(':')[-1].replace('\n', '')
                            Subdepth_set = Subdepth.split(',')
                            total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_set)
                            #if vcf_calculate:
                            #    SNP_type[Chr][genome_ID] += Subdepth
                            Map_length[Chr] += 1
                            if total_sub_depth > 0:
                                if int(Subdepth_set[0]) == 0:
                                    # not reference ALT
                                    SNP_noref[Chr] += 1
                                for num_allels in range(0, Total_alleles):
                                    allels = allels_set[num_allels]
                                    Subdepth_alleles = int(Subdepth_set[num_allels])
                                    if vcf_calculate:
                                        # count with REF, count with ALT, genome ID with REF, genome ID with ALT, genotype
                                        if Subdepth_alleles >= 1:
                                            SNP_count_genome_count[Chr_position][0][num_allels] += 1
                                            SNP_count_genome_count[Chr_position][-1] += allels
                                            SNP_type[Chr][genome_ID] += '1'
                                    if allels in Allels:
                                        Allels_frq[Allels[allels]] += Subdepth_alleles
                                    else:
                                        pass
                            genome_ID += 1
                            if vcf_calculate:
                                SNP_count_genome_count[Chr_position][-1] += '\t'
                        if vcf_calculate:
                            # check any minor ALT, vcf snp
                            Depth_position_snp = sum(Allels_frq)
                            SNP.setdefault(Chr_position, Allels_frq)
                            if not any(Depth_position_snp == allels for allels in Allels_frq):
                                # not fixed SNPs
                                SNP_gene_list.add(Chr)
                            else:
                                SNP_gene_list_fixed_position.add(Chr_position)
                        else:
                            # vcf coverage
                            SNP.setdefault(Chr_position, Allels_frq)
        if vcf_calculate:
            SNP_type_sum = dict()
            for Chr in SNP_type:
                SNP_type_Chr_all = SNP_type[Chr]
                if SNP_type_Chr_all!=[]:
                    SNP_type_Chr_uniq = set(SNP_type_Chr_all)
                    SNP_type_chr_cal = []
                    for SNP_type_Chr_sub in SNP_type_Chr_uniq:
                        SNP_type_chr_cal.append(SNP_type_Chr_all.count(SNP_type_Chr_sub))
                    SNP_type_chr_cal.sort(reverse=True)
                    SNP_type_sum.setdefault(Chr,':'.join(str(numbers) for numbers in SNP_type_chr_cal))
            return [SNP,SNP_count, Total, SNP_count_genome_count,SNP_type_sum, SNP_gene_list,SNP_gene_list_fixed_position]
        else:
            # delete wrong genes, more than 50 non-reference per 1kb
            for Chr in SNP_noref:
                if SNP_noref[Chr] / Map_length[Chr] >= Diff_gene_cutoff:
                    all_wrong_gene.append(Chr)
            return [SNP, SNP_count, Total]


def dN_dS_ratio(seq1,seq2,method = 'ML'):
    if method=='ML':
        return cal_dn_ds(CodonSeq(seq1), CodonSeq(seq2), method='ML')
    elif method == 'NG86':
        return cal_dn_ds(CodonSeq(seq1), CodonSeq(seq2), method='NG86')
    elif method == 'LWL85':
        return cal_dn_ds(CodonSeq(seq1), CodonSeq(seq2), method='LWL85')
    elif method == 'YN00':
        return cal_dn_ds(CodonSeq(seq1), CodonSeq(seq2), method='YN00')


def filtersnp(SNP_gene_temp,SNP_gene_all,Chr,all_output_filter_list,sample_name):
    Gene_length = 1000
    if Chr in Mapping:
        Gene_length = Mapping[Chr]
    if Chr == 'all':
        Gene_length = 0
        for allchr in Mapping:
            Gene_length += Mapping[allchr]
    SNP_gene_temp.cov = float(SNP_gene_temp.cov) / float(Gene_length)
    SNP_gene_temp.mutpositioncal()
    if SNP_gene_temp.cov >= 0.8 or Chr == 'all':
        temp_depth = SNP_gene_temp.depth
        temp_depth.sort()
        temp_std = statistics.stdev(temp_depth)
        temp_mean = statistics.mean(temp_depth)
        temp_range = [temp_mean - 1.5 * temp_std, temp_mean + 1.5 * temp_std]
        temp_depth_filter = []
        for position in range(0, len(temp_depth)):
            a = temp_depth[position]
            if (temp_range[0] <= a and a <= temp_range[1]):
                temp_depth_filter.append(a)
            else:
                pass
        SNP_gene_temp.depth = statistics.mean(temp_depth_filter)
    if Chr != 'all':
        if SNP_gene_temp.cov < 0.8:
            # not enough coverage
            SNP_gene_all.deleteposition(Chr,0)
            SNP_gene_all.deletemutposition(Chr, 0)
            SNP_gene_all.deleteSNP_pair(SNP_gene_temp)
            all_output_filter_list.append('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (
            sample_name, Chr, Gene_length, SNP_gene_temp.cov, SNP_gene_temp.mutpositionsum, 'not covered'))
            return [Chr,'lowcoverage']
        if SNP_gene_temp.mutpositionsum >= Gene_length*Diff_gene_cutoff:
            # a different gene
            SNP_gene_all.deletemutposition(Chr,0)
            SNP_gene_all.deleteSNP_pair(SNP_gene_temp)
            all_output_filter_list.append('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (
                sample_name, Chr, Gene_length, SNP_gene_temp.cov, SNP_gene_temp.mutpositionsum, 'different gene'))
            return [Chr,'diffgene']
        return ['','']


def sum_snp(SNP_gene_temp, SNP_gene_all, Chr, sample_name, Total,SNP_type_chr_cal):
    N_temp = SNP_gene_temp.NSratio[0]
    S_temp = SNP_gene_temp.NSratio[1]
    N_S_sum = N_temp + S_temp
    Gene_length = 1000
    if Chr in Mapping:
        Gene_length = Mapping[Chr]
    if Chr == 'all':
        Gene_length = 0
        for allchr in SNP_gene_all.position:
            Gene_length += Mapping.get(allchr,0)
    new_line = '%s\t%s\t%s\t%s\t%s\t%s' % (sample_name, Chr, Total,
                                                                   '%.3f' % (SNP_gene_temp.depth),
                                                                   '%.3f' % (SNP_gene_temp.cov), Gene_length)
    if N_S_sum > 0:
        minor_feq = statistics.mean(SNP_gene_temp.minor_freq)
        new_line += '\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (N_S_sum,
                                                             SNP_gene_temp.mutpositionsum,
                                                             '%.5f' % (minor_feq),
                                                                  SNP_gene_temp.mutpositiontype.count('non_fixed'),
                                                                  SNP_gene_temp.mutpositiontype.count('fixed'),
                                                                  SNP_type_chr_cal,
                                                             SNP_gene_temp.mutpositionsum1[0],
                                                             SNP_gene_temp.mutpositionsum1[1],
                                                             SNP_gene_temp.mutpositionsum1[2],
                                                             SNP_gene_temp.mutpositionsum2[0],
                                                             SNP_gene_temp.mutpositionsum2[1],
                                                             SNP_gene_temp.mutpositionsum2[2])
        # cancel dNdS
        #SNP_gene_temp.dN_dS(SNP_gene_all,1) # normalized
        try:
            SNP_gene_temp.NSratiosum = SNP_gene_temp.NSratio[0] / SNP_gene_temp.NSratio[1]
        except (ZeroDivisionError, TypeError):
            if SNP_gene_temp.NSratio[0] == 0:
                SNP_gene_temp.NSratiosum = 'observe_None'
            else:
                SNP_gene_temp.NSratiosum = 'observe_N_only'
        new_line += ('\t%s:%s\t%s' % (N_temp, S_temp, SNP_gene_temp.NSratiosum))
        #new_line += '\t%s\t%s' % (SNP_gene_temp.expectNSratio, SNP_gene_temp.dNdS)
        #SNP_gene_temp.dN_dS(SNP_gene_all,0) # non-nomalized
        #new_line += '\t%s\t%s' % (SNP_gene_temp.expectNSratio, SNP_gene_temp.dNdS)
        for pair in SNP_gene_temp.SNP_pair:
            pair_freq = SNP_gene_temp.SNP_pair_freq[pair]#/SNP_gene_temp.mutpositionsum/SNP_gene_temp.depth
            pair_N = SNP_gene_temp.SNP_pair[pair][0]
            pair_S = SNP_gene_temp.SNP_pair[pair][1]
            new_line += ('\t%s\t%s:%s' % ('%.5f'%pair_freq, pair_N, pair_S))
    new_line += '\n'
    return new_line


def freq_call(vcf_file,cov_file,Not_pass=dict()):
    Output =  dict()
    Output2 = dict()
    Output3 = dict()
    all_SNP_gene_temp = dict()
    all_wrong_gene = []
    # global mutation frequency
    Chr = 'all'
    SNP_gene_all = SNP_gene()
    SNP_gene_all.init(Chr)
    print('%s summarize SNPs for vcf file %s' % (datetime.now(), vcf_file))
    try:
        # only output genes that have SNPs
        SNP, SNP_count, Total,SNP_count_genome_count,SNP_type_sum,SNP_gene_list,SNP_gene_list_fixed_position_all = freq_call_sub(vcf_file,all_wrong_gene,True,set())
        SNP_all, SNP_count_all, Total = freq_call_sub(cov_file,all_wrong_gene,False,SNP_gene_list)
        print('%s deleting non-qualified genes %s' % (datetime.now(), ' '.join(all_wrong_gene)))
        sample_name = os.path.split(cov_file)[1].split('.fasta')[0].split('.fastq')[0].split(args.fa)[0]
        #Depth = dict()
        for Chr_position in SNP_all:
            Depth_position = sum(SNP_all[Chr_position])
            Chr = Chr_position.split('--')[0]
            if Chr not in all_wrong_gene:
                if Chr not in all_SNP_gene_temp:
                    SNP_gene_temp = SNP_gene()
                    SNP_gene_temp.init(Chr)
                    all_SNP_gene_temp.setdefault(Chr,SNP_gene_temp)
                    Output.setdefault(Chr,[])
                    Output2.setdefault(Chr, [])
                    Output3.setdefault(Chr, [])
                SNP_gene_temp = all_SNP_gene_temp[Chr]
                #N_S.setdefault(Chr, [0, 0])
                position = Chr_position.split('--')[1]
                SNP_gene_temp.addposition(Chr,position, Depth_position)
                SNP_gene_all.addposition(Chr,position, Depth_position)
                lines = '%s\t%s\t%s\t%s\t%s\t%s\t%s' % (Chr,position,
                        Depth_position, SNP_all[Chr_position][0],
                        SNP_all[Chr_position][1],
                        SNP_all[Chr_position][2],
                        SNP_all[Chr_position][3])
                Output[Chr].append(lines + '\n')
                if Chr_position in SNP:
                    mutation = 'non_fixed'
                    position = int(position)
                    if Chr_position in SNP_gene_list_fixed_position_all:
                        mutation = 'fixed'
                    # count how many genomes have a SNP
                    Chr_position_genome_count = SNP_count_genome_count.get(Chr_position, [[0], ''])
                    # a SNP
                    # cancel dNdS
                    # calculate N or S
                    # refSNP_pair_sum_all = Ref_NSratio.get(Chr, 'None')
                    refSNP_pair_sum_all = 'None'
                    # calculate Allele frequency
                    Major_ALT, Minor_ALT = ALT_freq(SNP[Chr_position])
                    REF = Major_ALT[0]
                    Ref_frq = Major_ALT[1]
                    AllALT_frq = (Depth_position - Ref_frq)/Depth_position
                    lines += '\t%s' % (AllALT_frq)
                    SNP_gene_temp.addalt(AllALT_frq)
                    SNP_gene_all.addalt(AllALT_frq)
                    if refSNP_pair_sum_all != 'None':
                        #  observed NS ratio calculated
                        refSNP_condon_start = refSNP_pair_sum_all[-1]
                        codon_start = position - 1 - int((position - 1) % 3) + refSNP_condon_start
                        if codon_start <= position - 1:
                            Ref_seq_chr = Ref_seq[Chr]
                            SNP_seq_chr = Ref_seq_chr
                            Ref_seq_chr = causeSNP(Ref_seq_chr, position, REF)
                            Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                            if len(Ref_seq_codon) == 3:
                                Ref_seq_aa = translate(Ref_seq_codon)[0]
                                for minor in Minor_ALT:
                                    ALT = minor[0]
                                    ALT_frq = minor[1]
                                    SNP_seq_chr = causeSNP(SNP_seq_chr, position, ALT)
                                    SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                                    SNP_seq_aa = translate(SNP_seq_codon)[0]
                                    temp_NorS = NorS(Ref_seq_aa, SNP_seq_aa)
                                    SNP_pair = transitions(REF, ALT)
                                    SNP_gene_temp.addSNP_pair(SNP_pair, N_S_set[temp_NorS],
                                                              ALT_frq, 1, Depth_position)
                                    SNP_gene_temp.addmutposition(Chr, position,mutation)
                                    lines += '\t%s\t%s\t%s:%s:%s\t%s:%s:%s\t%s\t%s' % (mutation,
                                                                                       '%s' % (
                                                                                       Chr_position_genome_count[0]),
                                                                                       REF, Ref_frq, Ref_seq_aa,
                                                                                       ALT, ALT_frq, SNP_seq_aa,
                                                                                       temp_NorS,
                                                                                       Chr_position_genome_count[-1])
                                    SNP_gene_all.addSNP_pair(SNP_pair, N_S_set[temp_NorS],
                                                             ALT_frq, 1, Depth_position)
                                    SNP_gene_all.addmutposition(Chr, position,mutation)
                    else:
                        codon_start = position - 1 - int((position - 1) % 3)
                        if codon_start <= position - 1:
                            Ref_seq_chr = Ref_seq[Chr]
                            SNP_seq_chr = Ref_seq_chr
                            Ref_seq_chr = causeSNP(Ref_seq_chr, position, REF)
                            Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                            if len(Ref_seq_codon) == 3:
                                Ref_seq_aa = translate(Ref_seq_codon)[0]
                                for minor in Minor_ALT:
                                    ALT = minor[0]
                                    ALT_frq = minor[1]
                                    SNP_seq_chr = causeSNP(SNP_seq_chr, position, ALT)
                                    SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                                    SNP_seq_aa = translate(SNP_seq_codon)[0]
                                    temp_NorS = NorS(Ref_seq_aa, SNP_seq_aa)
                                    SNP_pair = transitions(REF, ALT)
                                    SNP_gene_temp.addSNP_pair(SNP_pair, N_S_set[temp_NorS],
                                                              ALT_frq, 1, Depth_position)
                                    SNP_gene_temp.addmutposition(Chr, position,mutation)
                                    lines += '\t%s\t%s\t%s:%s:%s\t%s:%s:%s\t%s\t%s' % (mutation,
                                        '%s' % (Chr_position_genome_count[0]),
                                        REF, Ref_frq, Ref_seq_aa,
                                        ALT, ALT_frq, SNP_seq_aa, temp_NorS, Chr_position_genome_count[-1])
                                    SNP_gene_all.addSNP_pair(SNP_pair, N_S_set[temp_NorS],
                                                             ALT_frq, 1, Depth_position)
                                    SNP_gene_all.addmutposition(Chr, position,mutation)
                    Output2[Chr].append(lines + '\n')
                    Output3[Chr].append('%s\t%s\t%s\t%s\n' % (SNP[Chr_position][0],
                                                              SNP[Chr_position][1],
                                                              SNP[Chr_position][2],
                                                              SNP[Chr_position][3]))
        # global mutation frequency
        if all_SNP_gene_temp!= dict():
            all_output_list = []
            all_output_list2 = []
            all_output_filter_list = []
            for Chr in all_SNP_gene_temp:
                SNP_gene_temp = all_SNP_gene_temp[Chr]
                filter_result = filtersnp(SNP_gene_temp, SNP_gene_all, Chr,all_output_filter_list,sample_name)
                Not_pass.setdefault(filter_result[0],filter_result[1])
            Not_pass.pop('', None)
            for Chr in Not_pass:
                all_SNP_gene_temp.pop(Chr, None)
                Output2.pop(Chr, None)
                Output3.pop(Chr, None)
                if Not_pass[Chr]=='lowcoverage':
                    Output.pop(Chr, None)
            print('%s deleting non-qualified genes %s'%(datetime.now(),' '.join(Not_pass)))
            filtersnp(SNP_gene_all, SNP_gene_all, 'all',all_output_filter_list,sample_name)
            new_line = sum_snp(SNP_gene_all, SNP_gene_all, 'all', sample_name, Total,SNP_type_sum.get('all',''))
            all_output_list.append(new_line)
            all_output_list2.append(new_line)
            Output_list = []
            Output2_list = []
            Output3_list = []
            for Chr in all_SNP_gene_temp:
                SNP_gene_temp = all_SNP_gene_temp[Chr]
                new_line = sum_snp(SNP_gene_temp,SNP_gene_all,Chr,sample_name,Total,SNP_type_sum.get(Chr,''))
                all_output_list.append(new_line)
                if SNP_gene_temp.mutpositionsum > 0:
                    all_output_list2.append(new_line)
                Output_list+=Output[Chr]
                Output2_list += Output2[Chr]
                Output3_list += Output3[Chr]
            all_output.write(''.join(all_output_list))
            all_output2.write(''.join(all_output_list2))
            all_output_filter.write(''.join(all_output_filter_list))
            Output_list = sorted(Output_list, key=operator.itemgetter(1))
            Output_list = sorted(Output_list, key=operator.itemgetter(0))
            Output2_list = sorted(Output2_list, key=operator.itemgetter(1))
            Output2_list = sorted(Output2_list, key=operator.itemgetter(0))
            foutput = open(cov_file + '.frq', 'w')
            foutput.write('#CHR\tPOS\tDEP\tA\tT\tG\tC\n')
            foutput.write(''.join(Output_list))
            foutput.close()
            foutput = open(cov_file + '.frq.snp', 'w')
            foutput.write('#CHR\tPOS\tDEP\tA\tT\tG\tC\tALT_frq\tmutation_type\tgenome_REF_ALT_count\tREF:count:REF_aa\tALT:count:ALT_aa\tN_or_S\tgenotype_allgenomes\n')
            foutput.write(''.join(Output2_list))
            foutput.close()
            #foutput = open(cov_file + '.frq.snp.clean', 'w')
            #foutput.write('#A\tT\tG\tC\n')
            #foutput.write(''.join(Output3_list))
            #foutput.close()
    except (IOError, FileNotFoundError):
        pass
    if args.strainfinder!= 'None':
        strain_finder(cov_file + '.frq.snp')
    return Not_pass


def freq_call_sub_dynamics(vcf_file):
    SNP_count = dict()
    Set_subsample = dict()
    print('%s calculate SNP dynamics for %s by 10 iterations' % (datetime.now(), vcf_file))
    # get total number of samples
    for lines in open(vcf_file, 'r'):
        Total = 0
        if not lines.startswith("#"):
            lines_set = lines.split('\t')
            if Total == 0:
                Total =  len(lines_set)
                Set1 = range(1, Total-9)
                Set2 = range(9, Total)
                break
    for i in Set1:
        for repeat_time in range(0, 10):
            temp_list = set()
            while len(temp_list) < i:
                temp_list.add(random.choice(Set2))
            Set_subsample.setdefault(i,set())
            Set_subsample[i].add(temp_list)
    Set_subsample.setdefault(Total-9, [Set2])
    for lines in open(vcf_file, 'r'):
        if not lines.startswith("#"):
            lines_set = lines.split('\t')
            Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
            Chr = lines_set[0]
            if "INDEL" not in lines_set[7] and (lines_set[6] != 'LowQual' or Depth >= 100):
                # Depth >= 10 for genome or Depth >= 100 for metagenomes
                # new SNPs on Chr
                for i in Set_subsample:
                    Chr_i = '%s\t%s' % (Chr, i)
                    SNP_count.setdefault(Chr_i, 0)
                    for subset in Set_subsample[i]:
                        temp_10 = 0
                        for j in subset:
                            for Subdepth_all in lines_set[j]:
                                Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                                if sum(int(Subdepth_sub) for Subdepth_sub in Subdepth) > 0:
                                    # one new SNP
                                    temp_10 += 1
                                    break
                        if temp_10 > 10:
                            print('wrong calculation')
                    SNP_count[Chr_i] += temp_10
    return SNP_count


def SNP_dynamics(output_files, SNP_outputfile):
    SNP_dynamics_pair = dict()
    SNP_dynamics_pair_output = []
    for files in output_files:
        metadata = os.path.split(files)[1].split('.flt.snp.vcf')[0]
        SNP_count = freq_call_sub_dynamics(files)
        print('%s output calculate SNP dynamics'
              % (datetime.now()))
        for Chr_i in SNP_count:
            metadata_Chr_i = '%s\t%s'%(metadata,Chr_i)
            SNP_dynamics_pair.setdefault(metadata_Chr_i, [0,0]) #SNP number, subvcffile number
            SNP_dynamics_pair[metadata_Chr_i][0] += SNP_count[Chr_i]
            SNP_dynamics_pair[metadata_Chr_i][1] += 1
    for metadata_Chr_i in SNP_dynamics_pair:
        SNP_dynamics_pair_output.append('%s\t%s\n' % (metadata_Chr_i,
                                                      SNP_dynamics_pair[metadata_Chr_i][0]/SNP_dynamics_pair[metadata_Chr_i][1]))
    SNP_outputfile.write(''.join(SNP_dynamics_pair_output))


################################################### Programme #######################################################
# load database length file
Mapping = dict()
reference_database = os.path.split(args.db)[-1]
print('reference database set as %s'%(reference_database))
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
    print('%s calculate expected N S ratio for %s' % (datetime.now(), args.db))
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
                    temp_NorS = NorS(translate(codon)[0], translate(new_codon)[0])
                    SNP_pair = transitions(REF, ALT)
                    SNP_gene_temp.addSNP_pair(SNP_pair, N_S_set[temp_NorS], 1, 1)
    # load database seq
    Ref_seq = dict()
    Ref_NSratio = dict()
    # cancel dNdS
    for record in SeqIO.parse(args.db, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        Ref_seq.setdefault(record_id,record_seq)
    #    Ref_NSratio.setdefault(record_id,
    #                           expectNS(record_id,record_seq))
    #foutput = open(args.db+'.ref.NS.ratio','w')
    #foutput_list = []
    #for Ref in Ref_NSratio:
    #    foutput_list.append('%s\t%s\t\n'%(Ref,Ref_NSratio[Ref]))
    #foutput.write(''.join(foutput_list))
    #foutput.close()
    print('%s output reference NS ratio into %s' % (datetime.now(),args.db+'.ref.NS.ratio'))
    print('%s load metadata' % (datetime.now()))
    # load metadata
    if args.meta != 'None':
        Meta_meta = load_meta(args.meta)
    if args.taxa != 'None':
        Taxa_meta = load_meta(args.taxa)
    # all bam file call vcf
    tempbamoutput = glob.glob(os.path.join(args.r, 'summary/vcf/%s*.flt.snp.vcf*')%(reference_database))
    print('%s merge files and call SNPs %s' % (datetime.now(), 'summary/vcf/%s*.flt.snp.vcf*')%(reference_database))
    if tempbamoutput == [] and args.meta != 'None' and args.taxa != 'None':
        output_files = glob.glob(os.path.join(args.r, 'bwa/0/*%s*.sorted.bam')%(args.fa))
        # grouping
        Sample_group = dict()
        Group_sample = dict()
        if args.meta != 'None':
            group_samples(output_files,Meta_meta)
        if args.taxa != 'None':
            group_samples(output_files,Taxa_meta)
        if Sample_group == dict():
            Group_sample.setdefault('all',output_files)
        else:
            for files in Sample_group:
                groups = Sample_group[files]
                Group_sample.setdefault(groups,[])
                Group_sample[groups].append(files)
        cmds = ''
        for groups in Group_sample:
            sub_samples = Group_sample[groups]
            total_sample = len(sub_samples)
            tempbamoutput = os.path.join(args.r, 'summary/vcf/%s.%s.flt.snp.vcf'%(reference_database,groups))
            print('%s merge files and call SNPs %s' % (datetime.now(), os.path.join(args.r, 'summary/vcf/%s.%s.flt.snp.vcf'%(reference_database,groups))))
            if total_sample <= 100:
                cmds += ('%s mpileup --threads %s -a FMT/AD -q30 -B -Ou -d3000 -f %s %s  | %s call --ploidy 1 -A --threads %s -m > %s' % (
                    args.bcf, min(args.th, 40), args.db, ' '.join(sub_samples), args.bcf, min(args.th, 40), tempbamoutput.replace('.flt.snp.vcf','.raw.vcf')))
                cmds += ('\n%s filter --threads %s -s LowQual -e \'DP>100\' %s > %s \n' % (
                    args.bcf,min(args.th,40), tempbamoutput.replace('.flt.snp.vcf','.raw.vcf'),
                    tempbamoutput.replace('.flt.snp.vcf','.flt.vcf')))
                cmds += ('\n%s view -v snps %s > %s \nrm -rf %s\n' % (
                    args.bcf, tempbamoutput.replace('.flt.snp.vcf','.flt.vcf'), tempbamoutput,
                    tempbamoutput.replace('.flt.snp.vcf','.flt.vcf')))
            else:
                for i in range(1,int(total_sample/100)):
                    cmds += ('%s mpileup --threads %s -a FMT/AD -q30 -B -Ou -d3000 -f %s %s | %s call --ploidy 1 -A --threads %s -m > %s.%s' % (
                            args.bcf, min(args.th, 40), args.db, ' '.join(sub_samples[(i-1)*100:(i*100)]), args.bcf, min(args.th, 40),
                             tempbamoutput.replace('.flt.snp.vcf', '.raw.vcf'),i))
                    cmds += ('\n%s filter --threads %s -s LowQual -e \'DP>=100\' %s.%s > %s.%s \n' % (
                        args.bcf, min(args.th, 40), tempbamoutput.replace('.flt.snp.vcf', '.raw.vcf'), i,
                        tempbamoutput.replace('.flt.snp.vcf','.flt.vcf'),i))
                    cmds += ('\n%s view -v snps %s.%s > %s.%s \nrm -rf %s.%s\n' % (
                        args.bcf, tempbamoutput.replace('.flt.snp.vcf','.flt.vcf'), i,
                        tempbamoutput,i,
                        tempbamoutput.replace('.flt.snp.vcf', '.flt.vcf'), i))
                cmds += (
                            '%s mpileup --threads %s -a FMT/AD -q30 -B -Ou -d3000 -f %s %s | %s call --ploidy 1 -A --threads %s -m > %s.%s' % (
                        args.bcf, min(args.th, 40), args.db, ' '.join(sub_samples[(i * 100):total_sample]), args.bcf,
                        min(args.th, 40),
                        tempbamoutput.replace('.flt.snp.vcf', '.raw.vcf'), i+1))
                cmds += ('\n%s filter --threads %s -s LowQual -e \'DP>100\' %s.%s > %s.%s \n' % (
                    args.bcf, min(args.th, 40), tempbamoutput.replace('.flt.snp.vcf', '.raw.vcf'), i+1,
                    tempbamoutput.replace('.flt.snp.vcf','.flt.vcf'),i+1))
                cmds += ('\n%s view -v snps %s.%s > %s.%s \nrm -rf %s.%s\n' % (
                    args.bcf, tempbamoutput.replace('.flt.snp.vcf','.flt.vcf'),i+1,
                    tempbamoutput, i+1,
                    tempbamoutput.replace('.flt.snp.vcf','.flt.vcf'),i+1))
        os.system(cmds)
        fscripts = open('%s.summary.sh'%(reference_database),'w')
        fscripts.write('#!/bin/bash\nsource ~/.bashrc\npy37\n%s'%(''.join(cmds)))
        fscripts.close()
    # calculate average coverage
    try:
        all_output = open((os.path.join(args.r, 'summary/%s.all.bam.cov.sum.snp.NS.ratio')%(reference_database)), 'r')
    except (IOError,FileNotFoundError):
        all_output = open((os.path.join(args.r, 'summary/%s.all.bam.cov.sum.snp.NS.ratio')%(reference_database)), 'w')
        all_output2 = open((os.path.join(args.r, 'summary/%s.all.bam.cov.sum.snp.NS.ratio.snp.only') % (reference_database)), 'w')
        all_output_filter = open((os.path.join(args.r, 'summary/%s.all.bam.cov.sum.nopass.gene')%(reference_database)), 'w')
        all_output.write(
                'metagenome\tgenome\tnumber_sample\tdepth\tcoverage\tgene_length\ttotal_SNP\ttotal_SNP_position\tavg_minor_freq\tsequence_type_count\tSNP_dis_0\tSNP_dis_1\tSNP_dis_2\tSNP_codon_0\tSNP_codon_1\tSNP_codon_2\tN:S\tobserved_ratio\texpected_ratio_norm\tdNdS_norm\texpected_ratio\tdNdS\tA-T_freq\tA-T_N:S' +
            '\tA-C_freq\tA-C_N:S\tG-C_freq\tG-C_N:S\tG-T_freq\tG-T_N:S\tA-G_freq\tA-G_N:S\tG-A_freq\tG-A_N:S\n')
        all_output2.write(
                'metagenome\tgenome\tnumber_sample\tdepth\tcoverage\tgene_length\ttotal_SNP\ttotal_SNP_position\tavg_minor_freq\tsequence_type_count\tSNP_dis_0\tSNP_dis_1\tSNP_dis_2\tSNP_codon_0\tSNP_codon_1\tSNP_codon_2\tN:S\tobserved_ratio\texpected_ratio_norm\tdNdS_norm\texpected_ratio\tdNdS\tA-T_freq\tA-T_N:S' +
            '\tA-C_freq\tA-C_N:S\tG-C_freq\tG-C_N:S\tG-T_freq\tG-T_N:S\tA-G_freq\tA-G_N:S\tG-A_freq\tG-A_N:S\n')
        # calculate allele frequency
        output_files = glob.glob(os.path.join(args.r, 'bwa/0/*%s*.flt.snp.vcf'%(args.fa)))
        i = 0
        print('%s summarizing allele frequency and SNPs for %s' % (datetime.now(), os.path.join(args.r, 'bwa/0')))
        ref_notpass = dict()
        for vcf_file in output_files:
            cov_file = vcf_file.replace('.flt.snp.vcf', '.raw.vcf')
            i += 1
            if i%100 == 0:
                print('%s summarizing allele frequency and SNPs for %s files' % (datetime.now(), i))
            freq_call(vcf_file, cov_file, ref_notpass)
        all_output.close()
        all_output2.close()
        all_output_filter.close()
    if args.meta != 'None' and args.taxa != 'None':
        # calculate average coverage
        try:
            all_output = open((os.path.join(args.r, 'summary/%s.all.merged.bam.cov.sum.snp.NS.ratio')%(reference_database)), 'r')
        except (IOError, FileNotFoundError):
            all_output = open((os.path.join(args.r, 'summary/%s.all.merged.bam.cov.sum.snp.NS.ratio')%(reference_database)), 'w')
            all_output2 = open(
                (os.path.join(args.r, 'summary/%s.all.merged.bam.cov.sum.snp.NS.ratio.snp.only') % (reference_database)), 'w')
            all_output_filter = open((os.path.join(args.r, 'summary/%s.all.merged.bam.cov.sum.nopass.gene')%(reference_database)), 'w')
            all_output.write(
                'metagenome\tgenome\tnumber_sample\tdepth\tcoverage\tgene_length\ttotal_SNP\ttotal_SNP_position\tavg_minor_freq\tnon_fixed_count\tfixed_count\tsequence_type_count\tSNP_dis_0\tSNP_dis_1\tSNP_dis_2\tSNP_codon_0\tSNP_codon_1\tSNP_codon_2\tN:S\tobserved_ratio\texpected_ratio_norm\tdNdS_norm\texpected_ratio\tdNdS\tA-T_freq\tA-T_N:S' +
            '\tA-C_freq\tA-C_N:S\tG-C_freq\tG-C_N:S\tG-T_freq\tG-T_N:S\tA-G_freq\tA-G_N:S\tG-A_freq\tG-A_N:S\n')
            all_output2.write(
                'metagenome\tgenome\tnumber_sample\tdepth\tcoverage\tgene_length\ttotal_SNP\ttotal_SNP_position\tavg_minor_freq\tnon_fixed_count\tfixed_count\tsequence_type_count\tSNP_dis_0\tSNP_dis_1\tSNP_dis_2\tSNP_codon_0\tSNP_codon_1\tSNP_codon_2\tN:S\tobserved_ratio\texpected_ratio_norm\tdNdS_norm\texpected_ratio\tdNdS\tA-T_freq\tA-T_N:S' +
            '\tA-C_freq\tA-C_N:S\tG-C_freq\tG-C_N:S\tG-T_freq\tG-T_N:S\tA-G_freq\tA-G_N:S\tG-A_freq\tG-A_N:S\n')
            # calculate allele frequency
            output_files = glob.glob(os.path.join(args.r, 'summary/vcf/*%s*.flt.snp.vcf*'%(args.fa)))
            i = 0
            print('%s summarizing allele frequency and SNPs for %s' % (datetime.now(), os.path.join(args.r, 'summary/vcf/')))
            for vcf_file in output_files:
                cov_file = vcf_file.replace('.flt.snp.vcf', '.raw.vcf')
                i += 1
                if i % 100 == 0:
                    print('%s summarizing allele frequency and SNPs for %s files' % (datetime.now(), i))
                freq_call(vcf_file, cov_file)
            all_output.close()
            all_output2.close()
            all_output_filter.close()
        print('%s start calculate SNP dynamics into %s' % (datetime.now(),
                                                           os.path.join(args.r, 'summary/%s.all.snp.dynamics.txt'%(reference_database))))
        try:
            foutput = open(os.path.join(args.r, 'summary/%s.all.snp.dynamics.txt'%(reference_database)), 'r')
        except (IOError,FileNotFoundError):
            foutput = open(os.path.join(args.r, 'summary/%s.all.snp.dynamics.txt')%(reference_database), 'w')
            foutput.write('habitat\treference\tsample_number\ttotal_snps\n')
            output_files = glob.glob(os.path.join(args.r, 'summary/vcf/*%s*.flt.snp.vcf*'%(args.fa)))
            SNP_dynamics(output_files, foutput)
            foutput.close()

