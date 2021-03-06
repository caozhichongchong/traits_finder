import os
from Bio import SeqIO
import argparse
import glob
from datetime import datetime
import statistics
import random
import subprocess
import random
import operator


############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-t",
                    help="trait name", type=str, default='ARG',metavar='trait name')
parser.add_argument("-m",
                    help="mapping file of traits to function", type=str,
                    default='Butyrate.pro.mapping.txt',
                    metavar='Butyrate.pro.mapping.txt')
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
parser.add_argument("-i",
                        help="input folder of genomes", type=str,
                        default='.',metavar='current dir (.)')
parser.add_argument("-fa",
                        help="input format of genome sequence",
                        type=str, default='.fa', metavar='.fasta, .fna, .fastq or .fa')
# optional input setup
parser.add_argument("--r",
                    help="input directory or folder of your previous results by Traits_WGD.py",
                    type=str, default='Result',metavar='Result')
parser.add_argument("--s",
                    help="input directory or folder of your previous results by traits summary",
                    type=str, default='None',metavar='summary')
parser.add_argument("--meta",
                        help="metadata  of metagenomes", type=str,
                        default='None',
                        metavar='/scratch/users/anniz44/scripts/1MG/metadata/all_MGD_GMD_metagenome.metadata.txt')
parser.add_argument("--taxa",
                      help="metadata  of genomes", type=str,
                      default='None',
                      metavar='/scratch/users/anniz44/scripts/1MG/metadata/GTDB_taxon_CG_GMC.brief.habitat.species.selected.all')
# optional search parameters
parser.add_argument('--g',
                        help="Optional: gene-level HGT finding; --g F (default function-level); --g T for gene-level",
                        metavar=['T', 'F'], action='store', default='F', type=str)
parser.add_argument('--th',
                        help="Optional: set the thread number assigned for running XXX (default 1)",
                        metavar="1 or more", action='store', default=1, type=int)
# requirement for software calling
parser.add_argument('--u','--usearch',
                        help="Optional: use two-step method for blast search,"+
                             " \'None\' for using one step, \'usearch\' for using two-step \
                             (complete path to usearch if not in PATH), (default: \'None\')",
                        metavar="None or usearch",
                        action='store', default='None', type=str)
parser.add_argument('--dm', '--diamond',
                      help="Optional: use two-step method for blast search," +
                           " \'None\' for using one step, \'diamond\' for using two-step \
                           (complete path to diamond if not in PATH), (default: \'None\')",
                      metavar="None or diamond",
                      action='store', default='None', type=str)
parser.add_argument('--hs',
                      help="Optional: use two-step method for blast search," +
                           " \'None\' for using one step, \'hs-blastn\' for using two-step \
                           (complete path to hs-blastn if not in PATH), (default: \'None\')",
                      metavar="None or hs-blastn",
                      action='store', default='None', type=str)
parser.add_argument('--bp',
                    help="Optional: complete path to blastp or blastn if not in PATH,",
                    metavar="/usr/local/bin/blast",
                    action='store', default='blast', type=str)
parser.add_argument('--mf', '--mafft',
                      help="Optional: complete path to mafft if not in PATH,",
                      metavar="/usr/local/bin/mafft",
                      action='store', default='None', type=str)
parser.add_argument('--ft', '--fasttree',
                      help="Optional: complete path to fasttree if not in PATH,",
                      metavar="/usr/local/bin/fasttree",
                      action='store', default='None', type=str)
parser.add_argument('--ani', '--fastani',
                      help="Optional: complete path to fastANI if not in PATH,",
                      metavar="/usr/local/bin/fastANI",
                      action='store', default='None', type=str)
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


################################################## Definition ########################################################
args = parser.parse_args()
Cutoff_16S=0.97
Cutoff_HGT=0.99
Cutoff_aa=0.8
Cutoff_extended=0.8 # outside is 60% if the gene is 1000 bp
Cutoff_extended2=0.99 # outside is 99% if the gene is 1000 bp
Cutoff_fastani = 0.95
Cutoff_fastani2 = 0.998
Hit_length = 0.9
Cutoff_fastani_hit = 0.2
Cutoff_fastani_hit2 = 0.95
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


# set output
script_i = 0
script_i_max = int(args.th)
if args.s == 'None':
    input_dir = os.path.join(args.r,'summary')
    result_dir = os.path.join(args.r, 'HGT')
else:
    input_dir = args.s
    result_dir = os.path.join(args.s, '../HGT')
try:
    os.mkdir(result_dir)
except OSError:
    pass
try:
    os.mkdir(result_dir + '/sub_fun_summary')
except OSError:
    pass
try:
    os.mkdir(result_dir + '/sub_fun')
except OSError:
    pass
if args.ani != 'None':
    try:
        os.mkdir(os.path.join(result_dir, '../fastANI'))
    except OSError:
        pass
workingdir=os.path.abspath(os.path.dirname(__file__))
try:
    os.system('rm -rf HGT_subscripts')
    os.mkdir('HGT_subscripts')
except OSError:
    pass

################################################### new class #########################################################
__metaclass__ = type


class HGT_function:
    # create a class to store HGT_function
    'a class to store HGT_function'

    def init(self, function,type,cutoff,range16S_same,outputfile,clusterfile):
        self.function = function
        self.type = type
        self.cutoff = cutoff
        self.range16S_same = range16S_same
        self.Diff_16S_min = Species_cutoff
        self.mge_to_genome = 0
        self.mge_to_mge = 0
        self.sameCluster_16S_Set = set()
        self.diffCluster_16S_Set = set()
        self.outputfile = outputfile
        self.outputfile_list = []
        self.outputdiff_list = []
        self.outputsame_list = []
        self.outputmge_list = []
        if clusterfile == "/dev/null":
            self.output_file_diff = clusterfile
            self.output_file_same = clusterfile
            self.output_file_mge = clusterfile
        else:
            self.output_file_diff = clusterfile + '.diff.cluster'
            self.output_file_same = clusterfile + '.same.cluster'
            self.output_file_mge = clusterfile + '.mge.cluster'
        #output_file = open(self.outputfile,'w')
        #output_file.write('function\ttype_cutoff\tgenome_pair\tid_gene\tid_16S\n')
        #output_file.close()
        self.Same_genome_set = set()
        self.Diff_genome_set = set()

    def addsame16Scluster(self, cluster):
        self.sameCluster_16S_Set.add(cluster)

    def adddiff16Scluster(self, cluster):
        self.diffCluster_16S_Set.add(cluster)

    def adddiffgenome_set(self, genome_pair):
        self.Diff_genome_set.add(genome_pair)

    def addsamegenome_set(self, genome_pair):
        self.Same_genome_set.add(genome_pair)

    def setDiff_16S_min(self, newlow):
        self.Diff_16S_min = min(self.Diff_16S_min,newlow)

    def addmge_to_genome(self):
        self.mge_to_genome += 1

    def addmge_to_mge(self):
        self.mge_to_mge += 1

    def addoutput(self, lines):
        self.outputfile_list.append(lines)

    def adddiff(self, lines):
        self.outputdiff_list.append(lines)

    def addsame(self, lines):
        self.outputsame_list.append(lines)

    def addmge(self, lines):
        self.outputmge_list.append(lines)

    def writeoutput(self):
        if self.outputfile_list!=[]:
            output_file = open(self.outputfile, 'a')
            output_file.write(''.join(self.outputfile_list))
            self.outputfile_list = []
            output_file.close()
        if self.outputdiff_list != []:
            output_file = open(self.output_file_diff, 'a')
            output_file.write(''.join(self.outputdiff_list))
            self.outputdiff_list = []
            output_file.close()
        if self.outputsame_list != []:
            output_file = open(self.output_file_same, 'a')
            output_file.write(''.join(self.outputsame_list))
            self.outputsame_list = []
            output_file.close()
        if self.outputmge_list != []:
            output_file = open(self.output_file_mge, 'a')
            output_file.write(''.join(self.outputmge_list))
            self.outputmge_list = []
            output_file.close()


class Tree_list:
    # create a class to store HGT_function
    'a class to store HGT_function'

    def init(self, node1, node2):
        self.parent = node1
        self.node = node2
        self.nodelist = []
        # small
        self.leftchild = ''
        self.leftchildid = 0
        # big
        self.rightchild = ''
        self.rightchildid = 0

    def addnode(self, newnode, ID, hit_length):
        # not finished
        if ID >= Cutoff_fastani:
            if hit_length >= Cutoff_fastani_hit:
                # same species
                self.nodelist.append(newnode)
        elif self.leftchild == '':
            self.leftchild = newnode
            self.leftchildid = ID
        elif ID > self.leftchildid:
            self.rightchild = newnode
            self.rightchildid = ID


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
        self.mutposition = set()
        self.NSratio = [0,0]
        self.protein = ''

    def addmutposition(self, position):
        self.mutposition.add(position)

    def addprotein(self, aa):
        self.protein += aa

    def mutpositioncal(self):
        self.mutpositionsum1 = {0: 0,
                         1: 0,
                         2: 0}
        self.mutpositionsum2 = {0: 0,
                                1: 0,
                                2: 0}
        self.mutposition = list(self.mutposition)
        self.mutposition.sort()
        Total = len(self.mutposition)
        if Total > 1:
            for i in range(0,len(self.mutposition)-1):
                self.mutpositionsum1[abs(self.mutposition[i+1] - self.mutposition[i]) % 3] += 1
                self.mutpositionsum2[(self.mutposition[i]) % 3] += 1
            self.mutpositionsum2[(self.mutposition[-1]) % 3] += 1

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
        refSNP_pair_sum_all = Ref_NSratio.get(self.gene, 'None')
        refSNP_pair_sum = refSNP_pair_sum_all[0]
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


def run_vcf(Group_sample):
    cmds = ''
    for groups in Group_sample:
        sub_samples = Group_sample[groups]
        total_sample = len(sub_samples)
        tempbamoutput = os.path.join(args.r, 'summary/vcf/%s.flt.vcf' % (groups))
        print('%s merge files and call SNPs %s' % (
            datetime.now(), os.path.join(args.r, 'summary/vcf/%s.flt.vcf' % (groups))))
        if total_sample <= 100:
            cmds += (
                    '%s mpileup --threads %s -a FMT/AD -B -Ou -d3000 -f %s %s  | %s call --ploidy 1 -A --threads %s -m > %s' % (
                args.bcf, min(args.th, 40), args.db, ' '.join(sub_samples), args.bcf, min(args.th, 40),
                tempbamoutput.replace('.flt.vcf', '.raw.vcf')))
            cmds += ('\n%s filter --threads %s -s LowQual -e \'DP>100\' %s > %s \n' % (
                args.bcf, min(args.th, 40), tempbamoutput.replace('.flt.vcf', '.raw.vcf'), tempbamoutput))
            cmds += ('\n%s view -v snps %s > %s \n' % (
                args.bcf, tempbamoutput, tempbamoutput.replace('.flt.vcf', '.flt.snp.vcf')))
        else:
            for i in range(1, int(total_sample / 100)):
                cmds += (
                        '%s mpileup --threads %s -a FMT/AD -q30 -B -Ou -d3000 -f %s %s | %s call --ploidy 1 -A --threads %s -m > %s.%s' % (
                    args.bcf, min(args.th, 40), args.db, ' '.join(sub_samples[(i - 1) * 100:(i * 100)]), args.bcf,
                    min(args.th, 40),
                    tempbamoutput.replace('.flt.vcf', '.raw.vcf'), i))
                cmds += ('\n%s filter --threads %s -s LowQual -e \'DP>=100\' %s.%s > %s.%s \n' % (
                    args.bcf, min(args.th, 40), tempbamoutput.replace('.flt.vcf', '.raw.vcf'), i, tempbamoutput, i))
                cmds += ('\n%s view -v snps %s.%s > %s.%s \n' % (
                    args.bcf, tempbamoutput, i, tempbamoutput.replace('.flt.vcf', '.flt.snp.vcf'), i))
            cmds += (
                    '%s mpileup --threads %s -a FMT/AD -q30 -B -Ou -d3000 -f %s %s | %s call --ploidy 1 -A --threads %s -m > %s.%s' % (
                args.bcf, min(args.th, 40), args.db, ' '.join(sub_samples[(i * 100):total_sample]), args.bcf,
                min(args.th, 40),
                tempbamoutput.replace('.flt.vcf', '.raw.vcf'), i + 1))
            cmds += ('\n%s filter --threads %s -s LowQual -e \'DP>100\' %s.%s > %s.%s \n' % (
                args.bcf, min(args.th, 40), tempbamoutput.replace('.flt.vcf', '.raw.vcf'), i + 1, tempbamoutput,
                i + 1))
            cmds += ('\n%s view -v snps %s.%s > %s.%s \n' % (
                args.bcf, tempbamoutput, i + 1, tempbamoutput.replace('.flt.vcf', '.flt.snp.vcf'), i + 1))
    fscripts = open('vcf_group_run.sh', 'w')
    fscripts.write('#!/bin/bash\nsource ~/.bashrc\npy37\n%s' % (''.join(cmds)))
    fscripts.close()
    return cmds


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
    Total = int(len(record_seq)/3)
    temp_SNP_gene = SNP_gene()
    temp_SNP_gene.init(record_name)
    for i in range(0, (Total - position)):
        codon = record_seq[(i * 3 + position):((i + 1) * 3 + position)]
        codon_NSratio = codontable_NSratio[codon]
        temp_SNP_gene.addprotein(codontable[codon])
        for pair in codon_NSratio.SNP_pair:
            temp_SNP_gene.addSNP_pair(pair, 0, codon_NSratio.SNP_pair[pair][0], Total)
            temp_SNP_gene.addSNP_pair(pair, 1, codon_NSratio.SNP_pair[pair][1], Total)
    temp_SNP_gene.sum_SNP_pair()
    return [temp_SNP_gene.SNP_pair_sum,temp_SNP_gene.protein,position]


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

def freq_call_sub(vcf_file):
        SNP = dict()
        SNP_count = dict()
        Total = 0
        for lines in open(vcf_file, 'r'):
            if not lines.startswith("#"):
                lines_set = lines.split('\t')
                Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
                Chr = lines_set[0]
                SNP_count.setdefault(Chr, 0)
                Position = int(lines_set[1])
                Chr_position = '%s--%s' % (Chr, Position)
                Allels_frq = [0, 0, 0, 0]
                SNP.setdefault(Chr_position, Allels_frq)
                Allels_frq = SNP[Chr_position]
                if Total == 0:
                    Total = len(lines_set) - 9
                if "INDEL" not in lines_set[7] and (lines_set[6] != 'LowQual' or Depth >= 100):
                    # Depth >= 10 for genome or Depth >= 100 for metagenomes
                    # new SNPs on Chr
                    SNP_count[Chr] += 1
                    allels_set = [lines_set[3]]
                    if '.' not in lines_set[4]:
                        allels_set += lines_set[4].split(',')
                    for Subdepth_all in lines_set[9:]:
                        Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                        if Subdepth != ['0', '0']:
                            for num_allels in range(0, len(allels_set)):
                                allels = allels_set[num_allels]
                                Subdepth_alleles = Subdepth[num_allels]
                                if allels in Allels:
                                    Allels_frq[Allels[allels]] += int(Subdepth_alleles)
                                else:
                                    pass
                    SNP[Chr_position] = Allels_frq
        return [SNP, SNP_count, Total]


def freq_call(vcf_file,cov_file):
    Output =  []
    Output2 = []
    Output3 = []
    all_SNP_gene_temp = dict()
    try:
        SNP, SNP_count, Total = freq_call_sub(vcf_file)
        SNP_all, SNP_count_all, Total = freq_call_sub(cov_file)
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
                    SNP_gene_temp.addmutposition(position)
                    # calculate N or S
                    refSNP_pair_sum_all = Ref_NSratio.get(Chr, 'None')
                    # calculate Allele frequency
                    Major_ALT, Minor_ALT = ALT_freq(SNP[Chr_position])
                    REF = Major_ALT[0]
                    Ref_frq = Major_ALT[1]
                    if refSNP_pair_sum_all != 'None':
                        #  expected NS ratio calculated
                        refSNP_condon_start = refSNP_pair_sum_all[-1]
                        codon_start = position - 1 - int((position - 1)%3) + refSNP_condon_start
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
                                    temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                                    SNP_pair = transitions(REF, ALT)
                                    SNP_gene_temp.addSNP_pair(SNP_pair, N_S_set[temp_NorS],
                                                              ALT_frq,
                                                              Depth_position_snp)
                                    SNP_gene_temp.addmutposition(position)
                                    lines += '\t%s:%s:%s\t%s:%s:%s\t%s' % (REF, Ref_frq, Ref_seq_aa,
                                                                           ALT, ALT_frq, SNP_seq_aa, temp_NorS)
                    else:
                        for minor in Minor_ALT:
                            ALT = minor[0]
                            ALT_frq = minor[1]
                            lines += '\t%s:%s:None\t%s:%s:None\tNone' % (REF, Ref_frq,
                                                                   ALT, ALT_frq)
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
                    # Output.pop(Depth[Chr][2][position], None)
                    # Output2.pop(Depth[Chr][2][position], None)
                    # Output3.pop(Depth[Chr][2][position], None)
            SNP_gene_temp.depth = statistics.mean(temp_depth_filter)
            SNP_gene_temp.cov = float(SNP_gene_temp.cov) / float(Gene_length)
            SNP_gene_temp.mutpositioncal()
            new_line = '%s\t%s\t%s\t%s\t%s\t%s\t%s:%s:%s\t%s:%s:%s' % (sample_name, Chr, Total,
                                                       '%.3f' % (SNP_gene_temp.depth),
                                                       '%.3f' % (SNP_gene_temp.cov), Gene_length,
                                               SNP_gene_temp.mutpositionsum1[0],SNP_gene_temp.mutpositionsum1[1],
                                               SNP_gene_temp.mutpositionsum1[2],SNP_gene_temp.mutpositionsum2[0],
                                            SNP_gene_temp.mutpositionsum2[1],SNP_gene_temp.mutpositionsum2[2])
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
        #foutput = open(cov_file + '.frq.snp.position.diff', 'w')
        #foutput.write('\n'.join(Position_diff_new))
        #foutput.close()
    except (IOError, FileNotFoundError):
        pass
    if args.strainfinder!= 'None':
        strain_finder(cov_file + '.frq.snp')

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
                                print(Subdepth)
                                if Subdepth != ['0', '0']:
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


def Calculate_length(file_name):
    DB_length=set()
    try:
        for lines in open(file_name + '.length', 'r'):
            DB_length.add(float(str(lines.split('\t')[-1]).replace('\n','')))
    except (IOError,FileNotFoundError):
        Fasta_name = open(file_name, 'r')
        f = open(file_name + '.length', 'w')
        for record in SeqIO.parse(Fasta_name, 'fasta'):
            f.write(str(record.id) + '\t'  + str(
                len(record.seq)) + '\n')
            DB_length.add(len(str(record.seq)))
        f.close()
    DB_length_min = min(DB_length)
    if args.dbf == 1:
        DB_length_out = [DB_length_min, DB_length_min/3.0]
    else:
        DB_length_out = [DB_length_min*3.0, DB_length_min]
    return DB_length_out


def split_string_last(input_string,substring):
    last_loci = input_string.rfind(substring)
    if last_loci > -1:
        return input_string[0 : last_loci]
    else:
        return input_string


def checkfile(filename,i):
    try:
        f1 = open(filename,'r')
        if os.path.getsize(filename) > 0:
            for lines in f1:
                try:
                    lines.split('\t',maxsplit=i+1)[i]
                    return 'not empty'
                except (IOError, FileNotFoundError):
                    return 'wrong content by spliting %s \\t' % (str(i))
                break
        else:
            return 'empty'
    except (IOError,FileNotFoundError):
        return 'non-existed'


def genome_com(genome1, genome2):
    if genome1 == genome2:
        return 'same'
    elif genome1 < genome2:
        return '%s-%s'%(genome1,genome2)
    else:
        return 'skip'


def usearch_16S_load(input_file):
    Checkoutput = checkfile(input_file, 2)
    if Checkoutput == 'not empty':
        for lines in open(input_file,'r'):
            line_set = lines.split('\t',maxsplit=4)
            genome_16S = genome_com(line_set[0], line_set[1])
            if genome_16S not in ['same','skip']:
                ID_16S.setdefault(genome_16S, float(line_set[2])/100.0)
    else:
        print('file %s is %s' % (input_file, Checkoutput))


def self_clustering(input_usearch, output_uc):
    # for sequences without any hits, remove from clustering
    output_uc = open(output_uc, 'w')
    clusters = dict()
    cluster_num = 0
    Checkoutput = checkfile(input_usearch, 2)
    if Checkoutput == 'not empty':
        for lines in open(input_usearch):
                line_set = lines.split('\t',maxsplit=4)
                Gene1 = line_set[0]
                Gene2 = line_set[1]
                ID = float(line_set[2])/100.0
                genome_16S = genome_com(Gene1, Gene2)
                if genome_16S not in ['same', 'skip']:
                    ID_16S.setdefault(function_com(Gene1, Gene2), ID)
                if Gene1 not in clusters and Gene2 not in clusters:
                    cluster_num += 1
                    if ID < (Cutoff_16S):
                        # diff cluster
                        clusters.setdefault(Gene1, cluster_num)
                        output_uc.write('*\t%s\t*\t*\t*\t*\t*\t*\t%s\t*\n' % (cluster_num, Gene1))
                        cluster_num += 1
                        clusters.setdefault(Gene2, cluster_num)
                        output_uc.write('*\t%s\t*\t*\t*\t*\t*\t*\t%s\t*\n' % (cluster_num, Gene2))
                    else:
                        # same cluster
                        clusters.setdefault(Gene1, cluster_num)
                        output_uc.write('*\t%s\t*\t*\t*\t*\t*\t*\t%s\t*\n' % (cluster_num, Gene1))
                        clusters.setdefault(Gene2, cluster_num)
                        output_uc.write('*\t%s\t*\t*\t*\t*\t*\t*\t%s\t*\n' % (cluster_num, Gene2))
                else:
                    if ID < (Cutoff_16S):
                        # diff cluster
                        if Gene1 not in clusters:
                            cluster_num += 1
                            clusters.setdefault(Gene1, cluster_num)
                            output_uc.write('*\t%s\t*\t*\t*\t*\t*\t*\t%s\t*\n' % (cluster_num, Gene1))
                        elif Gene2 not in clusters:
                            cluster_num += 1
                            clusters.setdefault(Gene2, cluster_num)
                            output_uc.write('*\t%s\t*\t*\t*\t*\t*\t*\t%s\t*\n' % (cluster_num, Gene2))
                        else:
                            # both genes are set
                            pass
                    else:
                        # same cluster
                        if Gene1 not in clusters:
                            Gene2_cluster = clusters[Gene2]
                            clusters.setdefault(Gene1, Gene2_cluster)
                            output_uc.write('*\t%s\t*\t*\t*\t*\t*\t*\t%s\t*\n' % (Gene2_cluster, Gene1))
                        elif Gene2 not in clusters:
                            Gene1_cluster = clusters[Gene1]
                            clusters.setdefault(Gene2, Gene1_cluster)
                            output_uc.write('*\t%s\t*\t*\t*\t*\t*\t*\t%s\t*\n' % (Gene1_cluster, Gene2))
                        else:
                            # both genes are set
                            pass
    else:
        print('file %s is %s' % (input_usearch,Checkoutput))
    output_uc.close()


def loci_seq(record_name):
    loci_last = record_name.rfind('_')
    return [int(record_name[record_name.rfind('_', 0, loci_last) + 1:loci_last]),
            int(record_name[loci_last + 1:])]


def function_load(input_file,type_fasta):
    Function_Set = dict()
    Checkoutput = checkfile(input_file, 8)
    if Checkoutput == 'not empty':
        for lines in open(input_file, 'r'):
            line_set = lines.split('\t', maxsplit=3)
            gene = line_set[1]
            if args.g == 'F':
                function = line_set[0].replace("(", "").replace(")", "").replace(".", "_").replace(" ", "_").replace(
                    "-", "_")
            else:
                function = line_set[2]  # reference gene name
            if type_fasta == 'dna':
                    loci_new = loci_seq(gene)
                    gene = gene[0:gene.rfind('_', 0, (gene.rfind('_') - 1))]
                    # query gene
                    Function_Set.setdefault(gene, [[], []])
                    loci_set = [int(loci_new[0]), int(loci_new[1])]
                    if loci_set not in Function_Set[gene][-1]:
                        Function_Set[gene][-1].append(loci_set)
                        Function_Set[gene][0].append([function, loci_set])
            else:
                    Function_Set.setdefault(gene, function)
    else:
        print('file %s is %s' % (input_file, Checkoutput))
    return Function_Set


def compare_loci(loci_new,loci_ref):
    if min(loci_new[0],loci_new[1]) <= min(loci_ref[0],loci_ref[1]) and\
        max(loci_new[0], loci_new[1]) >= max(loci_ref[0], loci_ref[1]):
        return True
    else:
        return False


def function_find(Function_Set, Genome, type_fasta):
    if not Genome.startswith("reference"):
        # query genes
        if type_fasta == 'aa':
            return Function_Set.get(Genome)
        else:
            loci_last_1 = Genome.rfind('_')
            loci_last_2 = Genome.rfind('_',0,loci_last_1)
            #loci_last_3 = Genome.rfind('_', 0, loci_last_2)
            #loci_last_4 = Genome.rfind('_', 0, loci_last_3)
            loci_new = [int(Genome[loci_last_2 + 1:loci_last_1]),
                    int(Genome[loci_last_1 + 1:])]
            Genome_name = Genome[0:loci_last_2]
            for functions in Function_Set.get(Genome_name)[0]:
                    if type_fasta == 'dna':
                        if loci_new == functions[-1]:
                            return functions[0]
                    else:
                        loci_ref = functions[-1]
                        if compare_loci(loci_new, loci_ref):
                            return functions[0]
            return "reference"
    else:
        return "reference"


# deduplicate input_fasta
def deduplicate(input_fasta,Function_Set,type_fasta):
    unique_set = dict()
    record_list = dict()
    i = 0
    for record in SeqIO.parse(input_fasta, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        if record_seq not in unique_set:
            function = function_find(Function_Set, record_id, type_fasta)
            new_name = '%s-%s' % (function, i)
            unique_set.setdefault(record_seq, [new_name,len(record_seq)])
            i += 1
            record_list.setdefault(new_name, [record_id])
        else:
            new_name = unique_set[record_seq][0]
            if new_name.startswith('reference'):
                function = function_find(Function_Set, record_id, type_fasta)
                new_new_name = new_name.replace('reference',function)
                unique_set[record_seq][0] = new_new_name
                record_list[new_new_name] = record_list.pop(new_name)
                new_name = new_new_name
            record_list[new_name].append(str(record.id))
    fout1 = open(input_fasta + '.unique','w')
    fout1_list = []
    fout2 = open(input_fasta + '.unique_list','w')
    fout2_list = []
    fout3 = open(input_fasta + '.unique_length', 'w')
    fout3_list = []
    for record_seq in unique_set:
        record_id = unique_set[record_seq][0]
        record_len = unique_set[record_seq][1]
        if not record_id.startswith('reference'):
            OUtput = 1
            if args.db == '/scratch/users/anniz44/scripts/database/AHR.aa.db':
                # select function
                #if not any(subfun in record_id for subfun in ['pksN','pksL','pksM','fldH','porB','fldB','acdA']):
                if not any(
                            subfun in record_id for subfun in ['pksN', 'pksL', 'pksM']):
                    OUtput = 0
            if OUtput == 1:
                fout1_list.append('>%s\n%s\n'%(record_id,record_seq))
                fout3_list.append('%s\t%s\t\n' % (record_id, record_len))
    fout1.write(''.join(fout1_list))
    fout3.write(''.join(fout3_list))
    for newrecord in record_list:
        if not newrecord.startswith('reference'):
            OUtput = 1
            if args.db == '/scratch/users/anniz44/scripts/database/AHR.aa.db':
                # select function
                if not any(subfun in newrecord for subfun in ['pksN', 'pksL', 'pksM']):
                #if not any(subfun in newrecord for subfun in ['pksN', 'pksL', 'pksM', 'fldH', 'porB', 'fldB', 'acdA']):
                    OUtput = 0
            if OUtput == 1:
                for oldrecord in record_list[newrecord]:
                    fout2_list.append('%s\t%s\t\n'%(newrecord,oldrecord))
    fout2.write(''.join(fout2_list))
    fout1.close()
    fout2.close()
    fout3.close()


def function_com(function1, function2):
    if function1 < function2:
        return function1 + '-' + function2
    elif function1 == function2:
        return function1
    else:
        return function2 + '-' + function1

def cluster_fastani_pair(ID, hit_length):
    if ID >= Cutoff_fastani:
        #if hit_length >= Cutoff_fastani_hit2 and ID >= Cutoff_fastani2:
            #return 'same_strain'
        if hit_length >= Cutoff_fastani_hit:
            return 'same_species'
        else:
            return 'mix'
    elif ID <= 0.92:
        return 'diff'
    else:
        return 'mix'

def cluster_fastani(fastani_out):
    Checkoutput = checkfile(fastani_out + '.sum', 1)
    if Checkoutput == 'not empty':
        for lines in open(fastani_out + '.sum', 'r'):
            line_set = lines.split('\t', maxsplit=2)
            ID_16S.setdefault(line_set[0], float(line_set[1]))
    else:
        print('file %s is %s' % (fastani_out + '.sum', Checkoutput))
    # read cluster results
    Clusters = dict()
    Strains = dict()
    Clusters_seqs = dict()
    Clusters_com = dict()
    Max_ID_len = 0
    Min_ID_len = 0
    Checkoutput = checkfile(fastani_out + '.cluster', 1)
    if Checkoutput != 'not empty':
        print('%s start cluster for %s' % (datetime.now(), fastani_out))
        cluster = 1
        strain = 1
        Checkoutput2 = checkfile(fastani_out, 2)
        #os.system('sed -i \"s/%s//g\" %s' % ('_final.scaffolds.fasta', fastani_out))
        #os.system('sed -i \"s/%s//g\" %s' % ('_genomic.fna', fastani_out))
        os.system('sed -i \"s/%s//g\" %s' % (args.fa, fastani_out))
        os.system('sed -i \"s/%s/_/g\" %s' % ('-', fastani_out))
        if Checkoutput2 == 'not empty':
            for lines in open(fastani_out, 'r'):
                #line_set = lines.replace('\r', '').replace('\n', '').replace('_final.scaffolds.fasta', '').replace('_genomic.fna', '').replace('args.fa', '').replace('-', '').split('\t')
                line_set = lines.replace('\r', '').replace('\n', '').split('\t')
                genome_16S = genome_com(os.path.split(line_set[0])[1], os.path.split(line_set[1])[1])
                if genome_16S not in ['same', 'skip']:
                    ID = float(line_set[2]) / 100.0
                    hit_length = float(line_set[3]) / float(line_set[4])
                    ID_16S.setdefault(genome_16S, ID)
                    if cluster % 100000 == 0:
                        print('%s clustering %s clusters for %s' % (datetime.now(), cluster, fastani_out))
                    Genome1, Genome2 = genome_16S.split('-')
                    cluster1 = Clusters_seqs.get(Genome1, 0)
                    cluster2 = Clusters_seqs.get(Genome2, 0)
                    compare_pair = cluster_fastani_pair(ID, hit_length)
                    length1 = len(Genome1)
                    length2 = len(Genome2)
                    Max_ID_len = max(Max_ID_len, length1, length2)
                    if Min_ID_len == 0:
                        Min_ID_len = min(length1, length2)
                    Min_ID_len = min(Min_ID_len, length1, length2)
                    if cluster1 != 0 or cluster2 != 0:
                        #if compare_pair == 'same_strain' and compare_pair == 'same_species':
                        if compare_pair == 'same_species':
                            if cluster1 == 0:
                                Clusters_seqs[Genome1] = cluster2
                            elif cluster2 == 0:
                                Clusters_seqs[Genome2] = cluster1
                            else:
                                cluster_new = min(cluster1, cluster2)
                                Clusters_seqs[Genome1] = cluster_new
                                Clusters_seqs[Genome2] = cluster_new
                        elif compare_pair == 'diff':
                            cluster += 1
                            if cluster1 == 0:
                                Clusters_seqs[Genome1] = cluster
                                cluster1 = cluster
                            elif cluster2 == 0:
                                Clusters_seqs[Genome2] = cluster
                                cluster2 = cluster
                            elif cluster1 == cluster2:
                                print('wrong cluster for %s and %s' % (Genome1,Genome2))
                            cluster_set = genome_com(cluster1, cluster2)
                            if cluster_set not in ['same', 'skip']:
                                Clusters_com.setdefault(cluster_set, [])
                                Clusters_com[cluster_set].append(ID)
                        else:
                            pass
                    else:
                        cluster += 1
                        Clusters_seqs.setdefault(Genome1, cluster)
                        Clusters.setdefault(cluster, [])
                        Clusters[cluster].append(Genome1)
                        if compare_pair!= 'same':
                            cluster += 1
                            cluster_set = genome_com(cluster, cluster-1)
                            if cluster_set not in ['same', 'skip']:
                                Clusters_com.setdefault(cluster_set, [])
                                Clusters_com[cluster_set].append(ID)
                        Clusters_seqs.setdefault(Genome2, cluster)
            print('%s output cluster result for %s' % (datetime.now(), fastani_out + '.cluster'))
            cluster_list = []
            for genome in Clusters_seqs:
                cluster = Clusters_seqs[genome]
                Clusters.setdefault(cluster, [])
                Clusters[cluster].append(genome)
                cluster_list.append('%s\t%s\t\t\n' % (genome, cluster))
            fout = open(fastani_out + '.cluster', 'w')
            fout.write(''.join(cluster_list))
            fout.close()
            print('%s output cluster sum result for %s' % (datetime.now(), fastani_out + '.sum'))
            cluster_list = []
            for cluster_set in Clusters_com:
                cluster_list.append(
                    '%s\t%s\t\t\n' % (cluster_set, '%.5f' % (statistics.mean(Clusters_com[cluster_set]))))
            fout = open(fastani_out + '.sum', 'w')
            fout.write(''.join(cluster_list))
            fout.close()
        else:
            print('file %s is %s\n load 16S results instead' % (fastani_out, Checkoutput2))
            return run_compare(f16s, Function_Set_dna, Cutoff_16S, 0.6, 'dna', 'T')
    else:
        print('%s load cluster for %s' % (datetime.now(), fastani_out + '.cluster'))
        for lines in open(fastani_out + '.cluster', 'r'):
            line_set = lines.split('\t')
            cluster = line_set[1]
            record_name = line_set[0]
            Clusters.setdefault(cluster, [])
            Clusters[cluster].append(record_name)
            length = len(record_name)
            Max_ID_len = max(Max_ID_len, length)
            if Min_ID_len == 0:
                Min_ID_len = length
            Min_ID_len = min(Min_ID_len, length)
            Clusters_seqs.setdefault(record_name, cluster)
    return [Clusters_seqs, Clusters, Max_ID_len, Min_ID_len]


def subprocess_cmd(command):
    process = subprocess.Popen(command,stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()


def run_fastani_sub(genome1,genome2_list):
    number = 0
    command = ''
    result_set = dict()
    f1 = open('Filelist.temp1.txt', 'w')
    f1.write('%s\n' % (genome1))
    f1.close()
    for genome2 in genome2_list:
        result_set.setdefault(genome2, ['mix'])
        f2 = open('Filelist.temp2.%s.txt' %(number), 'w')
        f2.write('%s\n' % (genome2))
        f2.close()
        command += ('%s --rl Filelist.temp1.txt --ql Filelist.temp2.%s.txt -o Filelist.temp.%s.output 2> /dev/null;' %
                  (args.ani, number, number))
        number += 1
    subprocess_cmd(command)
    os.system('rm -rf Filelist.temp.all.output')
    os.system('cat Filelist.temp.*.output > Filelist.temp.all.output')
    os.system('sed -i -e \"s/%s//g\"  -e \"s/%s//g\" Filelist.temp.all.output' %(args.fa,'-'))
    for lines in open('Filelist.temp.all.output','r'):
        line_set = lines.replace('\r', '').replace('\n', '').split('\t')
        ID = float(line_set[2]) / 100.0
        genome2 = line_set[1]
        if ID >= Cutoff_fastani and float(line_set[3]) / float(line_set[4]) >= Cutoff_fastani_hit:
            result_set[genome2]=['same', lines]
        elif ID <= 0.92:
            result_set[genome2] = ['diff',lines]
        else:
            result_set[genome2] = ['mix']
    return result_set


def ani_ref(flist_list, Genome_list, fastani_out):
    print('%s build new references output to %s' % (datetime.now(), Genome_list + '.reference'))
    Total_genome = len(flist_list)
    Set_num = max(min(100,int(Total_genome/100)),int(Total_genome/10))
    if Set_num == 0:
        Set_num = int(Total_genome/2)
    Ref_set = random.sample(flist_list, Set_num)
    Nonuniq = []
    Same_list = []
    fastani_out_list_temp = set()
    Total = len(Ref_set)
    genome_num = 0
    Run_set = []
    for i in range(Total):
        for j in range(i + 1, Total):
            genome1 = Ref_set[i]
            genome2 = Ref_set[j]
            if genome1 not in Nonuniq and genome2 not in Nonuniq:
                genome_num += 1
                if len(Run_set) <= args.th:
                    Run_set.append(genome2)
                else:
                    result_fastani = run_fastani_sub(genome1, Run_set)
                    Run_set = []
                    for genome_run in result_fastani:
                        if result_fastani[genome_run][0] != 'diff':
                            Nonuniq.append(genome_run)
                            if result_fastani[genome_run][0] == 'same':
                                fastani_out_list_temp.add(result_fastani[genome_run][1])
                                Same_list.append(genome_run)
                        else:
                            fastani_out_list_temp.add(result_fastani[genome_run][1])
                if genome_num % 100 == 0:
                    print('%s run fastANI for %s pairs and kicked out %s ref genomes' % (
                        datetime.now(), genome_num, len(Nonuniq)))
                    flist = open(fastani_out, 'a')
                    flist.write(''.join(fastani_out_list_temp))
                    flist.close()
                    fastani_out_list_temp = set()
            else:
                break
    Ref_set = list(set(Ref_set) - set(Nonuniq))
    flist_list = list(set(flist_list) - set(Same_list))
    Total = len(Ref_set)
    while Total < Set_num and Total*2 < Total_genome:
        New_se1 = random.sample(flist_list, len(Nonuniq))
        Nonuniq = []
        Run_set = []
        for genome1 in Ref_set:
            for genome2 in New_se1:
                if genome1 != genome2:
                    if genome1 not in Nonuniq and genome2 not in Nonuniq:
                        genome_num += 1
                        if len(Run_set) <= args.th:
                            Run_set.append(genome2)
                        else:
                            result_fastani = run_fastani_sub(genome1, Run_set)
                            Run_set = []
                            for genome_run in result_fastani:
                                if result_fastani[genome_run][0] != 'diff':
                                    Nonuniq.append(genome_run)
                                    if result_fastani[genome_run][0] == 'same':
                                        fastani_out_list_temp.add(result_fastani[genome_run][1])
                                        Same_list.append(genome_run)
                                else:
                                    fastani_out_list_temp.add(result_fastani[genome_run][1])
                        if genome_num % 100 == 0:
                                print('%s run fastANI for %s pairs and kicked out %s ref genomes' % (
                                    datetime.now(), genome_num, len(Nonuniq)))
                                flist = open(fastani_out, 'a')
                                flist.write(''.join(fastani_out_list_temp))
                                flist.close()
                                fastani_out_list_temp = set()
                    else:
                        break
        Ref_set = list(set(Ref_set) - set(Nonuniq))
        flist_list = list(set(flist_list) - set(Same_list))
        Total = len(Ref_set)
    flist = open(Genome_list + '.reference', 'a')
    flist.write('\n'.join(Ref_set)+'\n')
    flist.close()
    flist = open(fastani_out, 'a')
    flist.write(''.join(fastani_out_list_temp))
    flist.close()
    return Ref_set


def run_fastani_genome(Total,Total_genome,fastani_out,flist_list,Ref_set,last='N'):
    # run fastANI
    print('%s run fastANI for %s ref and %s query and output to %s' % (
    datetime.now(), Total, Total_genome, fastani_out))
    genome_num = 1
    Nonuniq = []
    fastani_out_list_temp = set()
    Run_set = []
    for genome1 in Ref_set:
        for genome2 in flist_list:
            if genome1 != genome2 and genome1 not in Nonuniq and genome2 not in Nonuniq:
                genome_num += 1
                if len(Run_set) <= args.th:
                    Run_set.append(genome2)
                else:
                    result_fastani = run_fastani_sub(genome1, Run_set)
                    Run_set = []
                    for genome_run in result_fastani:
                        if result_fastani[genome_run][0] == 'same':
                            fastani_out_list_temp.add(result_fastani[genome_run][1])
                            Nonuniq.append(genome_run)
                        elif result_fastani[genome_run][0] == 'diff' or last == 'T':
                            try:
                                fastani_out_list_temp.add(result_fastani[genome_run][1])
                            except (IOError, FileNotFoundError):
                                pass
                if genome_num % 100 == 0:
                    print('%s run fastANI for %s pairs' % (datetime.now(), genome_num))
                    flist = open(fastani_out, 'a')
                    flist.write(''.join(fastani_out_list_temp))
                    flist.close()
                    fastani_out_list_temp = set()
            else:
                pass
    flist = open(fastani_out, 'a')
    flist.write(''.join(fastani_out_list_temp))
    flist.close()
    return Nonuniq


def run_fastani(fastani_out,input_dir,input_format):
    Checkoutput = checkfile(fastani_out, 2)
    if Checkoutput == 'non-existed':
        Genome_list = 'Filelist.txt'
        flist_list = []
        try:
            for lines in open(Genome_list, 'r'):
                flist_list.append(lines.split('\n')[0])
        except (IOError, FileNotFoundError):
            print('%s generate genome list %s' % (datetime.now(), Genome_list))
            if input_format != 'None':
                for root in glob.glob(os.path.join(input_dir, '*')):
                    list_fasta1 = glob.glob(os.path.join(root, '*' + input_format))
                    if list_fasta1 != []:
                        for genomefile in list_fasta1:
                            flist_list.append(str(genomefile))
            flist = open(Genome_list, 'w')
            flist.write('\n'.join(flist_list)+'\n')
            flist.close()
        # pick up reference genomes
        print('%s pick up references from %s' % (datetime.now(), Genome_list + '.reference'))
        Ref_set = []
        try:
            # load reference genomes
            flist = open(Genome_list + '.reference', 'r')
            for lines in flist:
                Ref_set.append(lines.split('\n')[0])
            flist.close()
        except (IOError, FileNotFoundError):
            # build reference genomes
            Ref_set = (flist_list, Genome_list, fastani_out)
        # start to run fastANI
        flist_list = list(set(flist_list) - set(Ref_set))
        Total_genome = len(flist_list)
        Total = len(Ref_set)
        while Total_genome > Total*2:
            Nonuniq = run_fastani_genome(Total, Total_genome, fastani_out, flist_list, Ref_set)
            flist_list = list(set(flist_list) - set(Ref_set) - set(Nonuniq))
            Total_genome = len(flist_list)
            Ref_set = (flist_list, Genome_list, fastani_out)
            Total = len(Ref_set)
    # run the last fastANI
    if Ref_set != [] and flist_list!=[]:
        Nonuniq = run_fastani_genome(Total, Total_genome, fastani_out, flist_list, Ref_set,'T')
    # load fastani result
    if ID_16S == dict():
        print('%s cluster fastani output %s' % (datetime.now(), fastani_out))
        result = cluster_fastani(fastani_out)
        return result


def run_compare(input_fasta, Function_Set, cutoff1, cutoff2,type_fasta,clustering='F'):
    # deduplicate input_fasta
    if clustering != 'T':
        try:
            f1 = open("%s.unique_length" % (input_fasta), 'r')
        except (IOError,FileNotFoundError):
            print('%s deduplicate %s' % (datetime.now(), input_fasta))
            deduplicate(input_fasta,Function_Set,type_fasta)
        input_fasta = input_fasta + '.unique'
    filesize = int(os.path.getsize(input_fasta))
    if filesize <= 1E+9 and args.u != 'None':
        # smaller than 1000Mb
        try:
            f1 = open("%s.sorted" % (input_fasta), 'r')
        except (IOError,FileNotFoundError):
            os.system('%s -sortbylength %s -fastaout %s.sorted' % (args.u, input_fasta, input_fasta))
        input_fasta = input_fasta + '.sorted'
    try:
        f1 = open("%s.%s.usearch.txt" % (input_fasta, cutoff2), 'r')
    except (IOError,FileNotFoundError):
        print('%s Running usearch for %s' % (datetime.now(), input_fasta))
        if filesize <= 3E+7 and args.u != 'None':
            # smaller than 30Mb
            try:
                f1 = open("%s.udb" % (input_fasta), 'r')
            except (IOError,FileNotFoundError):
                os.system("%s -makeudb_usearch %s -output %s.udb\n"
                          % (args.u, input_fasta, input_fasta))
            if type_fasta == 'aa':
                os.system(
                    "%s  -usearch_global %s -db %s.udb  -id %s -maxaccepts 0 -maxrejects 0 -blast6out %s.%s.usearch.txt  -threads %s\n"
                    % (args.u, input_fasta, input_fasta, cutoff2, input_fasta, cutoff2, str(args.th)))
            else:
                os.system(
                    "%s  -usearch_global %s -db %s.udb  -strand both -id %s -maxaccepts 0 -maxrejects 0 -blast6out %s.%s.usearch.txt  -threads %s\n"
                    % (args.u, input_fasta, input_fasta, cutoff2, input_fasta, cutoff2, str(args.th)))
        elif type_fasta != 'aa' and args.hs != 'None':
            print('%s Using hs-blastn instead of usearch because the input file is larger than 2GB\n'%(datetime.now()))
            try:
                f1 = open("%s.counts.obinary" % (input_fasta), 'r')
            except (IOError,FileNotFoundError):
                os.system('%s -in %s -input_type fasta -dbtype nucl' %
                          (os.path.join(os.path.split(args.bp)[0], 'makeblastdb'),input_fasta))
                os.system('windowmasker -in %s -infmt blastdb -mk_counts -out %s.counts' %
                          (input_fasta, input_fasta))
                os.system('windowmasker -in %s.counts -sformat obinary -out %s.counts.obinary -convert' %
                          (input_fasta, input_fasta))
                os.system('%s index %s' % (args.hs, input_fasta))
            os.system(
                "%s align -db %s -window_masker_db %s.counts.obinary -query %s -out %s.%s.usearch.txt -outfmt 6 -evalue 1 -perc_identity %s -num_threads %s\n" \
                % (args.hs, input_fasta,
                   input_fasta, input_fasta, input_fasta, cutoff2,
                   cutoff2, str(min(int(args.th),40))))
        elif args.dm != 'None' and type_fasta == 'aa':
            print('%s Using diamond instead of usearch because the input file is larger than 2GB\n'%(datetime.now()))
            try:
                f1 = open("%s.dmnd" % (input_fasta), 'r')
            except (IOError,FileNotFoundError):
                os.system('%sdiamond makedb --in %s -d %s.dmnd' %
                          (split_string_last(args.dm, 'diamond'),input_fasta,input_fasta))
            os.system(
                "%sdiamond blastp  --query  %s  --db  %s.dmnd --out %s.%s.usearch.txt --outfmt 6 --id %s --evalue 1 --max-target-seqs 0 --threads %s\n" \
                % (split_string_last(args.dm, 'diamond'), input_fasta,
                   input_fasta, input_fasta, cutoff2,
                   cutoff2, str(min(int(args.th), 40))))
        else:
            if type_fasta == 'aa':
                print('Input file %s is %sMb, too large for usearch\nPlease provide diamond using --dm' % (
                    input_fasta, filesize/1E+7))
            else:
                print('Input file %s is %sMb, too large for usearch\nPlease provide hs-blastn using --hs' % (
                input_fasta, filesize/1E+7))
    if clustering == 'T':
        try:
            f1 = open("%s.uc" % (input_fasta),'r')
        except (IOError,FileNotFoundError):
            print('%s Running usearch cluster for %s' % (datetime.now(),input_fasta))
            # smaller than 1G
            if int(os.path.getsize(input_fasta)) <= 1E+9 and args.u != 'None':
                os.system('%s -sort length -cluster_fast %s -id %s -centroids %s.cluster.aa -uc %s.uc -threads %s'
                          % (args.u, input_fasta, cutoff1, input_fasta, input_fasta, args.th))
            else:
                print('%s We roughly clustered the 16S by 97% identity\n'%(datetime.now()))
                self_clustering('%s.%s.usearch.txt' % (input_fasta, cutoff2), input_fasta + '.uc')
        # load 16S usearch result
        if ID_16S == dict():
            usearch_16S_load("%s.%s.usearch.txt" % (input_fasta, cutoff2))
        print('%s finish running usearch cluster for %s' % (datetime.now(), input_fasta))
        # read cluster results
        Clusters = dict()
        Clusters_seqs = dict()
        Checkoutput = checkfile(input_fasta + '.uc', -2)
        Max_ID_len = 0
        Min_ID_len = 0
        if Checkoutput == 'not empty':
            for lines in open(input_fasta + '.uc', 'r'):
                line_set = lines.split('\t')
                cluster = line_set[1]
                record_name = line_set[-2].split(' ', maxsplit=2)[0]
                Clusters.setdefault(cluster, [])
                Clusters[cluster].append(record_name)
                length = len(record_name)
                Max_ID_len = max(Max_ID_len, length)
                if Min_ID_len == 0:
                    Min_ID_len = length
                Min_ID_len = min(Min_ID_len, length)
                Clusters_seqs.setdefault(record_name, cluster)
        else:
            print('file %s is %s' % (input_fasta + '.uc', Checkoutput))
        return [Clusters_seqs,Clusters,Max_ID_len,Min_ID_len]


def compare_16S(Genome1, Genome2, cutoff,HGT_function1,HGT_function2):
    if Genome1 == 'mge' or Genome2 == 'mge':
        return ["mge"]
    else:
        Genome_set = genome_com(Genome1, Genome2)
        if Genome_set not in ['same','skip']:
            if Genome_set in ID_16S:
                ID = float(ID_16S[Genome_set])
                cluster1 = cluster_16S[0][Genome1]
                cluster2 = cluster_16S[0][Genome2]
                if ID < cutoff:
                    # calculate diff 16S cluster
                    # count diff genome pairs
                    HGT_function1.adddiffgenome_set(Genome_set)
                    HGT_function1.adddiff16Scluster(cluster1)
                    HGT_function1.adddiff16Scluster(cluster2)
                    HGT_function2.adddiffgenome_set(Genome_set)
                    HGT_function2.adddiff16Scluster(cluster1)
                    HGT_function2.adddiff16Scluster(cluster2)
                    return [True, Genome_set, ID]
                else:
                    # calculate same 16S cluster
                    # count same genome pairs
                    HGT_function1.addsamegenome_set(Genome_set)
                    HGT_function1.addsame16Scluster(cluster1)
                    HGT_function2.addsamegenome_set(Genome_set)
                    HGT_function2.addsame16Scluster(cluster1)
                    return [False, Genome_set,ID]
            else:
                return ['16S missing']
        else:
            return ['16S missing']


def compare_fastani(Genome1, Genome2, cutoff,HGT_function1,HGT_function2):
    # calculate total number of 16S
    if Genome1 == 'mge' or Genome2 == 'mge':
        return ["mge"]
    else:
        Genome_set = genome_com(Genome1, Genome2)
        if Genome_set not in ['same','skip']:
            cluster1 = cluster_16S[0].get(Genome1, 0)
            cluster2 = cluster_16S[0].get(Genome2, 0)
            if cluster1 != 0 and cluster2 != 0:
                if cluster1 == cluster2:
                    # calculate same 16S cluster
                    # count same genome pairs
                    HGT_function1.addsamegenome_set(Genome_set)
                    HGT_function1.addsame16Scluster(cluster1)
                    HGT_function2.addsamegenome_set(Genome_set)
                    HGT_function2.addsame16Scluster(cluster1)
                    return [False, Genome_set, cutoff]
                else:
                    cluster_set = genome_com(cluster1, cluster2)
                    ID = float(ID_16S.get(cluster_set,0.8))
                    # calculate diff 16S cluster
                    # count diff genome pairs
                    HGT_function1.adddiffgenome_set(Genome_set)
                    HGT_function1.adddiff16Scluster(cluster1)
                    HGT_function1.adddiff16Scluster(cluster2)
                    HGT_function2.adddiffgenome_set(Genome_set)
                    HGT_function2.adddiff16Scluster(cluster1)
                    HGT_function2.adddiff16Scluster(cluster2)
                    return [True, Genome_set, ID]
            else:
                return ['ANI missing']
        else:
            return ['ANI missing']


def function_pair(Genome1,Genome2):
    function1 = split_string_last(Genome1,'-')
    function2 = split_string_last(Genome2,'-')
    return function_com(function1, function2)


def extract_dna(dna_file,gene_list,output_fasta,type_fasta,script_i):
    output_file = open(output_fasta,'a')
    for record in SeqIO.parse(open(dna_file, 'r'), 'fasta'):
        record_id = str(record.id)
        if record_id in gene_list:
            output_file.write('>%s\n%s\n' %(record_id,str(record.seq)))
    output_file.close()
    if args.mf != 'None':
        try:
            f1 = open("%s.align.nwk" % (output_fasta), 'r')
        except (IOError,FileNotFoundError):
            f1 = open("%s.align.nwk" % (output_fasta), 'w')
            output_script_file = open(('HGT_subscripts/HGTalign.%s.sh')%(int(script_i % script_i_max)), 'a')
            script_i += 1
            output_script_file.write("#!/bin/bash\n")
            output_script_file.write('python %s/remove.duplicated.seq.py -i %s \n' % (workingdir,output_fasta))
            if 'dna' in type_fasta:
                output_script_file.write(
                    "%s --nuc --adjustdirection --quiet --nofft --maxiterate 0 --retree 1 --thread %s %s.dereplicated.id.fasta > %s.align\n"
                    % (args.mf, str(args.th), output_fasta, output_fasta))
                if  args.ft != 'None':
                    output_script_file.write("%s -nt -quiet -fastest -nosupport %s.align > %s.align.nwk\n"
                                      % (args.ft, output_fasta, output_fasta))
            else:
                output_script_file.write(
                    "%s --amino --quiet --retree 1 --maxiterate 0 --nofft --thread %s %s.dereplicated.id.fasta > %s.align\n"
                    % (args.mf, str(args.th), output_fasta, output_fasta))
                if args.ft != 'None':
                    output_script_file.write("%s -quiet -fastest -nosupport %s.align > %s.align.nwk\n"
                                      % (args.ft, output_fasta, output_fasta))
            output_script_file.close()
    return script_i


def add_gene_and_function(Diff_gene_set,Function,Gene):
    Diff_gene_set.setdefault(Function, set())
    Diff_gene_set[Function].add(Gene)


def find_genome(Genome1):
    if Genome1 in mapping:
        return mapping[Genome1]
    elif Genome1.startswith('mge'):
        mapping.setdefault(Genome1, 'mge')
        return 'mge'
    else:
        for i in range(cluster_16S[-1],len(Genome1)):
            if Genome1[i] in ['.', '_']:
                Candidate = Genome1[0:i]
                if Candidate in cluster_16S[0]:
                    mapping.setdefault(Genome1, Candidate)
                    return Candidate
                if i > cluster_16S[-2]:
                    mapping.setdefault(Genome1, 'None')
                    return 'None'
        if Genome1 in cluster_16S[0]:
            mapping.setdefault(Genome1, Genome1)
            return Genome1
        mapping.setdefault(Genome1, 'None')
        return 'None'


def unique_list_load(input_fasta):
    unique_list = dict()
    unique_length = dict()
    for lines in open(input_fasta + '.unique_list','r'):
        loci_last = lines.rfind('\t')
        loci_last_2 = lines.rfind('\t', 0, (loci_last - 1))
        new_name = lines[0:loci_last_2]
        new_genome = lines[loci_last_2+1:loci_last]
        new_genome2 = find_genome(new_genome)
        if new_genome2 != 'None':
            unique_list.setdefault(new_name,[])
            unique_list[new_name].append([new_genome,new_genome2])
    for lines in open(input_fasta + '.unique_length','r'):
        loci_last = lines.rfind('\t')
        loci_last_2 = lines.rfind('\t', 0, (loci_last - 1))
        new_name = lines[0:loci_last_2]
        gene_length = int(lines[loci_last_2+1:loci_last])
        unique_length.setdefault(new_name, gene_length)
    return [unique_list,unique_length]


def com_extended(HGT_dna_extended,HGT_dna):
    HGT_dna_extended.diffCluster_16S_Set = list(
        set(HGT_dna_extended.diffCluster_16S_Set).intersection(HGT_dna.diffCluster_16S_Set))
    HGT_dna_extended.Diff_16S_min = max(HGT_dna_extended.Diff_16S_min, HGT_dna.Diff_16S_min)
    HGT_dna_extended.Diff_genome_set = list(
        set(HGT_dna_extended.Diff_genome_set).intersection(HGT_dna.Diff_genome_set))
    HGT_dna_extended.sameCluster_16S_Set = list(
        set(HGT_dna_extended.sameCluster_16S_Set).intersection(HGT_dna.sameCluster_16S_Set))
    HGT_dna_extended.Same_genome_set = list(
        set(HGT_dna_extended.Same_genome_set).intersection(HGT_dna.Same_genome_set))


def HGT_finder_sum(type_fasta,input_folder,input_prefix,cutoff,cutoff_hit_length,
                   script_i,output_file1,input_fasta,DB_length_min):
    # Setup function list
    Function_list = dict()
    Sequence_list = []
    # load record id mapping
    unique_list_all = unique_list_load(input_fasta)
    unique_list = unique_list_all[0]
    unique_length = unique_list_all[1]
    all_usearch = glob.glob(os.path.join(input_folder, input_prefix))
    line_num = 1
    for files in all_usearch:
        Checkoutput = checkfile(files, 2)
        if Checkoutput == 'not empty':
            Diff_gene_set = dict()
            for lines in open(files,'r'):
                try:
                    line_set = lines.split('\t',maxsplit=4)
                    ID = float(line_set[2]) / 100.0
                    if ID >= cutoff:
                        newGene1 = line_set[0]
                        newGene2 = line_set[1]
                        min_gene_length = float(min(unique_length[newGene1], unique_length[newGene2]))
                        if min_gene_length >= DB_length_min:
                            hit_length = float(line_set[3]) / min_gene_length
                            if hit_length >= cutoff_hit_length:
                                # function count
                                Function = function_pair(newGene1, newGene2)
                                # setup HGT finder class for this function
                                if Function not in Function_list:
                                    HGT_function_temp = HGT_function()
                                    HGT_function_temp.init(Function, type_fasta, cutoff, '%.3f-1.000' % ((Species_cutoff)),
                                                           os.path.join(result_dir,'sub_fun_summary/%s.%s.%.2f.identity.summary.txt'
                                                                        % (Function,type_fasta,cutoff)),
                                                           os.path.join(result_dir + '/sub_fun',
                                                                        "%s.%s.%s.function" %
                                                                        (Function, args.t, type_fasta))
                                                           )
                                    Function_list.setdefault(Function,HGT_function_temp)
                                HGT_function_temp = Function_list[Function]
                                # unique sequence count
                                HGT_function_temp_seq = HGT_function()
                                HGT_function_temp_seq.init('None', type_fasta, 1,
                                                           '%.3f-1.000' % ((Species_cutoff)),
                                                           '/dev/null','/dev/null')
                                if newGene1 == newGene2:
                                    # setup HGT finder class for this sequence
                                    if newGene1 not in Function_list:
                                        HGT_function_temp_seq = HGT_function()
                                        HGT_function_temp_seq.init(newGene1, type_fasta, 1,
                                                               '%.3f-1.000' % ((Species_cutoff)),
                                                                   '/dev/null',
                                                                   os.path.join(result_dir + '/sub_fun',
                                                                                "%s.%s.%s.seq" %
                                                                                (newGene1, args.t, type_fasta))
                                                                   )
                                        Function_list.setdefault(newGene1, HGT_function_temp_seq)
                                        Sequence_list.append(newGene1)
                                    HGT_function_temp_seq = Function_list[newGene1]
                                # process all same genes in newGene1 and newGene2
                                for Gene1_set in unique_list.get(newGene1,[]):
                                    for Gene2_Set in unique_list.get(newGene2,[]):
                                        Gene1 = Gene1_set[0]
                                        Gene2 = Gene2_Set[0]
                                        # find Genome name
                                        Genome1 = Gene1_set[1]
                                        Genome2 = Gene2_Set[1]
                                        Gene_pair = genome_com(Gene1, Gene2)
                                        # cluster = int(os.path.split(files)[-1].split('.fasta.sorted')[0].split('.')[-1])
                                        # not the same gene
                                        if Gene_pair not in ['same', 'skip']:
                                            # record line number and output files
                                            if line_num % 100000 == 0:
                                                for Function in Function_list:
                                                    HGT_function_temp_output = Function_list[Function]
                                                    HGT_function_temp_output.writeoutput()
                                                print('%s HGT_finder processing %s lines' % (datetime.now(), line_num))
                                            # compare gene ID to 16S ID
                                            if args.ani != 'None':
                                                compare_result = compare_fastani(Genome1, Genome2, Cutoff_fastani,HGT_function_temp,HGT_function_temp_seq)
                                            else:
                                                compare_result = compare_16S(Genome1, Genome2, Cutoff_16S,HGT_function_temp,HGT_function_temp_seq)
                                            if compare_result[0] != '16S missing':
                                                line_num += 1
                                                if compare_result[0] != 'mge':
                                                    try:
                                                        if Genome1 != Genome2:
                                                            # count genome pairs
                                                            Genome_pair = compare_result[1]
                                                            if compare_result[0]:
                                                                # different 16S clusters
                                                                # output blast results into diff.cluster
                                                                HGT_function_temp.adddiff(
                                                                    '%s\t%s\t%s\t%s' % (Function, Gene1, Gene2, lines))
                                                                HGT_function_temp_seq.adddiff(
                                                                    '%s\t%s\t%s\t%s' % (Function, Gene1, Gene2, lines))
                                                                if args.mf != 'None':
                                                                    add_gene_and_function(Diff_gene_set, Function, newGene1)
                                                                    add_gene_and_function(Diff_gene_set, Function, newGene2)
                                                                # calculate lowest 16S similarity for same gene in diff 16S clusters
                                                                ID = float(line_set[2]) / 100.0
                                                                lowest_id = compare_result[2]
                                                                HGT_function_temp.setDiff_16S_min(lowest_id)
                                                                HGT_function_temp.addoutput(
                                                                    '%s\t%s_%s\t%s\t%.3f\t%.3f\n'
                                                                    % (Function, type_fasta, cutoff,
                                                                       Genome_pair, ID, lowest_id))
                                                                HGT_function_temp_seq.setDiff_16S_min(lowest_id)
                                                                HGT_function_temp_seq.addoutput(
                                                                    '%s\t%s_%s\t%s\t%.3f\t%.3f\n'
                                                                    % (newGene1, type_fasta, 1,
                                                                       Genome_pair, ID, lowest_id))
                                                            else:
                                                                # same 16S cluster
                                                                # output blast results into same.cluster
                                                                HGT_function_temp.addsame(
                                                                    '%s\t%s\t%s\t%s' % (Function, Gene1, Gene2, lines))
                                                                HGT_function_temp_seq.addsame(
                                                                    '%s\t%s\t%s\t%s' % (Function, Gene1, Gene2, lines))
                                                    except KeyError:
                                                        # missing 16S
                                                        pass
                                                else:
                                                    # mge clusters
                                                    # output blast results into mge.cluster
                                                    HGT_function_temp.addmge(
                                                        '%s\t%s\t%s\t%s' % (Function, Gene1, Gene2, lines))
                                                    HGT_function_temp_seq.addmge(
                                                        '%s\t%s\t%s\t%s' % (Function, Gene1, Gene2, lines))
                                                    # calculate MGEs
                                                    if Genome1 == "mge" and Genome2 == "mge":
                                                        # mge to mge
                                                        HGT_function_temp.addmge_to_mge()
                                                        HGT_function_temp_seq.addmge_to_mge()
                                                    else:
                                                        # mge to genome
                                                        HGT_function_temp.addmge_to_genome()
                                                        HGT_function_temp_seq.addmge_to_genome()
                except IndexError:
                    print('file %s is %s' % (files, 'wrong content by spliting %s \\t' % ('4')))
                    print(lines)
    print('%s HGT_finder processing %s lines' % (datetime.now(), line_num))
    # summarize all functions
    print('%s output HGT results' % (datetime.now()))
    for Function in Function_list:
        HGT_function_temp = Function_list[Function]
        if DNA_function_list != dict() and type_fasta == 'dna_extended':
            # screen out dna_extended pairs that are 99% as dna pairs
            print('%s curating %s dna_extended results' % (datetime.now(),Function))
            HGT_function_temp_dna = DNA_function_list.get(Function,'None')
            if HGT_function_temp_dna != 'None':
                com_extended(HGT_function_temp,HGT_function_temp_dna)
            else:
                pass
                #print('%s dna extended has HGT output but dna has no HGT output' % (Function))
        # summarize diff
        total_combination = []
        for clusters in HGT_function_temp.diffCluster_16S_Set:
            total_16S = len(cluster_16S[1][clusters])
            total_combination.append(total_16S)
        total_combination_sum = sum(total_combination)
        total_pair_diff = total_combination_sum * (total_combination_sum - 1) / 2.0
        for total_16S in total_combination:
            total_pair_diff -= total_16S * (total_16S - 1) / 2.0
        range16S_diff = '%.3f-%.3f' % ((HGT_function_temp.Diff_16S_min), (Species_cutoff))
        hit_pair_diff = len(HGT_function_temp.Diff_genome_set)
        try:
            # number of hits / total number of combination of random 2 genomes
            percentage_diff_pair = float(hit_pair_diff / total_pair_diff)
        except ZeroDivisionError:
            percentage_diff_pair = 0
        # summarize same
        total_combination = []
        total_pair_same = 0
        for clusters in HGT_function_temp.sameCluster_16S_Set:
            total_16S = len(cluster_16S[1][clusters])
            total_combination.append(total_16S)
        for total_16S in total_combination:
            total_pair_same += total_16S * (total_16S - 1) / 2.0
        hit_pair_same = len(HGT_function_temp.Same_genome_set)
        try:
            # number of hits / total number of combination of random 2 genomes
            percentage_same_pair = float(hit_pair_same / total_pair_same)
        except ZeroDivisionError:
            percentage_same_pair = 0
        # summarize diff to same ratio
        try:
            diff_same_ratio = "%.3f" % (percentage_diff_pair / percentage_same_pair)
        except ZeroDivisionError:
            diff_same_ratio = 0
        # output positive results
        if hit_pair_diff + hit_pair_same + HGT_function_temp.mge_to_genome + HGT_function_temp.mge_to_mge > 0:
            if Function in Sequence_list:
                # unique sequences
                Result = [Function, type_fasta, "%.2f" % (1),
                          range16S_diff, "%.1f" % (hit_pair_diff), total_pair_diff, "%.3f" % percentage_diff_pair,
                          HGT_function_temp.range16S_same,
                          "%.1f" % (hit_pair_same), total_pair_same, "%.3f" % percentage_same_pair, diff_same_ratio,
                          (HGT_function_temp.mge_to_genome), (HGT_function_temp.mge_to_mge)]
            else:
                # general function
                Result = [Function, type_fasta, "%.2f" % (cutoff),
                          range16S_diff, "%.1f" % (hit_pair_diff), total_pair_diff, "%.3f" % percentage_diff_pair, HGT_function_temp.range16S_same,
                          "%.1f" % (hit_pair_same), total_pair_same, "%.3f" % percentage_same_pair, diff_same_ratio,
                          (HGT_function_temp.mge_to_genome), (HGT_function_temp.mge_to_mge)]
            for i in range(0, len(Result)):
                Result[i] = str(Result[i])
            output_file1.write('\t'.join(Result) + '\n')
            # output the last lines
            HGT_function_temp.writeoutput()
    # merge all sub_fun_summary
    sum_output = os.path.join(result_dir, 'all.%s.identity.summary.txt'%(type_fasta))
    fsum = open(sum_output,'a')
    fsum.write('function\ttype_cutoff\tgenome_pair\tid_gene\tid_16S\n')
    fsum.close()
    os.system('cat %s >> %s' % (
        os.path.join(result_dir,'sub_fun_summary/*.%s.*.identity.summary.txt' % (type_fasta)),
        os.path.join(result_dir, 'all.%s.identity.summary.txt'%(type_fasta))
    ))
    # extract sequences for alignment
    if args.mf != 'None' and Diff_gene_set != dict():
        print('%s extract sequences' % (datetime.now()))
        for Function in Diff_gene_set:
            if Diff_gene_set[Function] != []:
                script_i = extract_dna(input_fasta+'.unique', Diff_gene_set[Function],
                                       os.path.join(result_dir + '/sub_fun',
                                                    "%s.%s.%s.diff.cluster.fasta" %
                                                    (Function, args.t, type_fasta)),
                                       type_fasta, script_i)
    if type_fasta == 'dna':
        return [script_i,Function_list]
    else:
        return script_i


################################################### Programme #######################################################
# set input
f16s = os.path.join(args.s, args.t + '.all.16S.fasta')
faa = os.path.join(args.s, args.t + '.all.traits.aa.fasta')
fdna = os.path.join(args.s, args.t + '.all.traits.dna.fasta')
fdna_500 = glob.glob(os.path.join(args.s, args.t + '.all.traits.dna.extra*.fasta'))[0]
ID_16S = dict()

# set cutoff
Species_cutoff = Cutoff_16S
if args.ani != 'None':
    Species_cutoff = Cutoff_fastani

# load traits search file for functions
Function_Set_dna = function_load(os.path.join(args.s, args.t + '.all.traits.dna.txt'),'dna')
Function_Set_aa = function_load(os.path.join(args.s, args.t + '.all.traits.aa.txt'),'aa')

# load gene length
DB_length=Calculate_length(args.db)

# comparing genes
print('%s comparing and clustering 16S' % (datetime.now()))
print('%s load metadata' % (datetime.now()))
# load metadata
if args.taxa != 'None':
    Taxa_meta = load_meta(args.taxa)
    # all bam file call vcf
    tempbamoutput = glob.glob(os.path.join(args.r, 'summary/vcf/*.flt.vcf*'))
    print('%s merge files and call SNPs %s' % (datetime.now(), os.path.join(args.r, 'summary/vcf/')))
    if tempbamoutput == []:
        output_files = glob.glob(os.path.join(args.r, 'bwa/0/*.sorted.bam'))
        # grouping
        Sample_group = dict()
        Group_sample = dict()
        if args.taxa != 'None':
            group_samples(output_files, Taxa_meta)
        if Sample_group == dict():
            Group_sample.setdefault('all', output_files)
        else:
            for files in Sample_group:
                groups = Sample_group[files]
                Group_sample.setdefault(groups, [])
                Group_sample[groups].append(files)
        cmds = run_vcf(Group_sample)
        os.system(cmds)
if args.ani == 'None':
    # use 16S
    cluster_16S = run_compare(f16s,Function_Set_dna, Cutoff_16S,0.6,'dna','T')
else:
    # use fastANI results
    cluster_16S = run_fastani(os.path.join(result_dir, '../fastANI/all.fastani.out'),args.i,args.fa)
print ('the range of length of 16S ID is %s to %s' %(cluster_16S[-1],cluster_16S[-2]))
print('%s comparing and clustering DNA' % (datetime.now()))
run_compare(fdna, Function_Set_dna,Cutoff_HGT,Cutoff_HGT,'dna','F')
print('%s comparing and clustering DNA extended' % (datetime.now()))
run_compare(faa, Function_Set_aa, Cutoff_aa,Cutoff_aa,'aa','F')
print('%s comparing and clustering AA' % (datetime.now()))
run_compare(fdna_500, Function_Set_dna, Cutoff_extended,Cutoff_extended,'dna_extended','F')

# load pre-mapping
print('%s loading pre-mapping file' % (datetime.now()))
mapping = dict()
mapping_file = os.path.join(result_dir ,'mapping.genome.16S.txt')
mapping_file_output = 0
try:
    for lines in open(mapping_file,'r'):
        lines_set = lines.split('\t',maxsplit=3)
        mapping.setdefault(lines_set[0],lines_set[1]) # gene_ID, 16S_ID
except (IOError,FileNotFoundError):
    mapping_file_output = 1

# calculate index for HGT
all_output_file = os.path.join(result_dir, 'HGT.summary.dna.%s.aa.%s.16S.%s.txt'
                                   % (Cutoff_HGT, Cutoff_aa, Species_cutoff))
all_output = open(all_output_file,'w')
all_output.write('function_name\ttype\tcutoff\trange16S_diff\thit_pair_diff\ttotal_pair_diff\t'+
                      'percentage_diff_pair\trange16S_same\thit_pair_same\ttotal_pair_same\tpercentage_same_pair\t'+
                      'diff_same_ratio\tmge_to_genome\tmge_to_mge\n')
all_output.close()
# summarize potential HGT
DNA_function_list = dict()
if glob.glob(os.path.join(result_dir + '/sub_fun', "*.%s.%s.*.cluster" %
                                                               (args.t,'dna'))) == []:
    print('%s summarize potential HGT of %s trait with cutoff of %s' % (datetime.now(), 'dna',  Cutoff_HGT))
    all_output = open(all_output_file, 'a')
    Result = HGT_finder_sum('dna', args.s,
                                  os.path.split(fdna)[-1] + '*.unique*.usearch.txt', Cutoff_HGT, Hit_length,script_i,
                   all_output,fdna,DB_length[0])
    script_i = Result[0]
    DNA_function_list = Result[1]
    Result = []
    all_output.close()
if glob.glob(os.path.join(result_dir + '/sub_fun', "*.%s.%s.*.cluster" %
                                                       (args.t, 'aa'))) == []:
    print('%s summarize potential HGT of %s trait with cutoff of %s' % (datetime.now(), 'aa', Cutoff_aa))
    all_output = open(all_output_file, 'a')
    script_i = HGT_finder_sum('aa',args.s,
                              os.path.split(faa)[-1] + '*.unique*.usearch.txt',Cutoff_aa,Hit_length,script_i,
                   all_output,faa,DB_length[1])
    all_output.close()
if glob.glob(os.path.join(result_dir + '/sub_fun', "*.%s.%s.*.cluster" %
                                                               (args.t,'dna_extended'))) == []:
    print('%s summarize potential HGT of %s trait with cutoff of %s' % (datetime.now(), 'extended dna', Cutoff_extended))
    all_output = open(all_output_file, 'a')
    script_i = HGT_finder_sum('dna_extended',args.s,
                              os.path.split(fdna_500)[-1] + '*.unique*.usearch.txt',Cutoff_extended,Hit_length,script_i,
                   all_output,fdna_500,DB_length[0])
    all_output.close()
    print('%s summarize potential HGT of %s trait with cutoff of %s' % (datetime.now(), 'extended dna', Cutoff_HGT))
    all_output = open(all_output_file, 'a')
    script_i = HGT_finder_sum('dna_extended',args.s,
                              os.path.split(fdna_500)[-1] + '*.unique*.usearch.txt',Cutoff_extended2,Hit_length,script_i,
                   all_output,fdna_500,DB_length[0])
    all_output.close()
# collect all bash files
list_of_files = glob.glob('HGT_subscripts/HGTalign.*.sh')
f1 = open("HGTalign.sh", 'w')
f1.write("#!/bin/bash\nsource ~/.bashrc\n")
for file_name in list_of_files:
    f1.write("jobmit %s HGTalign\n" % (file_name))
f1.close()

# output mapping file
if mapping_file_output == 1:
    fout = open (mapping_file,'w')
    fout_set = []
    for genome in mapping:
        fout_set.append(genome+'\t'+mapping[genome]+'\t\n')
    fout.write(''.join(fout_set))
    fout.close()

# dN dS ratio
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
    codontable_NSratio.setdefault(codon, SNP_gene_temp)
    for position in range(0, 3):
        REF = codon[position]
        for ALT in Allels_order:
            if ALT != REF:
                new_codon = causeSNP(codon, position, ALT)
                temp_NorS = dnORds(translate(codon)[0], translate(new_codon)[0])
                SNP_pair = transitions(REF, ALT)
                SNP_gene_temp.addSNP_pair(SNP_pair, N_S_set[temp_NorS], 1, 1)
# load database seq
Ref_seq = dict()
Ref_NSratio = dict()
for record in SeqIO.parse(args.db, 'fasta'):
    record_id = str(record.id)
    record_seq = str(record.seq)
    Ref_seq.setdefault(record_id, record_seq)
    Ref_NSratio.setdefault(record_id,
                           expectNS(record_id, record_seq))

# print(Ref_NSratio)
# calculate average coverage
try:
    all_output = open((os.path.join(args.r, 'summary/all.bam.cov.sum.snp.NS.ratio')), 'r')
except (IOError, FileNotFoundError):
    all_output = open((os.path.join(args.r, 'summary/all.bam.cov.sum.snp.NS.ratio')), 'w')
    all_output.write(
        'metagenome\tgenome\tnumber_sample\tdepth\tcoverage\tgene_length\tSNP_dis_0:1:2\tSNP_codon_0:1:2\tN:S\tobserved_ratio\texpected_ratio\tdNdS\tA-T\tfreq\tN:S\tNSratio' +
        '\tA-C\tfreq\tN:S\tNSratio\tG-C\tfreq\tN:S\tNSratio\tG-T\tfreq\tN:S\tNSratio\tA-G\tfreq\tN:S\tNSratio\tG-A\tfreq\tN:S\tNSratio\n')
    # calculate allel frequency
    output_files = glob.glob(os.path.join(args.r, 'bwa/0/*.flt.snp.vcf'))
    i = 0
    print('%s summarizing allele frequency and SNPs for %s' % (datetime.now(), os.path.join(args.r, 'bwa/0')))
    for vcf_file in output_files:
        cov_file = vcf_file.replace('.flt.snp.vcf', '.raw.vcf')
        i += 1
        if i % 100 == 0:
            print('%s summarizing allele frequency and SNPs for %s files' % (datetime.now(), i))
        freq_call(vcf_file, cov_file)
    all_output.close()
# calculate average coverage
try:
    all_output = open((os.path.join(args.r, 'summary/all.merged.bam.cov.sum.snp.NS.ratio')), 'r')
except (IOError, FileNotFoundError):
    all_output = open((os.path.join(args.r, 'summary/all.merged.bam.cov.sum.snp.NS.ratio')), 'w')
    all_output.write(
        'metagenome\tgenome\tnumber_sample\tdepth\tcoverage\tgene_length\tSNP_dis_0:1:2\tSNP_codon_0:1:2\tN:S\tobserved_ratio\texpected_ratio\tdNdS\tA-T\tfreq\tN:S\tNSratio' +
        '\tA-C\tfreq\tN:S\tNSratio\tG-C\tfreq\tN:S\tNSratio\tG-T\tfreq\tN:S\tNSratio\tA-G\tfreq\tN:S\tNSratio\tG-A\tfreq\tN:S\tNSratio\n')
    # calculate allel frequency
    output_files = glob.glob(os.path.join(args.r, 'summary/vcf/*.flt.snp.vcf*'))
    i = 0
    print('%s summarizing allele frequency and SNPs for %s' % (datetime.now(), os.path.join(args.r, 'summary/vcf/')))
    for vcf_file in output_files:
        cov_file = vcf_file.replace('.flt.snp.vcf', '.raw.vcf')
        i += 1
        if i % 100 == 0:
            print('%s summarizing allele frequency and SNPs for %s files' % (datetime.now(), i))
        freq_call(vcf_file, cov_file)
    all_output.close()
print('%s start calculate SNP dynamics into %s' % (datetime.now(),
                                                   os.path.join(args.r, 'summary/all.snp.dynamics.txt')))
try:
    foutput = open(os.path.join(args.r, 'summary/all.snp.dynamics.txt'), 'r')
except (IOError, FileNotFoundError):
    foutput = open(os.path.join(args.r, 'summary/all.snp.dynamics.txt'), 'w')
    foutput.write('habitat\treference\tsample_number\ttotal_snps\n')
    output_files = glob.glob(os.path.join(args.r, 'summary/vcf/*.flt.snp.vcf*'))
    SNP_dynamics(output_files, foutput)
    foutput.close()

