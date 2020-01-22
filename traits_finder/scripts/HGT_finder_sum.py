import os
from Bio import SeqIO
import argparse
import glob
from datetime import datetime

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
# optional input setup
parser.add_argument("--r",
                    help="input directory or folder of your previous results by Traits_WGD.py",
                    type=str, default='Result',metavar='Result')
parser.add_argument("--s",
                    help="input directory or folder of your previous results by traits summary",
                    type=str, default='None',metavar='summary')
# optional search parameters
parser.add_argument('--g',
                        help="Optional: gene-level HGT finding; --g T (default: function-level; --g F)",
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

################################################## Definition ########################################################
args = parser.parse_args()
Cutoff_16S=0.97
Cutoff_HGT=0.99
Cutoff_aa=0.8
Cutoff_extended=0.8 # outside is 60% if the gene is 1000 bp
Cutoff_extended2=0.99 # outside is 99% if the gene is 1000 bp
Hit_length = 0.9
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

    def init(self, function,type,cutoff,range16S_same,outputfile):
        self.function = function
        self.type = type
        self.cutoff = cutoff
        self.range16S_same = range16S_same
        self.Diff_16S_min = Cutoff_16S
        self.mge_to_genome = 0
        self.mge_to_mge = 0
        self.sameCluster_16S_Set = set()
        self.diffCluster_16S_Set = set()
        self.outputfile = outputfile
        self.outputfile_list = []
        output_file = open(outputfile,'w')
        output_file.write('function\tgenome_pair\tid_gene\tid_16S\n')
        output_file.close()
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

    def writeoutput(self):
        output_file = open(self.outputfile, 'a')
        output_file.write(''.join(self.outputfile_list))
        self.outputfile_list = []
        output_file.close()


################################################### Function ########################################################
def Calculate_length(file_name):
    DB_length=set()
    try:
        for lines in open(file_name + '.length', 'r'):
            DB_length.add(float(str(lines.split('\t')[-1]).replace('\n','')))
    except (IOError):
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
    return input_string[0 : input_string.rfind(substring)]


def checkfile(filename,i):
    try:
        f1 = open(filename,'r')
        if os.path.getsize(filename) > 0:
            for lines in f1:
                try:
                    lines.split('\t',maxsplit=i+1)[i]
                    return 'not empty'
                except IndexError:
                    return 'wrong content by spliting %s \\t' % (str(i))
                break
        else:
            return 'empty'
    except IOError:
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
            loci_last_2 = Genome.rfind('_',0,(Genome.rfind('_')-1))
            loci_last_3 = Genome.rfind('_', 0, loci_last_2)
            loci_last_4 = Genome.rfind('_', 0, loci_last_3)
            loci_new = [int(Genome[loci_last_4 + 1:loci_last_3]),
                    int(Genome[loci_last_3 + 1:loci_last_2])]
            Genome_name = Genome[0:loci_last_4]
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
            fout1_list.append('>%s\n%s\n'%(record_id,record_seq))
            fout3_list.append('%s\t%s\t\n' % (record_id, record_len))
    fout1.write(''.join(fout1_list))
    fout3.write(''.join(fout3_list))
    for newrecord in record_list:
        if not newrecord.startswith('reference'):
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


# merge HGT_finder here (cutoff_set)
def run_compare(input_fasta, Function_Set, cutoff1, cutoff2,type_fasta,clustering='F'):
    # deduplicate input_fasta
    if clustering != 'T':
        try:
            f1 = open("%s.unique_length" % (input_fasta), 'r')
        except IOError:
            print('%s deduplicate %s' % (datetime.now(), input_fasta))
            deduplicate(input_fasta,Function_Set,type_fasta)
        input_fasta = input_fasta + '.unique'
    filesize = int(os.path.getsize(input_fasta))
    if filesize <= 1E+9 and args.u != 'None':
        # smaller than 1000Mb
        try:
            f1 = open("%s.sorted" % (input_fasta), 'r')
        except IOError:
            os.system('%s -sortbylength %s -fastaout %s.sorted' % (args.u, input_fasta, input_fasta))
        input_fasta = input_fasta + '.sorted'
    try:
        f1 = open("%s.%s.usearch.txt" % (input_fasta, cutoff2), 'r')
    except IOError:
        print('%s Running usearch for %s' % (datetime.now(), input_fasta))
        if filesize <= 3E+7 and args.u != 'None':
            # smaller than 30Mb
            try:
                f1 = open("%s.udb" % (input_fasta), 'r')
            except IOError:
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
        elif type_fasta == 'dna' and args.hs != 'None':
            print('%s Using hs-blastn instead of usearch because the input file is larger than 2GB\n'%(datetime.now()))
            try:
                f1 = open("%s.counts.obinary" % (input_fasta), 'r')
            except IOError:
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
            except IOError:
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
        except IOError:
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
                    record_name = line_set[-2].split(' ',maxsplit=2)[0]
                    Clusters.setdefault(cluster, [])
                    Clusters[cluster].append(record_name)
                    length = len(record_name)
                    Max_ID_len = max(Max_ID_len,length)
                    if Min_ID_len == 0:
                        Min_ID_len = length
                    Min_ID_len = min(Min_ID_len,length)
        else:
            print('file %s is %s' % (input_fasta + '.uc', Checkoutput))
        for cluster in Clusters:
            for record_name in Clusters[cluster]:
                Clusters_seqs.setdefault(record_name, str(cluster))
        return [Clusters_seqs,Clusters,Max_ID_len,Min_ID_len]


def compare_16S(Genome1, Genome2, cutoff):
    if Genome1 == 'mge' or Genome2 == 'mge':
        return ["mge"]
    else:
        Genome_set = genome_com(Genome1, Genome2)
        if Genome_set not in ['same','skip']:
            if Genome_set in ID_16S:
                if ID_16S[Genome_set] < cutoff:
                    return [True, Genome_set]
                else:
                    return [False, Genome_set]
            else:
                return ['16S missing']
        else:
            return ['16S missing']


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
        except IOError:
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


def HGT_finder_sum(type_fasta,input_folder,input_prefix,cutoff,cutoff_hit_length,script_i,output_file1,input_fasta,DB_length_min):
    # Setup function list
    Function_list = dict()
    # load record id mapping
    unique_list_all = unique_list_load(input_fasta)
    unique_list = unique_list_all[0]
    unique_length = unique_list_all[1]
    all_usearch = glob.glob(os.path.join(input_folder, input_prefix))
    line_num = 1
    for files in all_usearch:
        Checkoutput = checkfile(files, 2)
        if Checkoutput == 'not empty':
            Outputfiles = dict()
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
                                Function = function_pair(newGene1, newGene2)
                                # setup HGT finder class for this function
                                if Function not in Function_list:
                                    HGT_function_temp = HGT_function()
                                    HGT_function_temp.init(Function, type_fasta, cutoff, '%.3f-1.000' % ((Cutoff_16S)),
                                                           os.path.join(result_dir,'sub_fun_summary/%s.%s.%.2f.identity.summary.txt'
                                                                        % (Function,type_fasta,cutoff)))
                                    Function_list.setdefault(Function,HGT_function_temp)
                                HGT_function_temp = Function_list[Function]
                                # process all same genes to newGene1 and newGene2
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
                                                for Outputfilename in Outputfiles:
                                                    fout = open(Outputfilename, 'a')
                                                    fout.write(''.join(Outputfiles[Outputfilename]))
                                                    fout.close()
                                                Outputfiles = dict()
                                                HGT_function_temp.writeoutput()
                                                print('%s HGT_finder processing %s lines' % (datetime.now(), line_num))
                                            # compare gene ID to 16S ID
                                            compare_result = compare_16S(Genome1, Genome2, Cutoff_16S)
                                            if compare_result[0] != '16S missing':
                                                line_num += 1
                                                if compare_result[0] != 'mge':
                                                    try:
                                                        if Genome1 != Genome2:
                                                            # count genome pairs
                                                            Genome_pair = compare_result[1]
                                                            cluster1 = cluster_16S[0][Genome1]
                                                            cluster2 = cluster_16S[0][Genome2]
                                                            if compare_result[0]:
                                                                # different 16S clusters
                                                                # output blast results into diff.cluster
                                                                output_file_name = os.path.join(result_dir + '/sub_fun',
                                                                                                "%s.%s.%s.diff.cluster" %
                                                                                                (Function, args.t, type_fasta))
                                                                Outputfiles.setdefault(output_file_name, [])
                                                                Outputfiles[output_file_name].append(
                                                                    '%s\t%s\t%s\t%s' % (Function, Gene1, Gene2, lines))
                                                                if args.mf != 'None':
                                                                    add_gene_and_function(Diff_gene_set, Function, newGene1)
                                                                    add_gene_and_function(Diff_gene_set, Function, newGene2)
                                                                # calculate diff 16S clusters
                                                                ID = float(line_set[2]) / 100.0
                                                                # calculate total number of 16S
                                                                HGT_function_temp.adddiffgenome_set(Genome_pair)
                                                                HGT_function_temp.adddiff16Scluster(cluster1)
                                                                HGT_function_temp.adddiff16Scluster(cluster2)
                                                                # calculate lowest 16S similarity for same gene in diff 16S clusters
                                                                lowest_id = Cutoff_16S
                                                                if Genome_pair in ID_16S:
                                                                    lowest_id = ID_16S[Genome_pair]
                                                                HGT_function_temp.setDiff_16S_min(lowest_id)
                                                                HGT_function_temp.addoutput('%s\t%s\t%.3f\t%.3f\n'
                                                                                            % (Function, Genome_pair, ID, lowest_id))
                                                            else:
                                                                # same 16S cluster
                                                                # output blast results into same.cluster
                                                                output_file_name = os.path.join(result_dir + '/sub_fun',
                                                                                                "%s.%s.%s.same.cluster" %
                                                                                                (Function, args.t, type_fasta))
                                                                Outputfiles.setdefault(output_file_name, [])
                                                                Outputfiles[output_file_name].append('%s\t%s\t%s\t%s' % (Function, Gene1, Gene2, lines))
                                                                # calculate same 16S cluster
                                                                # count same genome pairs
                                                                HGT_function_temp.addsamegenome_set(Genome_pair)
                                                                HGT_function_temp.addsame16Scluster(cluster1)
                                                                HGT_function_temp.addsame16Scluster(cluster2)
                                                    except KeyError:
                                                        # missing 16S
                                                        pass
                                                else:
                                                    # mge clusters
                                                    # output blast results into mge.cluster
                                                    output_file_name = os.path.join(result_dir + '/sub_fun',
                                                                                    "%s.%s.%s.mge.cluster" %
                                                                                    (Function, args.t, type_fasta))
                                                    Outputfiles.setdefault(output_file_name, [])
                                                    Outputfiles[output_file_name].append(
                                                        '%s\t%s\t%s\t%s' % (Function, Gene1, Gene2, lines))
                                                    # calculate MGEs
                                                    if Genome1 == "mge" and Genome2 == "mge":
                                                        # mge to mge
                                                        HGT_function_temp.addmge_to_mge()
                                                    else:
                                                        # mge to genome
                                                        HGT_function_temp.addmge_to_genome()
                except IndexError:
                    print('file %s is %s' % (files, 'wrong content by spliting %s \\t' % ('2')))
                    print(lines)
    # output the last lines
    for Outputfilename in Outputfiles:
        fout = open(Outputfilename, 'a')
        fout.write(''.join(Outputfiles[Outputfilename]))
        fout.close()
    print('%s HGT_finder processing %s lines' % (datetime.now(), line_num))
    # summarize all functions
    print('%s output HGT results' % (datetime.now()))
    for Function in Function_list:
        HGT_function_temp = Function_list[Function]
        # summarize diff
        total_combination = []
        for clusters in HGT_function_temp.diffCluster_16S_Set:
            total_16S = len(cluster_16S[1][clusters])
            total_combination.append(total_16S)
        total_combination_sum = sum(total_combination)
        total_pair_diff = total_combination_sum * (total_combination_sum - 1) / 2.0
        for total_16S in total_combination:
            total_pair_diff -= total_16S * (total_16S - 1) / 2.0
        range16S_diff = '%.3f-%.3f' % ((HGT_function_temp.Diff_16S_min), (Cutoff_16S))
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
        # output
        Result = [Function, type_fasta, "%.2f" % (cutoff),
                  range16S_diff, "%.1f" % (hit_pair_diff), total_pair_diff, "%.3f" % percentage_diff_pair, HGT_function_temp.range16S_same,
                  "%.1f" % (hit_pair_same), total_pair_same, "%.3f" % percentage_same_pair, diff_same_ratio,
                  (HGT_function_temp.mge_to_genome), (HGT_function_temp.mge_to_mge)]
        for i in range(0, len(Result)):
            Result[i] = str(Result[i])
        output_file1.write('\t'.join(Result) + '\n')
        HGT_function_temp.writeoutput()
    # merge all sub_fun_summary
    os.system('cat %s > %s' % (
        os.path.join(result_dir, 'sub_fun_summary/*.identity.summary.txt'),
        os.path.join(result_dir, 'all.identity.summary.txt')
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
    return script_i


################################################### Programme #######################################################
# set input
f16s = os.path.join(args.s, args.t + '.all.16S.fasta')
faa = os.path.join(args.s, args.t + '.all.traits.aa.fasta')
fdna = os.path.join(args.s, args.t + '.all.traits.dna.fasta')
fdna_500 = glob.glob(os.path.join(args.s, args.t + '.all.traits.dna.extra*.fasta'))[0]
ID_16S = dict()

# load traits search file for functions
Function_Set_dna = function_load(os.path.join(args.s, args.t + '.all.traits.dna.txt'),'dna')
Function_Set_aa = function_load(os.path.join(args.s, args.t + '.all.traits.aa.txt'),'aa')

# load gene length
DB_length=Calculate_length(args.db)

# comparing genes
print('%s comparing and clustering 16S' % (datetime.now()))
cluster_16S = run_compare(f16s,Function_Set_dna, Cutoff_16S,0.6,'dna','T')
print ('the range of length of 16S ID is %s to %s' %(cluster_16S[-1],cluster_16S[-2]))
print('%s comparing and clustering DNA' % (datetime.now()))
run_compare(fdna, Function_Set_dna,Cutoff_HGT,Cutoff_HGT,'dna','F')
print('%s comparing and clustering DNA extended' % (datetime.now()))
run_compare(faa, Function_Set_aa, Cutoff_aa,Cutoff_aa,'aa','F')
print('%s comparing and clustering AA' % (datetime.now()))
run_compare(fdna_500, Function_Set_dna, Cutoff_extended,Cutoff_extended,'dna','F')

# load pre-mapping
print('%s loading pre-mapping file' % (datetime.now()))
mapping = dict()
mapping_file = os.path.join(result_dir ,'mapping.genome.16S.txt')
mapping_file_output = 0
try:
    for lines in open(mapping_file,'r'):
        lines_set = lines.split('\t',maxsplit=3)
        mapping.setdefault(lines_set[0],lines_set[1]) # gene_ID, 16S_ID
except IOError:
    mapping_file_output = 1

# calculate index for HGT
all_output_file = os.path.join(result_dir, 'HGT.summary.dna.%s.aa.%s.16S.%s.txt'
                                   % (Cutoff_HGT, Cutoff_aa, Cutoff_16S))
all_output = open(all_output_file,'w')
all_output.write('function_name\ttype\tcutoff\trange16S_diff\thit_pair_diff\ttotal_pair_diff\t'+
                      'percentage_diff_pair\trange16S_same\thit_pair_same\ttotal_pair_same\tpercentage_same_pair\t'+
                      'diff_same_ratio\tmge_to_genome\tmge_to_mge\n')
all_output.close()
# summarize potential HGT
if glob.glob(os.path.join(result_dir + '/sub_fun', "*.%s.%s.*.cluster" %
                                                               (args.t,'dna'))) == []:
    print('%s summarize potential HGT of %s trait with cutoff of %s' % (datetime.now(), 'dna',  Cutoff_HGT))
    all_output = open(all_output_file, 'a')
    script_i = HGT_finder_sum('dna', args.s,
                                  os.path.split(fdna)[-1] + '*.unique*.usearch.txt', Cutoff_HGT, Hit_length,script_i,
                   all_output,fdna,DB_length[0])
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
