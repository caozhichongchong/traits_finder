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
# optional input setup
parser.add_argument("--r",
                    help="input directory or folder of your previous results by Traits_WGD.py",
                    type=str, default='Result',metavar='Result')
parser.add_argument("--s",
                    help="input directory or folder of your previous results by traits summary",
                    type=str, default='None',metavar='summary')
# optional search parameters
parser.add_argument('--th',
                        help="Optional: set the thread number assigned for running XXX (default 1)",
                        metavar="1 or more", action='store', default=1, type=int)
# requirement for software calling
parser.add_argument('--u', '--usearch',
                      help="Necessary: use usearch for 16s and gene clustering," +
                           "if your gene fasta file is larger than 2GB, please also install hs-blastn",
                      metavar="usearch",
                      action='store', default='usearch', type=str)
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
Cutoff_extended=0.8
Cutoff_extended2=0.9
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
    os.mkdir(result_dir + '/sub_fun')
except OSError:
    pass

workingdir=os.path.abspath(os.path.dirname(__file__))
try:
    os.mkdir('HGT_subscripts')
except OSError:
    pass

################################################### Function ########################################################
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
        return '%s_%s'%(genome1,genome2)
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
                    ID_16S.setdefault(Gene1 + '-' + Gene2, ID)
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

# merge HGT_finder here (cutoff_set)
def run_compare(input_fasta,cutoff1,cutoff2,datatype,clustering='F'):
    try:
        f1 = open("%s.sorted" % (input_fasta), 'r')
    except IOError:
        os.system('%s -sortbylength %s -fastaout %s.sorted' % (args.u, input_fasta, input_fasta))
    input_fasta = input_fasta + '.sorted'
    try:
        f1 = open("%s.%s.usearch.txt" % (input_fasta, cutoff2), 'r')
    except IOError:
        print('%s Running usearch for %s' % (datetime.now(), input_fasta))
        if int(os.path.getsize(input_fasta)) <= 1E+8:
            # smaller than 100Mb
            try:
                f1 = open("%s.udb" % (input_fasta), 'r')
            except IOError:
                os.system("%s -makeudb_usearch %s -output %s.udb\n"
                          % (args.u, input_fasta, input_fasta))
            os.system(
                "%s  -usearch_global %s -db %s.udb  -strand both -id %s -maxaccepts 0 -maxrejects 0 -blast6out %s.%s.usearch.txt  -threads %s\n"
                % (args.u, input_fasta, input_fasta, cutoff2, input_fasta, cutoff2, str(args.th)))
        elif datatype =='dna':
            print('%s Using hs-blastn instead of usearch because the input file is larger than 2GB\n'%(datetime.now()))
            try:
                f1 = open("%s.counts.obinary" % (input_fasta), 'r')
            except IOError:
                os.system('makeblastdb -in %s -input_type fasta -dbtype nucl' %
                          (input_fasta))
                os.system('windowmasker -in %s -infmt blastdb -mk_counts -out %s.counts' %
                          (input_fasta, input_fasta))
                os.system('windowmasker -in %s.counts -sformat obinary -out %s.counts.obinary -convert' %
                          (input_fasta, input_fasta))
                os.system('%s index %s' % ('hs-blastn', input_fasta))
            os.system(
                "%s align -db %s -window_masker_db %s.counts.obinary -query %s -out %s.%s.usearch.txt -outfmt 6 -evalue 1 -perc_identity %s -num_threads %s\n" \
                % ('hs-blastn', input_fasta,
                   input_fasta, input_fasta, input_fasta, cutoff2,
                   cutoff2, str(min(int(args.th),40))))
        else:
            print('%s Using diamond instead of usearch because the input file is larger than 2GB\n'%(datetime.now()))
            try:
                f1 = open("%s.dmnd" % (input_fasta), 'r')
            except IOError:
                os.system('diamond makedb --in %s -d %s.dmnd' %
                          (input_fasta,input_fasta))
            os.system(
                "%s  --query  %s  --db  %s.dmnd --out %s.%s.usearch.txt --outfmt 6 --id %s --evalue 1 --max-target-seqs 0 --threads %s\n" \
                % ('diamond blastp', input_fasta,
                   input_fasta, input_fasta, cutoff2,
                   cutoff2, str(min(int(args.th), 40))))
    if clustering == 'T':
        try:
            f1 = open("%s.uc" % (input_fasta),'r')
        except IOError:
            print('%s Running usearch cluster for %s' % (datetime.now(),input_fasta))
            # smaller than 2G
            if int(os.path.getsize(input_fasta)) <= 2 * 1E+9:
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


def loci_seq(record_name):
    loci_last = record_name.rfind('_')
    return [int(record_name[record_name.rfind('_', 0, loci_last) + 1:loci_last]),
            int(record_name[loci_last + 1:])]


def function_load(input_file,type_fasta):
    Function_Set = dict()
    Checkoutput = checkfile(input_file, 8)
    if Checkoutput == 'not empty':
        if type_fasta == 'dna':
            for lines in open(input_file,'r'):
                    line_set = lines.split('\t',maxsplit=3)
                    function = line_set[0].replace("(","").replace(")","").replace(".","_").replace(" ","_")
                    gene = line_set[1]
                    loci_new = loci_seq(gene)
                    gene = gene[0:gene.rfind('_', 0, (gene.rfind('_') - 1))]
                    # query gene
                    Function_Set.setdefault(gene,[[],[]])
                    loci_set = [int(loci_new[0]),int(loci_new[1])]
                    if loci_set not in Function_Set[gene][-1]:
                        Function_Set[gene][-1].append(loci_set)
                        Function_Set[gene][0].append([function,loci_set])
        else:
            for lines in open(input_file, 'r'):
                line_set = lines.split('\t')
                function = line_set[0].replace("(", "").replace(")", "").replace(".", "_").replace(" ", "_")
                gene = line_set[1]
                # query gene
                Function_Set.setdefault(gene, function)
    else:
        print('file %s is %s' % (input_file, Checkoutput))
    return Function_Set


def find_genome(Genome1):
    if Genome1 in mapping:
        return mapping[Genome1]
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


def compare_16S(Genome1,Genome2,cutoff):
    if Genome1.startswith('mge') or Genome2.startswith('mge'):
        return "mge"
    else:
        Genome1 = find_genome(Genome1)
        Genome2 = find_genome(Genome2)
        if Genome1 != 'None' and Genome2 != 'None':
            if Genome1 < Genome2:
                Genome_set = Genome1 + '-' + Genome2
            else:
                Genome_set = Genome2 + '-' + Genome1
            if Genome_set in ID_16S:
                if ID_16S[Genome_set] < cutoff:
                    return True
                else:
                    return False
            else:
                return '16S missing'
        else:
            return '16S missing'


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


def function_com(function1, function2):
    if function1 == function2:
        return function1
    else:
        temp = [function1,function2]
        temp.sort()
        return '%s_%s'%(temp[0],temp[1])


def function_pair(Function_Set,Genome1,Genome2,type_fasta):
    function1 = function_find(Function_Set, Genome1, type_fasta)
    function2 = function_find(Function_Set, Genome2, type_fasta)
    if function2 == 'reference':
        function2 = function1
    if function1 == 'reference':
        function1 = function2
    return function_com(function1, function2)


def extract_dna(dna_file,gene_list,input_fasta,type_fasta,script_i):
    output_file = open(input_fasta,'a')
    for record in SeqIO.parse(open(dna_file, 'r'), 'fasta'):
        if str(record.id) in gene_list:
            output_file.write('>%s\n%s\n' %(str(record.id),str(record.seq)))
    output_file.close()
    if args.mf != 'None':
        try:
            f1 = open("%s.align.nwk" % (input_fasta), 'r')
        except IOError:
            f1 = open("%s.align.nwk" % (input_fasta), 'w')
            output_script_file = open(('HGT_subscripts/HGTalign.%s.sh')%(int(script_i % script_i_max)), 'a')
            script_i += 1
            output_script_file.write("#!/bin/bash\n")
            output_script_file.write('python %s/remove.duplicated.seq.py -i %s \n' % (workingdir,input_fasta))
            if 'dna' in type_fasta:
                output_script_file.write(
                    "%s --nuc --adjustdirection --quiet --nofft --maxiterate 0 --retree 1 --thread %s %s.dereplicated.id.fasta > %s.align\n"
                    % (args.mf, str(args.th), input_fasta, input_fasta))
                if  args.ft != 'None':
                    output_script_file.write("%s -nt -quiet -fastest -nosupport %s.align > %s.align.nwk\n"
                                      % (args.ft, input_fasta, input_fasta))
            else:
                output_script_file.write(
                    "%s --amino --quiet --retree 1 --maxiterate 0 --nofft --thread %s %s.dereplicated.id.fasta > %s.align\n"
                    % (args.mf, str(args.th), input_fasta, input_fasta))
                if args.ft != 'None':
                    output_script_file.write("%s -quiet -fastest -nosupport %s.align > %s.align.nwk\n"
                                      % (args.ft, input_fasta, input_fasta))
            output_script_file.close()
    return script_i


def add_gene_and_function(Diff_gene_set,Function,Gene):
    Diff_gene_set.setdefault(Function, [])
    Diff_gene_set[Function].append(Gene)


def compare_traits_16S(Function_Set,type_fasta,input_folder,input_prefix,cutoff,script_i):
    all_usearch = glob.glob(os.path.join(input_folder, input_prefix))
    line_num = 0
    for files in all_usearch:
        Diff_gene_set = dict()
        Checkoutput = checkfile(files, 2)
        if Checkoutput == 'not empty':
            Outputfiles = dict()
            for lines in open(files,'r'):
                try:
                    line_set = lines.split('\t',maxsplit=4)
                    if float(line_set[2])/100.0 >= cutoff:
                        Genome1 = line_set[0]
                        Genome2 = line_set[1]
                        Function = function_pair(Function_Set, Genome1, Genome2, type_fasta)
                        Genome_pair = genome_com(Genome1, Genome2)
                        #cluster = int(os.path.split(files)[-1].split('.fasta.sorted')[0].split('.')[-1])
                        # not the same gene
                        if Genome_pair not in ['same','skip'] and "reference" not in Genome1 and "reference" not in Genome2:
                            line_num += 1
                            if line_num % 1000000 == 0:
                                # output files
                                for Outputfilename in Outputfiles:
                                    fout = open(Outputfilename, 'a')
                                    fout.write(''.join(Outputfiles[Outputfilename]))
                                    fout.close()
                                Outputfiles = dict()
                                print('%s compare_traits_16S processing %s lines' % (datetime.now(), line_num))
                            compare_result = compare_16S(Genome1, Genome2, Cutoff_16S)
                            if compare_result != '16S missing':
                                if compare_result != 'mge':
                                    if compare_result:
                                        # different 16S clusters
                                        output_file_name = os.path.join(result_dir + '/sub_fun',
                                                               "%s.%s.%s.diff.cluster" %
                                                               (Function,args.t,type_fasta))
                                        Outputfiles.setdefault(output_file_name,[])
                                        Outputfiles[output_file_name].append(Function+'\t'+lines)
                                        # record diff gene set
                                        Gene1 = line_set[0]
                                        Gene2 = line_set[1]
                                        if Gene1 not in Diff_gene_set:
                                            add_gene_and_function(Diff_gene_set, Function, Gene1)
                                        if Gene2 not in Diff_gene_set:
                                            add_gene_and_function(Diff_gene_set, Function, Gene2)
                                    else:
                                        pass
                                        # same 16S cluster
                                        output_file_name = os.path.join(result_dir + '/sub_fun',
                                                               "%s.%s.%s.same.cluster" %
                                                               (Function,args.t,type_fasta))
                                        Outputfiles.setdefault(output_file_name, [])
                                        Outputfiles[output_file_name].append(Function + '\t'  + lines)
                                else:
                                    pass
                                    # mge clusters
                                    output_file_name = os.path.join(result_dir + '/sub_fun',
                                                           "%s.%s.%s.mge.cluster" %
                                                           (Function,args.t,type_fasta))
                                    Outputfiles.setdefault(output_file_name, [])
                                    Outputfiles[output_file_name].append(Function + '\t' + lines)
                        if "reference" in Genome1:
                            add_gene_and_function(Diff_gene_set, Function, Genome1)
                        if "reference" in Genome2:
                            add_gene_and_function(Diff_gene_set, Function, Genome2)
                except IndexError:
                    print('file %s is %s' % (files, 'wrong content by spliting %s \\t' % ('2')))
                    print(lines)
            # output remaining lines
            for Outputfilename in Outputfiles:
                fout = open(Outputfilename, 'a')
                fout.write(''.join(Outputfiles[Outputfilename]))
                fout.close()
            Outputfiles = dict()
        else:
            print('file %s is %s' % (files, Checkoutput))
        # extract sequences for alignment
        if args.mf != 'None' and Diff_gene_set !=dict():
            for Function in Diff_gene_set:
                if Diff_gene_set[Function]!=[]:
                    script_i = extract_dna(files.split('.fasta.sorted')[0]+'.fasta.sorted', Diff_gene_set[Function],
                                os.path.join(result_dir + '/sub_fun',
                                             "%s.%s.%s.diff.cluster.fasta" %
                                             (Function, args.t, type_fasta)),
                                type_fasta,script_i)
    return script_i


def HGT_finder_sum(function_name,type_fasta,cutoff,diff,same,mge,output_file1):
    # for each function
    #Diff_cluster=[]
    Diff_16S_min = Cutoff_16S
    #Same_cluster=[]
    Same_genome_set = []
    Diff_genome_set = []
    line_num = 0
    # set up output
    range16S_diff = ''
    hit_pair_diff = 0
    total_pair_diff = 0
    percentage_diff_pair = 0
    range16S_same = '%.3f-1.000' % ((Cutoff_16S))
    hit_pair_same = 0
    total_pair_same = 0
    percentage_same_pair = 0
    diff_same_ratio = 0
    mge_to_genome = 0
    mge_to_mge = 0
    # calculate diff 16S clusters
    Checkoutput = checkfile(diff, 4)
    if Checkoutput == 'not empty':
        output_file2 = open(os.path.join(result_dir,'%s.%s.%.2f.identity.summary.txt'
                                   % (function_name,type_fasta,cutoff)),'a')
        output_file2.write('genome_pair\tid_gene\tid_16S\n')
        Cluster_16S_Set = set()
        output_file2_lines = []
        for lines in open(diff,'r'):
            line_set = lines.split('\t',maxsplit=4)
            Gene1 = line_set[1]
            Gene2 = line_set[2]
            Gene_pair = genome_com(Gene1, Gene2)
            if Gene_pair not in ['same','skip'] and "reference" not in Gene_pair:
                try:
                    line_num += 1
                    if line_num % 100000 == 0:
                        print('%s HGT_finder_sum processing %s lines for file %s' % (datetime.now(), line_num, diff))
                    ID = float(line_set[3])/100.0
                    Genome1 = find_genome(Gene1, cluster_16S[0])
                    cluster1 = cluster_16S[0][Genome1]
                    Genome2 = find_genome(Gene2, cluster_16S[0])
                    cluster2 = cluster_16S[0][Genome2]
                    if Genome1 != Genome2:
                        # count genome pairs
                        Genome_pair = function_com(Genome1, Genome2)
                        if Genome_pair not in Diff_genome_set:
                            Diff_genome_set.append(Genome_pair)
                            # calculate total number of 16S
                            Cluster_16S_Set.add(cluster1)
                            Cluster_16S_Set.add(cluster2)
                        # calculate total number of gene clusters
                        #if line_set[1] not in Diff_cluster:
                        #    Diff_cluster.append(line_set[1])
                        # calculate lowest 16S similarity for same gene in diff 16S clusters
                        lowest_id = Cutoff_16S
                        if Genome1 < Genome2:
                            Genome_set = Genome1 + '-' + Genome2
                        else:
                            Genome_set = Genome2 + '-' + Genome1
                        if Genome_set in ID_16S:
                            lowest_id = ID_16S[Genome_set]
                        Diff_16S_min = min(float(lowest_id),Diff_16S_min)
                        output_file2_lines.append('%s\t%.3f\t%.3f\n'
                                               % (Genome_pair,(ID),lowest_id))
                except KeyError:
                    # missing 16S
                    pass
        # summarize diff
        total_combination = []
        for clusters in Cluster_16S_Set:
            total_16S = len(cluster_16S[1][clusters])
            total_combination.append(total_16S)
        total_pair_diff = sum(total_combination)*sum(total_combination) / 2.0
        for total_16S1 in total_combination:
            total_pair_diff -= total_16S1 * total_16S1 / 2.0
        range16S_diff = '%.3f-%.3f' % ((Diff_16S_min), (Cutoff_16S))
        hit_pair_diff = len(Diff_genome_set)
        try:
            # number of hits / total number of combination of random 2 genomes
            percentage_diff_pair = "%.3f" % (hit_pair_diff / total_pair_diff)
        except ZeroDivisionError:
            percentage_diff_pair = 0
        output_file2.write(''.join(output_file2_lines))
        output_file2.close()
    else:
        print('file %s is %s' % (diff, Checkoutput))
    # calculate same 16S cluster
    Checkoutput = checkfile(same, 3)
    if Checkoutput == 'not empty':
        Cluster_16S_Set = set()
        for lines in open(same,'r'):
            line_set = lines.split('\t',maxsplit=4)
            Gene1 = line_set[1]
            Gene2 = line_set[2]
            Gene_pair = function_com(Gene1, Gene2)
            if Gene_pair not in ['same','skip']:
                try:
                    line_num += 1
                    if line_num % 100000 == 0:
                        print('%s HGT_finder_sum processing %s lines for file %s' % (datetime.now(), line_num, same))
                    Genome1 = find_genome(Gene1, cluster_16S[0])
                    cluster1 = cluster_16S[0][Genome1]
                    Genome2 = find_genome(Gene2, cluster_16S[0])
                    cluster2 = cluster_16S[0][Genome2]
                    # calculate total number of 16S
                    if Genome1 != Genome2:
                        # count genome pairs
                        Genome_pair = function_com(Genome1, Genome2)
                        if Genome_pair not in Same_genome_set:
                            Same_genome_set.append(Genome_pair)
                            # calculate total number of 16S
                            Cluster_16S_Set.add(cluster1)
                            Cluster_16S_Set.add(cluster2)
                        # calculate total number of gene clusters
                        #if line_set[1] not in Same_cluster:
                        #    Same_cluster.append(line_set[1])
                except KeyError:
                    # missing 16S
                    pass
        # summarize same
        total_combination = []
        for clusters in Cluster_16S_Set:
            total_16S = len(cluster_16S[1][clusters])
            total_combination.append(total_16S)
        for total_16S1 in total_combination:
            total_pair_same += total_16S1 * total_16S1 / 2.0
        hit_pair_same = len(Same_genome_set)
        try:
            # number of hits / total number of combination of random 2 genomes
            percentage_same_pair = "%.3f" % (hit_pair_same / total_pair_same)
        except ZeroDivisionError:
            percentage_same_pair = 0
    else:
        print('file %s is %s' % (same, Checkoutput))
    try:
        # number of hits / total number of combination of random 2 genomes
        diff_same_ratio = "%.3f" % (percentage_same_pair / percentage_diff_pair)
    except ZeroDivisionError:
        diff_same_ratio = 0
    # calculate MGEs
    Checkoutput = checkfile(mge, 3)
    if Checkoutput == 'not empty':
        for lines in open(mge,'r'):
            line_set = lines.split('\t',maxsplit=4)
            Gene1 = line_set[1]
            Gene2 = line_set[2]
            Gene_pair = function_com(Gene1, Gene2)
            if Gene_pair not in ['same','skip']:
                line_num += 1
                if line_num % 100000 == 0:
                    print('%s HGT_finder_sum processing %s lines for file %s' % (datetime.now(), line_num, mge))
                if "mge_" in Gene1 and "mge_" in Gene2:
                    # mge to mge
                    mge_to_mge += 1
                else:
                    # mge to genome
                    mge_to_genome += 1
    else:
        print('file %s is %s' % (mge, Checkoutput))
    # output
    Result = [function_name, type_fasta, "%.2f" % (cutoff),
              range16S_diff, hit_pair_diff, total_pair_diff, percentage_diff_pair, range16S_same,
              hit_pair_same, total_pair_same, percentage_same_pair, diff_same_ratio, mge_to_genome, mge_to_mge]
    for i in range(0,len(Result)):
        Result[i]=str(Result[i])
    output_file1.write('\t'.join(Result)+'\n')


################################################### Programme #######################################################
# set input
f16s = os.path.join(args.s, args.t + '.all.16S.fasta')
faa = os.path.join(args.s, args.t + '.all.traits.aa.fasta')
fdna = os.path.join(args.s, args.t + '.all.traits.dna.fasta')
fdna_500 = glob.glob(os.path.join(args.s, args.t + '.all.traits.dna.extra*.fasta'))[0]
ID_16S = dict()

# load traits search file for functions
Function_Set_dna=function_load(os.path.join(args.s, args.t + '.all.traits.dna.txt'),'dna')
Function_Set_aa=function_load(os.path.join(args.s, args.t + '.all.traits.aa.txt'),'aa')

# comparing genes
print('%s comparing and clustering 16S' % (datetime.now()))
cluster_16S = run_compare(f16s, Cutoff_16S,0.6,'dna','T')
print ('the range of length of 16S ID is %s to %s' %(cluster_16S[-1],cluster_16S[-2]))
print('%s comparing and clustering DNA' % (datetime.now()))
run_compare(fdna, Cutoff_HGT,Cutoff_HGT,'dna')
print('%s comparing and clustering DNA extended' % (datetime.now()))
run_compare(faa, Cutoff_aa,Cutoff_aa,'aa')
print('%s comparing and clustering AA' % (datetime.now()))
run_compare(fdna_500, Cutoff_extended,Cutoff_extended,'dna')

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


# filter usearch output into same 16S cluster and diff 16S clusters
print('%s compare %s trait identity to 16S' %(datetime.now(), 'dna'))
script_i = compare_traits_16S(Function_Set_dna, 'dna',args.s,
                              os.path.split(fdna)[-1] + '*.usearch.txt',Cutoff_HGT,script_i)
print('%s compare %s trait identity to 16S' %(datetime.now(), 'aa'))
script_i = compare_traits_16S(Function_Set_aa, 'aa',args.s,
                              os.path.split(faa)[-1] + '*.usearch.txt',Cutoff_aa,script_i)
print('%s compare %s trait identity to 16S' %(datetime.now(), 'extended dna'))
script_i = compare_traits_16S(Function_Set_dna, 'dna_extended',args.s,
                              os.path.split(fdna_500)[-1] + '*.usearch.txt',Cutoff_extended,script_i)

# calculate index for HGT
# load each function
all_function_diff_same = glob.glob(os.path.join(result_dir + '/sub_fun','*.diff.cluster'))+\
glob.glob(os.path.join(result_dir + '/sub_fun','*.same.cluster'))
all_output_file = os.path.join(result_dir, 'HGT.summary.dna.%s.aa.%s.16S.%s.txt'
                                   % (Cutoff_HGT, Cutoff_aa, Cutoff_16S))
all_output = open(all_output_file,'w')
all_output.write('function_name\ttype\tcutoff\trange16S_diff\thit_pair_diff\ttotal_pair_diff\t'+
                      'percentage_diff_pair\trange16S_same\thit_pair_same\ttotal_pair_same\tpercentage_same_pair\t'+
                      'diff_same_ratio\tmge_to_genome\tmge_to_mge\n')
all_output.close()
all_function = []
for function_diff in all_function_diff_same:
    function_name = os.path.split(function_diff)[-1].split('.',maxsplit=2)[0]
    if function_name not in all_function:
        all_function.append(function_name)
print(all_function,all_function_diff_same)
for function_name in all_function:
    # for each function
    function_diff = os.path.join(result_dir + '/sub_fun',"%s.%s.%s.diff.cluster" % (function_name,args.t,'dna'))
    function_same = os.path.join(result_dir + '/sub_fun',"%s.%s.%s.same.cluster" % (function_name,args.t,'dna'))
    function_mge = os.path.join(result_dir + '/sub_fun',"%s.%s.%s.mge.cluster" % (function_name,args.t,'dna'))
    print('%s summarize potential HGT of %s %s trait with cutoff of %s' % (datetime.now(), 'dna', function_name, Cutoff_HGT))
    all_output = open(all_output_file, 'a')
    HGT_finder_sum(function_name, 'dna',
                   Cutoff_HGT,
                   function_diff,
                   function_same,
                   function_mge,
                   all_output)
    all_output.close()
    all_output = open(all_output_file, 'a')
    print('%s summarize potential HGT of %s %s trait with cutoff of %s' % (datetime.now(), 'aa', function_name, Cutoff_aa))
    HGT_finder_sum(function_name, 'aa',
                   Cutoff_aa,
                   function_diff.replace('.dna.diff.cluster', '.aa.diff.cluster'),
                   function_same.replace('.dna.same.cluster', '.aa.same.cluster'),
                   function_mge.replace('.dna.diff.cluster', '.aa.mge.cluster'),
                   all_output)
    all_output.close()
    all_output = open(all_output_file, 'a')
    print('%s summarize potential HGT of %s %s trait with cutoff of %s' % (datetime.now(), 'extended dna', function_name, Cutoff_extended))
    HGT_finder_sum(function_name, 'dna_extended',
                   Cutoff_extended,
                   function_diff.replace('.dna.diff.cluster', '.dna_extended.diff.cluster'),
                   function_same.replace('.dna.same.cluster', '.dna_extended.same.cluster'),
                   function_mge.replace('.dna.diff.cluster', '.dna_extended.mge.cluster'),
                   all_output)
    all_output.close()
    all_output = open(all_output_file, 'a')
    print('%s summarize potential HGT of %s %s trait with cutoff of %s' % (datetime.now(), 'extended dna', function_name, Cutoff_HGT))
    HGT_finder_sum(function_name, 'dna_extended',
                   Cutoff_extended2,
                   function_diff.replace('.dna.diff.cluster', '.dna_extended.diff.cluster'),
                   function_same.replace('.dna.same.cluster', '.dna_extended.same.cluster'),
                   function_mge.replace('.dna.diff.cluster', '.dna_extended.mge.cluster'),
                   all_output)
    all_output.close()
# collect all bash files
list_of_files = glob.glob('HGT_subscripts/HGTalign.*.sh')
f1 = open("HGTalign.sh", 'w')
f1.write("#!/bin/bash\nsource ~/.bashrc\n")
for file_name in list_of_files:
    f1.write("jobmit %s HGTalign big\n" % (file_name))
f1.close()

# output mapping file
if mapping_file_output == 1:
    fout = open (mapping_file,'w')
    fout_set = []
    for genome in mapping:
        fout_set.append(genome+'\t'+mapping[genome]+'\t\n')
    fout.write(''.join(fout_set))
fout.close()
