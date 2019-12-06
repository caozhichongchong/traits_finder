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
parser.add_argument("-m",
                    help="mapping file of traits to function", type=str,
                    default='Butyrate.pro.mapping.txt',
                    metavar='Butyrate.pro.mapping.txt')
# optional parameters
parser.add_argument("--fa",
                    help="input format of genome sequence", type=str,
                    default='.fna.add',metavar='.fasta, .fna or .fa')
parser.add_argument("--orf",
                    help="input format of genome orfs", type=str,
                    default='.genes.faa',metavar='.faa')
# optional input setup
parser.add_argument("--r",
                    help="input directory or folder of your previous results by Traits_WGD.py",
                    type=str, default='Result',metavar='Result')
parser.add_argument("--r16",
                    help="input directory or folder of your previous 16S sequences extracted by Traits_WGD.py",
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
                      help="Optional: use two-step method for blast search," +
                           " \'None\' for using one step, \'usearch\' or \'diamond\' for using two-step \
                           (complete path to usearch or diamond if not in PATH, \
                           please make sure the search tools can be directly called), (default: \'None\')",
                      metavar="None or usearch",
                      action='store', default='None', type=str)
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
fasta_format = args.fa
orfs_format = args.orf
Cutoff_16S=0.97
Cutoff_HGT=0.99
Cutoff_aa=0.9
Cutoff_extended=0.1

if args.s == 'None':
    input_dir = os.path.join(args.r,'summary')
    result_dir = os.path.join(args.r, 'HGT')
else:
    input_dir = args.s
    result_dir = os.path.join(args.s, '../HGT')
try:
    os.mkdir(result_dir + '/sub_fun')
except OSError:
    pass
workingdir=os.path.abspath(os.path.dirname(__file__))


################################################### Function ########################################################
def run_16S(input_fasta,cutoff=0.97):
    print('Running usearch cluster for %s' % (input_fasta))
    try:
        f1 = open("%s.uc" % (input_fasta),'r')
    except IOError:
        os.system('%s -sort length -cluster_fast %s -id %s -centroids %s.cluster.aa -uc %s.uc'
                  % (args.u, input_fasta, cutoff, input_fasta, input_fasta))
    try:
        f1 = open("%s.sorted" % (input_fasta), 'r')
    except IOError:
        os.system('%s -sortbylength %s -fastaout %s.sorted' % (args.u, input_fasta, input_fasta))
    input_fasta = input_fasta + '.sorted'
    try:
        f1 = open("%s.%s.usearch.txt" % (input_fasta, str(cutoff)), 'r')
    except IOError:
        os.system("%s -makeudb_usearch %s -output %s.udb\n"
                          % (args.u, input_fasta, input_fasta))
        os.system("%s  -ublast %s -db %s.udb  -strand both -id %s -evalue 1e-2 -accel 0.5 -blast6out %s.%s.usearch.txt  -threads %s\n"
                % (args.u, input_fasta, input_fasta, 0.6, input_fasta, 0.6, str(args.th)))
    print('finish running usearch cluster for %s' % (input_fasta))
    # read cluster results
    Clusters = dict()
    Clusters_seqs = dict()
    for lines in open(input_fasta + '.rc', 'r'):
        cluster = lines.split('\t')[1]
        record_name = lines.split('\t')[-2].split(' ')[0]
        if cluster not in Clusters:
            Clusters.setdefault(cluster, [record_name])
        else:
            Clusters[cluster].append(record_name)
    for cluster in Clusters:
        for record_name in Clusters[cluster]:
            Clusters_seqs.setdefault(record_name, str(cluster))
    return [Clusters_seqs,Clusters]


def function_load(input_file):
    Function_Set = dict()
    for lines in open(input_file,'r'):
        function = lines.split('\t')[0]
        gene = lines.split('\t')[1]
        Function_Set.setdefault(gene,function)
    return Function_Set


def find_genome(Genome1,cluster_16S_seqs):
    for i in range(0, len(Genome1.split('_')) + 1):
        if '_'.join(Genome1.split('_')[0:i]) in cluster_16S_seqs:
            return '_'.join(Genome1.split('_')[0:i])
    else:
        return 'None'


def compare_16S(Genome1,Genome2,cluster_16S_seqs):
    if 'reference' not in Genome1 and 'reference' not in Genome2:
        Genome1 = find_genome(Genome1,cluster_16S_seqs)
        Genome2 = find_genome(Genome2,cluster_16S_seqs)
    if Genome1 != 'None' and Genome2 != 'None':
        Cluster1 = cluster_16S_seqs[Genome1]
        Cluster2 = cluster_16S_seqs[Genome2]
        if Cluster1 != Cluster2:
            diff_16S.setdefault([Genome1,Genome2],0)
        return Cluster1 != Cluster2
    else:
        return '16S missing'


def function_pair(function1,function2):
    if function1 == function2:
        return function1
    else:
        temp = [function1,function2]
        temp.sort()
        return '_'.join(temp)


def compare_traits_16S(Function_Set,cluster_16S_seqs,type_fasta,input_folder,input_prefix,cutoff):
    all_usearch = glob.glob(os.path.join(input_folder, input_prefix))
    for files in open(all_usearch,'r'):
        for lines in files:
            if float(str(lines).split('\t')[2]) >= cutoff:
                Genome1 = lines.split('\t')[0]
                Genome2 = lines.split('\t')[1]
                compare_result= compare_16S(Genome1,Genome2,cluster_16S_seqs)
                Function=function_pair(Function_Set.get(Genome1),Function_Set.get(Genome2))
                cluster=int(os.path.split(files)[-1].split('.usearch.txt')[-1].split('.')[-1])
                if compare_result != '16S missing':
                    if compare_result:
                        # different 16S clusters
                        fout=open(os.path.join(result_dir + '/sub_fun',
                                               "%s.%s.%s.diff.cluster" %
                                               (Function,type_fasta,os.path.split(files)[-1])),'a')
                        fout.write(Function+'\t'+cluster+'\t'+lines)
                        fout.close()
                    else:
                        # same 16S cluster
                        fout = open(os.path.join(result_dir + '/sub_fun',
                                               "%s.%s.%s.same.cluster" %
                                               (Function,type_fasta,os.path.split(files)[-1])),'a')
                        fout.write(Function+'\t'+cluster+'\t'+lines)
                        fout.close()


def usearch_16S_load(input_file):
    for lines in open(input_file,'r'):
        if [lines.split('\t')[0],lines.split('\t')[1]] in diff_16S:
            diff_16S[[lines.split('\t')[0],lines.split('\t')[1]]]=float(lines.split('\t')[2])
        elif [lines.split('\t')[1],lines.split('\t')[0]] in diff_16S:
            diff_16S[[lines.split('\t')[1],lines.split('\t')[0]]]=float(lines.split('\t')[2])


def HGT_finder_sum(cluster_16S,function_name,type_fasta,cutoff,diff,same,output_file1):
    # for each function
    output_file1.write('function_name\ttype\tcutoff\tcluster_num_diff\t16S_range_diff\thit_diff\tgenome_num_diff\t'+
                      'percentage_diff\tcluster_num_same\t16S_range_same\thit_same\tgenome_num_same\tpercentage_same\n')
    Result=[function_name,type_fasta,cutoff,
            0,'',0,0,0,
            0,'%s-1.00'%(str(Cutoff_16S)),0,0,0]
    Diff_cluster=[]
    Diff_16S_min = Cutoff_16S
    Same_cluster=[]
    Diff_genome_set=[]
    Same_genome_set = []
    # calculate diff 16S clusters
    try:
        fdiff = open(diff,'r')
        output_file2 = open(os.path.join(result_dir,'%s.%s.%s.identity.summary.txt'
                                   % (function_name,type_fasta,cutoff)),'a')
        output_file2.write('genome_pair\tid_gene\tid_16S\n')
        if '.dna_extended.' in diff:
            output_file3 = open(os.path.join(result_dir,'%s.%s.%s.identity.summary.txt'
                                   % (function_name,type_fasta,'0.9')),'a')
            output_file3.write('genome_pair\tid_gene\tid_16S\n')
        for lines in fdiff:
            Gene1 = lines.split('\t')[2]
            Gene2 = lines.split('\t')[3]
            ID = float(lines.split('\t')[4])
            Genome1 = find_genome(Gene1, cluster_16S[0])
            cluster1 = cluster_16S[0][Genome1]
            Genome2 = find_genome(Gene2, cluster_16S[0])
            cluster2 = cluster_16S[0][Genome2]
            if Genome1 != Genome2:
                # count genome pairs
                Genome_pair = function_pair(Genome1, Genome2)
                if Genome_pair not in Diff_genome_set:
                    Diff_genome_set.append(Genome_pair)
                # calculate total number of 16S
                if cluster1 != cluster2:
                    Result[6] += len(cluster_16S[-1][cluster1])+len(cluster_16S[-1][cluster2])
                else:
                    print('Error of compare_traits_16S: %s and %s should have different 16S clusters'
                          %(Gene1,Gene2))
                # calculate total number of gene clusters
                if lines.split('\t')[1] not in Diff_cluster:
                    Diff_cluster.append(lines.split('\t')[1])
                # calculate lowest 16S similarity for same gene in diff 16S clusters
                lowest_id = Cutoff_16S
                if [Genome1, Genome2] in diff_16S:
                    lowest_id = diff_16S[Genome1, Genome2]
                    Diff_16S_min = min(float(lowest_id),Diff_16S_min)
                elif [Genome2, Genome1] in diff_16S:
                    Diff_16S_min = min(float(diff_16S[Genome2, Genome1]),Diff_16S_min)
                output_file2.write('%s\tid_gene\tid_16S\n'
                                       % (Genome_pair,ID,lowest_id))
                if '.dna_extended.' in diff:
                    if ID >= 0.9:
                        output_file3.write('%s\tid_gene\tid_16S\n'
                                           % (Genome_pair, ID, lowest_id))
        if '.dna_extended.' in diff:
            output_file3.close()
        output_file2.close()
    except IOError:
        pass
    # calculate same 16S cluster
    try:
        for lines in open(same,'r'):
            Gene1 = lines.split('\t')[2]
            Gene2 = lines.split('\t')[3]
            Genome1 = find_genome(Gene1, cluster_16S[0])
            cluster1 = cluster_16S[0][Genome1]
            Genome2 = find_genome(Gene2, cluster_16S[0])
            cluster2 = cluster_16S[0][Genome2]
            # calculate total number of 16S
            if Genome1 != Genome2:
                # count genome pairs
                Genome_pair = function_pair(Genome1, Genome2)
                if Genome_pair not in Same_genome_set:
                    Same_genome_set.append(Genome_pair)
                # calculate total number of 16S
                if cluster1 == cluster2:
                    Result[10] += len(cluster_16S[-1][cluster1]) + len(cluster_16S[-1][cluster2])
                else:
                    print('Error of compare_traits_16S: %s and %s should have the same 16S cluster'
                          % (Gene1, Gene2))
            # calculate total number of gene clusters
            if lines.split('\t')[1] not in Same_cluster:
                Same_cluster.append(lines.split('\t')[1])
    except IOError:
        pass
    # summarize
    Result[3] = len(Diff_cluster)
    Result[4] = '%s-%s'%(str(Diff_16S_min),str(Cutoff_16S))
    Result[5] = len(Diff_genome_set)
    Result[7] = Result[5]/Result[6]
    Result[8] = len(Same_cluster)
    Result[11] = len(Same_genome_set)
    Result[12] = Result[10] / Result[11]
    # output
    for i in range(0,len(Result)):
        Result[i]=str(Result[i])
    output_file1.write('\t'.join(Result)+'\n')


################################################### Programme #######################################################
f16s = os.path.join(args.s, args.t + '.all.16S.fasta')
faa = os.path.join(args.s, args.t + '.all.traits.aa.fasta')
fdna = os.path.join(args.s, args.t + '.all.traits.dna.fasta')
fdna_500 = glob.glob(os.path.join(args.s, args.t + '.all.traits.dna.extra*.fasta'))[0]
diff_16S = dict()

# load traits search file for functions
Function_Set_dna=function_load(os.path.join(args.s, args.t + '.all.traits.dna.txt'))
Function_Set_aa=function_load(os.path.join(args.s, args.t + '.all.traits.dna.txt'))

# cluster 16S, by 3% similarity
cluster_16S = run_16S(f16s, Cutoff_16S)

# filter usearch output into same 16S cluster and diff 16S clusters
compare_traits_16S(Function_Set_dna, cluster_16S[0],'dna',result_dir +
                              '/sub_sequences',os.path.split(fdna)[-1] + '*.usearch.txt',Cutoff_HGT)
compare_traits_16S(Function_Set_aa, cluster_16S[0],'aa',result_dir +
                             '/sub_sequences',os.path.split(faa)[-1] + '*.usearch.txt',Cutoff_aa)
compare_traits_16S(Function_Set_aa, cluster_16S[0],'dna_extended',result_dir +
                             '/sub_sequences',os.path.split(faa)[-1] + '*.usearch.txt',Cutoff_extended)
# load diff 16S usearch result
usearch_16S_load("%s.%s.usearch.txt" % (f16s,0.6))

# calculate index for HGT
# load each function
all_function_diff = glob.glob(os.path.join(result_dir + '/sub_fun','*.dna.diff.cluster'))
all_output = open(os.path.join(result_dir,'HGT.summary.dna.%s.aa.%s.16S.%stxt'
                                   % (Cutoff_HGT,Cutoff_aa,Cutoff_16S)),'w')
for function_diff in all_function_diff:
    # for each function
    function_same=function_diff.replace('.diff.cluster','.same.cluster')
    function_name = os.path.split(function_diff)[-1].split('.dna')[0]
    HGT_finder_sum(cluster_16S,function_name,'dna',
                   Cutoff_HGT,function_diff,function_same,
                   all_output)
    HGT_finder_sum(cluster_16S, function_name, 'aa',
                   Cutoff_HGT,
                   function_diff.replace('.dna.diff.cluster', '.aa.diff.cluster'),
                   function_same.replace('.dna.same.cluster', '.aa.same.cluster'),
                   all_output)
    HGT_finder_sum(cluster_16S, function_name, 'dna',
                   Cutoff_HGT,
                   function_diff.replace('.dna.diff.cluster', '.dna_extended.diff.cluster'),
                   function_same.replace('.dna.same.cluster', '.dna_extended.same.cluster'),
                   all_output)
all_output.close()

# check an essential and check an ARG (mcr-1)