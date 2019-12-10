import os
from Bio import SeqIO
import argparse
import glob

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-t",
                    help="trait name", type=str, default='ARG',metavar='trait name')
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
Cutoff_16S=1
Cutoff_HGT=0.99
Cutoff_aa=0.8
Cutoff_extended=0.5

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
        f1 = open("%s.sorted.uc" % (input_fasta),'r')
    except IOError:
        os.system('%s -sort length -cluster_fast %s -id %s -centroids %s.sorted.cluster.aa -uc %s.sorted.uc'
                  % (args.u, input_fasta, cutoff, input_fasta, input_fasta))
    try:
        f1 = open("%s.sorted" % (input_fasta), 'r')
    except IOError:
        os.system('%s -sortbylength %s -fastaout %s.sorted' % (args.u, input_fasta, input_fasta))
    input_fasta = input_fasta + '.sorted'
    try:
        f1 = open("%s.%s.usearch.txt" % (input_fasta, 0.6), 'r')
    except IOError:
        # need fix
        #os.system("%s -makeudb_usearch %s -output %s.udb\n"
        #                  % (args.u, input_fasta, input_fasta))
        #os.system("%s  -ublast %s -db %s.udb  -strand both -id %s -evalue 1 -acceptall -blast6out %s.%s.usearch.txt  -threads %s\n"
        #        % (args.u, input_fasta, input_fasta, 0.6, input_fasta, 0.6, str(args.th)))
        os.system('%s windowmasker -in %s -infmt blastdb -mk_counts -out %s.counts' %
                  ('hs-blastn',input_fasta,input_fasta))
        os.system('%s windowmasker -in %s.counts -sformat obinary -out %s.counts.obinary -convert' %
                  ('hs-blastn',input_fasta,input_fasta))
        os.system('%s index %s' % ('hs-blastn',input_fasta))
        os.system( "%s align -db %s -window_masker_db %s.counts.obinary -query %s -out %s.%s.usearch.txt -outfmt 6 -evalue 1 -perc_identity %s -num_threads %s\n" \
                    % ('hs-blastn', input_fasta,
                       input_fasta, input_fasta, input_fasta,0.6,
                       0.6, str(args.th)))
    print('finish running usearch cluster for %s' % (input_fasta))
    # read cluster results
    Clusters = dict()
    Clusters_seqs = dict()
    for lines in open(input_fasta + '.uc', 'r'):
        cluster = lines.split('\t')[1]
        record_name = lines.split('\t')[-2].split(' ')[0]
        if cluster not in Clusters:
            Clusters.setdefault(cluster, [record_name])
        elif record_name not in Clusters[cluster]:
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
        if gene not in Function_Set:
            Function_Set.setdefault(gene,[[function,[int(lines.split('\t')[7]),int(lines.split('\t')[8])]]])
        else:
            Function_Set[gene].append([function,[int(lines.split('\t')[7]),int(lines.split('\t')[8])]])
    return Function_Set


def find_genome(Genome1,cluster_16S_seqs):
    for i in range(0, len(Genome1.split('_')) + 1):
        if '_'.join(Genome1.split('_')[0:i]) in cluster_16S_seqs:
            return '_'.join(Genome1.split('_')[0:i])
        if '_'.join(Genome1.split('_')[0:i]).split('.')[0] in cluster_16S_seqs:
            return '_'.join(Genome1.split('_')[0:i]).split('.')[0]
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
            diff_16S.setdefault(Genome1+'_'+Genome2,0)
        return Cluster1 != Cluster2
    else:
        return '16S missing'


def loci_seq(record_name):
    loci1 = int(record_name.split('_')[-2])
    loci2 = int(record_name.split('_')[-1])
    return [loci1,loci2]


def compare_loci(loci_new,loci_ref):
    if min(loci_new[0],loci_new[1]) <= min(loci_ref[0],loci_ref[1]) and\
        max(loci_new[0], loci_new[1]) >= max(loci_ref[0], loci_ref[1]):
        return 1
    else:
        return 0


def function_find(Function_Set,Genome):
    loci_new=loci_seq(Genome)
    Genome_name='_'.join(Genome.split('_')[0:-2])
    for functions in Function_Set.get(Genome_name):
        loci_ref=functions[-1]
        if compare_loci(loci_new, loci_ref) == 1:
            return functions[0]

def function_com(function1, function2):
    if function1 == function2:
        return function1
    else:
        temp = [function1,function2]
        temp.sort()
        return '%s_%s'%(temp[0],temp[1])


def function_pair(Function_Set,Genome1,Genome2,type_fasta):
    if type_fasta == 'aa':
        function1=Function_Set.get(Genome1)[0][0]
        function2 = Function_Set.get(Genome2)[0][0]
    else:
        function1 = function_find(Function_Set, Genome1)
        function2 = function_find(Function_Set, Genome2)
    return function_com(function1, function2)


def compare_traits_16S(Function_Set,cluster_16S_seqs,type_fasta,input_folder,input_prefix,cutoff):
    all_usearch = glob.glob(os.path.join(input_folder, input_prefix))
    for files in all_usearch:
        for lines in open(files,'r'):
            if float(str(lines).split('\t')[2]) >= cutoff:
                Genome1 = lines.split('\t')[0]
                Genome2 = lines.split('\t')[1]
                # not the same gene
                if Genome1 != Genome2:
                    compare_result= compare_16S(Genome1,Genome2,cluster_16S_seqs)
                    Function=function_pair(Function_Set,Genome1,Genome2,type_fasta)
                    cluster=int(os.path.split(files)[-1].split('.fasta.sorted')[0].split('.')[-1])
                    if compare_result != '16S missing':
                        if compare_result:
                            # different 16S clusters
                            fout=open(os.path.join(result_dir + '/sub_fun',
                                                   "%s.%s.%s.diff.cluster" %
                                                   (Function,args.t,type_fasta)),'a')
                            fout.write(Function+'\t'+str(cluster)+'\t'+lines)
                            fout.close()
                        else:
                            # same 16S cluster
                            fout = open(os.path.join(result_dir + '/sub_fun',
                                                   "%s.%s.%s.same.cluster" %
                                                   (Function,args.t,type_fasta)),'a')
                            fout.write(Function+'\t'+str(cluster)+'\t'+lines)
                            fout.close()


def usearch_16S_load(input_file):
    for lines in open(input_file,'r'):
        print(lines)
        if lines.split('\t')[0]+'_'+lines.split('\t')[1] in diff_16S:
            diff_16S[lines.split('\t')[0]+'_'+lines.split('\t')[1]]=float(lines.split('\t')[2])
        elif lines.split('\t')[1]+'_'+lines.split('\t')[0] in diff_16S:
            diff_16S[lines.split('\t')[1]+'_'+lines.split('\t')[0]]=float(lines.split('\t')[2])


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
        Gene_set=[]
        for lines in fdiff:
            Gene1 = lines.split('\t')[2]
            Gene2 = lines.split('\t')[3]
            Gene_pair=function_com(Gene1, Gene2)
            if Gene_pair not in Gene_set:
                Gene_set.append(Gene_pair)
                ID = float(lines.split('\t')[4])
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
                    if Genome1 +'_'+Genome2 in diff_16S:
                        lowest_id = diff_16S[Genome1 +'_'+Genome2]
                        Diff_16S_min = min(float(lowest_id),Diff_16S_min)
                    elif Genome2 +'_'+Genome1 in diff_16S:
                        Diff_16S_min = min(float(diff_16S[Genome2 +'_'+Genome1]),Diff_16S_min)
                    output_file2.write('%s\t%s\t%s\n'
                                           % (Genome_pair,str(ID),lowest_id))
                    if '.dna_extended.' in diff:
                        if ID >= 0.9:
                            output_file3.write('%s\t%s\t%s\n'
                                               % (Genome_pair, str(ID), lowest_id))
        if '.dna_extended.' in diff:
            output_file3.close()
        output_file2.close()
    except IOError:
        pass
    # calculate same 16S cluster
    try:
        Gene_set = []
        for lines in open(same,'r'):
            Gene1 = lines.split('\t')[2]
            Gene2 = lines.split('\t')[3]
            Gene_pair = function_com(Gene1, Gene2)
            if Gene_pair not in Gene_set:
                Gene_set.append(Gene_pair)
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
                        if cluster1 == cluster2:
                            Result[11] += len(cluster_16S[-1][cluster1])
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
    try:
        # number of hits / total number of combination of random 2 genomes
        Result[7] = Result[5]/(Result[6]*(Result[6]-1)/2)
    except ZeroDivisionError:
        Result[7] = 0
    Result[8] = len(Same_cluster)
    Result[10] = len(Same_genome_set)
    try:
        # number of hits / total number of combination of random 2 genomes
        Result[12] = Result[10] /(Result[11]*(Result[11]-1)/2)
    except ZeroDivisionError:
        Result[12] = 0
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
Function_Set_aa=function_load(os.path.join(args.s, args.t + '.all.traits.aa.txt'))

# cluster 16S, by 3% similarity
cluster_16S = run_16S(f16s, Cutoff_16S)
print(cluster_16S)
# filter usearch output into same 16S cluster and diff 16S clusters
compare_traits_16S(Function_Set_dna, cluster_16S[0],'dna',result_dir +
                              '/sub_sequences',os.path.split(fdna)[-1] + '*.usearch.txt',Cutoff_HGT)
compare_traits_16S(Function_Set_aa, cluster_16S[0],'aa',result_dir +
                             '/sub_sequences',os.path.split(faa)[-1] + '*.usearch.txt',Cutoff_aa)
compare_traits_16S(Function_Set_dna, cluster_16S[0],'dna_extended',result_dir +
                             '/sub_sequences',os.path.split(fdna_500)[-1] + '*.usearch.txt',Cutoff_extended)
# load diff 16S usearch result
usearch_16S_load("%s.sorted.%s.usearch.txt" % (f16s,0.6))
print(diff_16S)
# calculate index for HGT
# load each function
all_function_diff_same = glob.glob(os.path.join(result_dir + '/sub_fun','*dna.diff.cluster'))+\
glob.glob(os.path.join(result_dir + '/sub_fun','*dna.same.cluster'))
all_output = open(os.path.join(result_dir,'HGT.summary.dna.%s.aa.%s.16S.%s.txt'
                                   % (Cutoff_HGT,Cutoff_aa,Cutoff_16S)),'w')
all_function = []
for function_diff in all_function_diff_same:
    # for each function
    function_diff = function_diff.replace('.same.cluster', '.diff.cluster')
    function_same=function_diff.replace('.diff.cluster','.same.cluster')
    function_name = os.path.split(function_diff)[-1].split('.')[0]
    print(function_name)
    if function_name not in all_function:
        HGT_finder_sum(cluster_16S,function_name,'dna',
                       Cutoff_HGT,function_diff,function_same,
                       all_output)
        HGT_finder_sum(cluster_16S, function_name, 'aa',
                       Cutoff_aa,
                       function_diff.replace('.dna.diff.cluster', '.aa.diff.cluster'),
                       function_same.replace('.dna.same.cluster', '.aa.same.cluster'),
                       all_output)
        HGT_finder_sum(cluster_16S, function_name, 'dna_extended',
                       Cutoff_extended,
                       function_diff.replace('.dna.diff.cluster', '.dna_extended.diff.cluster'),
                       function_same.replace('.dna.same.cluster', '.dna_extended.same.cluster'),
                       all_output)
        all_function.append(function_name)
all_output.close()

# check an essential and check an ARG (mcr-1)