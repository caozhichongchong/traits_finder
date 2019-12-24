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
Cutoff_extended=0.5
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
                    lines.split('\t')[i]
                    return 'not empty'
                except IndexError:
                    return 'wrong content by spliting %s \\t' % (str(i))
                break
        else:
            return 'empty'
    except IOError:
        return 'non-existed'


def usearch_16S_load(input_file):
    Checkoutput = checkfile(input_file, 2)
    if Checkoutput == 'not empty':
        for lines in open(input_file,'r'):
            if lines.split('\t')[0] + '_' + lines.split('\t')[1] not in ID_16S and \
                    lines.split('\t')[1] + '_' + lines.split('\t')[0] not in ID_16S:
                ID_16S.setdefault(lines.split('\t')[0] + '_' + lines.split('\t')[1], float(lines.split('\t')[2])/100.0)
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
                Gene1 = lines.split('\t')[0]
                Gene2 = lines.split('\t')[1]
                ID = float(lines.split('\t')[2])/100.0
                if Gene1 + '_' + Gene2 not in ID_16S and \
                        Gene2 + '_' + Gene1 not in ID_16S:
                    ID_16S.setdefault(Gene1 + '_' + Gene2,ID)
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


def run_16S(input_fasta,cutoff=0.97):
    try:
        f1 = open("%s.sorted" % (input_fasta), 'r')
    except IOError:
        os.system('%s -sortbylength %s -fastaout %s.sorted' % (args.u, input_fasta, input_fasta))
    input_fasta = input_fasta + '.sorted'
    try:
        f1 = open("%s.%s.usearch.txt" % (input_fasta, 0.6), 'r')
    except IOError:
        print('Running usearch for %s' % (input_fasta))
        if int(os.path.getsize(input_fasta)) <= 2 * 1E+9:
            # smaller than 2G
            os.system("%s -makeudb_usearch %s -output %s.udb\n"
                      % (args.u, input_fasta, input_fasta))
            os.system(
                "%s  -usearch_global %s -db %s.udb  -strand both -id %s -maxaccepts 0 -maxrejects 0 -blast6out %s.%s.usearch.txt  -threads %s\n"
                % (args.u, input_fasta, input_fasta, 0.6, input_fasta, 0.6, str(args.th)))
        else:
            print('Using hs-blastn instead of usearch because the input file is larger than 2GB\n')
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
                   input_fasta, input_fasta, input_fasta, 0.6,
                   0.6, str(min(int(args.th),40))))
    try:
        f1 = open("%s.uc" % (input_fasta),'r')
    except IOError:
        print('Running usearch cluster for %s' % (input_fasta))
        # smaller than 2G
        if int(os.path.getsize(input_fasta)) <= 2 * 1E+9:
            os.system('%s -sort length -cluster_fast %s -id %s -centroids %s.cluster.aa -uc %s.uc -threads %s'
                      % (args.u, input_fasta, cutoff, input_fasta, input_fasta, args.th))
        else:
            print('We roughly clustered the 16S by 99% identity\n')
            self_clustering('%s.%s.usearch.txt' % (input_fasta, 0.6), input_fasta + '.uc')
    # load 16S usearch result
    if ID_16S == dict():
        usearch_16S_load("%s.%s.usearch.txt" % (input_fasta, 0.6))
    print('finish running usearch cluster for %s' % (input_fasta))
    # read cluster results
    Clusters = dict()
    Clusters_seqs = dict()
    Checkoutput = checkfile(input_fasta + '.uc', -2)
    if Checkoutput == 'not empty':
        for lines in open(input_fasta + '.uc', 'r'):
                cluster = lines.split('\t')[1]
                record_name = lines.split('\t')[-2].split(' ')[0]
                if cluster not in Clusters:
                    Clusters.setdefault(cluster, [record_name])
                elif record_name not in Clusters[cluster]:
                    Clusters[cluster].append(record_name)
    else:
        print('file %s is %s' % (input_fasta + '.uc', Checkoutput))
    for cluster in Clusters:
        for record_name in Clusters[cluster]:
            Clusters_seqs.setdefault(record_name, str(cluster))
    return [Clusters_seqs,Clusters]


def function_load(input_file):
    Function_Set = dict()
    Checkoutput = checkfile(input_file, 8)
    if Checkoutput == 'not empty':
        for lines in open(input_file,'r'):
                function = lines.split('\t')[0]
                gene = lines.split('\t')[1]
                # query gene
                if gene not in Function_Set:
                    Function_Set.setdefault(gene,[[function,[int(lines.split('\t')[7]),int(lines.split('\t')[8])]]])
                else:
                    Function_Set[gene].append([function,[int(lines.split('\t')[7]),int(lines.split('\t')[8])]])
    else:
        print('file %s is %s' % (input_file, Checkoutput))
    return Function_Set


def find_genome(Genome1,cluster_16S_seqs):
    for i in range(0, len(Genome1.split('_')) + 1):
        if '_'.join(Genome1.split('_')[0:i]) in cluster_16S_seqs:
            return '_'.join(Genome1.split('_')[0:i])
        if '_'.join(Genome1.split('_')[0:i]).split('.')[0] in cluster_16S_seqs:
            return '_'.join(Genome1.split('_')[0:i]).split('.')[0]
    else:
        return 'None'


def compare_16S(Genome1,Genome2,cluster_16S_seqs,cutoff):
    if 'reference' not in Genome1 and 'reference' not in Genome2\
            and 'mge' not in Genome1 and 'mge' not in Genome2:
        Genome1 = find_genome(Genome1,cluster_16S_seqs)
        Genome2 = find_genome(Genome2,cluster_16S_seqs)
        if Genome1 != 'None' and Genome2 != 'None':
            if Genome1 + '_' + Genome2 in ID_16S:
                if ID_16S[Genome1 + '_' + Genome2] < cutoff:
                    return True
                else:
                    return False
            elif Genome2 + '_' + Genome1 in ID_16S:
                if ID_16S[Genome2 + '_' + Genome1] < cutoff:
                    return True
                else:
                    return False
            else:
                return '16S missing'
    elif "mge_" in Genome1 or "mge_" in Genome2:
        return "mge"
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


def function_find(Function_Set, Genome, type_fasta):
    if "reference" not in Genome:
        # query genes
        if type_fasta == 'aa':
            return Function_Set.get(Genome)[0][0]
        else:
            loci_new=loci_seq(Genome)
            Genome_name='_'.join(Genome.split('_')[0:-2])
            for functions in Function_Set.get(Genome_name):
                loci_ref=functions[-1]
                if compare_loci(loci_new, loci_ref) == 1:
                    return functions[0]
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
    if Function not in Diff_gene_set:
        Diff_gene_set.setdefault(Function, [Gene])
    elif Gene not in Diff_gene_set[Function]:
        Diff_gene_set[Function].append(Gene)


def compare_traits_16S(Function_Set,cluster_16S_seqs,type_fasta,input_folder,input_prefix,cutoff,script_i):
    all_usearch = glob.glob(os.path.join(input_folder, input_prefix))
    for files in all_usearch:
        Diff_gene_set = dict()
        Checkoutput = checkfile(files, 2)
        if Checkoutput == 'not empty':
            for lines in open(files,'r'):
                try:
                    if float(str(lines).split('\t')[2])/100.0 >= cutoff:
                        Genome1 = lines.split('\t')[0]
                        Genome2 = lines.split('\t')[1]
                        Function = function_pair(Function_Set, Genome1, Genome2, type_fasta)
                        cluster = int(os.path.split(files)[-1].split('.fasta.sorted')[0].split('.')[-1])
                        # not the same gene
                        if Genome1 != Genome2 and "reference" not in Genome1 and "reference" not in Genome2:
                            compare_result= compare_16S(Genome1,Genome2,cluster_16S_seqs,Cutoff_16S)
                            if compare_result != '16S missing':
                                if compare_result != 'mge':
                                    if compare_result:
                                        # different 16S clusters
                                        fout=open(os.path.join(result_dir + '/sub_fun',
                                                               "%s.%s.%s.diff.cluster" %
                                                               (Function,args.t,type_fasta)),'a')
                                        fout.write(Function+'\t'+str(cluster)+'\t'+lines)
                                        fout.close()
                                        # record diff gene set
                                        Gene1 = lines.split('\t')[0]
                                        Gene2 = lines.split('\t')[1]
                                        if Gene1 not in Diff_gene_set:
                                            add_gene_and_function(Diff_gene_set, Function, Gene1)
                                        if Gene2 not in Diff_gene_set:
                                            add_gene_and_function(Diff_gene_set, Function, Gene2)
                                    else:
                                        # same 16S cluster
                                        fout = open(os.path.join(result_dir + '/sub_fun',
                                                               "%s.%s.%s.same.cluster" %
                                                               (Function,args.t,type_fasta)),'a')
                                        fout.write(Function+'\t'+str(cluster)+'\t'+lines)
                                        fout.close()
                                else:
                                    # different 16S clusters
                                    fout=open(os.path.join(result_dir + '/sub_fun',
                                                           "%s.%s.%s.mge.cluster" %
                                                           (Function,args.t,type_fasta)),'a')
                                    fout.write(Function+'\t'+str(cluster)+'\t'+lines)
                                    fout.close()
                        if "reference" in Genome1:
                            add_gene_and_function(Diff_gene_set, Function, Genome1)
                        if "reference" in Genome2:
                            add_gene_and_function(Diff_gene_set, Function, Genome2)
                except IndexError:
                    print('file %s is %s' % (files, 'wrong content by spliting %s \\t' % ('2')))
                    print(lines)
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


def HGT_finder_sum(cluster_16S,function_name,type_fasta,cutoff,diff,same,mge,output_file1):
    # for each function
    Result=[function_name,type_fasta,"%.2f"%(cutoff),
            0,'',0,0,0,
            0,'%.3f-1.000'%((Cutoff_16S)),0,0,0,
            0,0]
    Diff_cluster=[]
    Diff_16S_min = Cutoff_16S
    Same_cluster=[]
    Same_genome_set = []
    Diff_genome_set = []
    # calculate diff 16S clusters
    Checkoutput = checkfile(diff, 4)
    if Checkoutput == 'not empty':
        output_file2 = open(os.path.join(result_dir,'%s.%s.%.2f.identity.summary.txt'
                                   % (function_name,type_fasta,cutoff)),'a')
        output_file2.write('genome_pair\tid_gene\tid_16S\n')
        if '.dna_extended.' in diff:
            output_file3 = open(os.path.join(result_dir,'%s.%s.%s.identity.summary.txt'
                                   % (function_name,type_fasta,'0.90')),'a')
            output_file3.write('genome_pair\tid_gene\tid_16S\n')
        Gene_set=[]
        Cluster_16S_Set=[]
        for lines in open(diff,'r'):
            Gene1 = lines.split('\t')[2]
            Gene2 = lines.split('\t')[3]
            Gene_pair=function_com(Gene1, Gene2)
            if Gene_pair not in Gene_set and "reference" not in Gene_pair:
                Gene_set.append(Gene_pair)
                ID = float(lines.split('\t')[4])/100.0
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
                        if cluster1 not in Cluster_16S_Set:
                            Result[6] += len(cluster_16S[-1][cluster1])
                            Cluster_16S_Set.append(cluster1)
                        if cluster1 != cluster2:
                            if cluster2 not in Cluster_16S_Set:
                                Result[6] += len(cluster_16S[-1][cluster2])
                                Cluster_16S_Set.append(cluster2)
                        else:
                            pass
                            #print('%s and %s have the same 16S clusters but < %s 16S similarity'
                            #      %(Gene1,Gene2,Cutoff_16S))
                    # calculate total number of gene clusters
                    if lines.split('\t')[1] not in Diff_cluster:
                        Diff_cluster.append(lines.split('\t')[1])
                    # calculate lowest 16S similarity for same gene in diff 16S clusters
                    lowest_id = Cutoff_16S
                    if Genome1 +'_'+Genome2 in ID_16S:
                        lowest_id = ID_16S[Genome1 +'_'+Genome2]
                    elif Genome2 +'_'+Genome1 in ID_16S:
                        lowest_id = ID_16S[Genome2 +'_'+Genome1]
                    Diff_16S_min = min(float(lowest_id),Diff_16S_min)
                    output_file2.write('%s\t%.3f\t%.3f\n'
                                           % (Genome_pair,(ID),lowest_id))
                    if '.dna_extended.' in diff:
                        if ID >= 0.9:
                            output_file3.write('%s\t%.3f\t%.3f\n'
                                               % (Genome_pair, (ID), lowest_id))
        if '.dna_extended.' in diff:
            output_file3.close()
        output_file2.close()
    else:
        print('file %s is %s' % (diff, Checkoutput))
    # calculate same 16S cluster
    Checkoutput = checkfile(same, 3)
    if Checkoutput == 'not empty':
        Gene_set = []
        Cluster_16S_Set = []
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
                        if cluster1 not in Cluster_16S_Set:
                            Result[11] += len(cluster_16S[-1][cluster1])
                            Cluster_16S_Set.append(cluster1)
                        if cluster1 == cluster2:
                            pass
                        else:
                            if cluster2 not in Cluster_16S_Set:
                                Result[11] += len(cluster_16S[-1][cluster2])
                                Cluster_16S_Set.append(cluster2)
                            #print('%s and %s have different 16S clusters but >= %s 16S similarity'
                            #      % (Gene1, Gene2, Cutoff_16S))
                    # calculate total number of gene clusters
                    if lines.split('\t')[1] not in Same_cluster:
                        Same_cluster.append(lines.split('\t')[1])
    else:
        print('file %s is %s' % (same, Checkoutput))
    # calculate MGEs
    Checkoutput = checkfile(mge, 3)
    if Checkoutput == 'not empty':
        Gene_set = []
        for lines in open(mge,'r'):
            Gene1 = lines.split('\t')[2]
            Gene2 = lines.split('\t')[3]
            Gene_pair = function_com(Gene1, Gene2)
            if Gene_pair not in Gene_set:
                Gene_set.append(Gene_pair)
                if "mge_" in Gene1 and "mge_" in Gene2:
                    # mge to mge
                    Result[14] += 1
                else:
                    # mge to genome
                    Result[13] += 1
    else:
        print('file %s is %s' % (mge, Checkoutput))
    # summarize
    Result[3] = len(Diff_cluster)
    Result[4] = '%.3f-%.3f'%((Diff_16S_min),(Cutoff_16S))
    Result[5] = len(Diff_genome_set)
    try:
        # number of hits / total number of combination of random 2 genomes
        Result[7] = "%.3f"%(Result[5]/(Result[6]*(Result[6]-1)/2))
    except ZeroDivisionError:
        Result[7] = 0
    Result[8] = len(Same_cluster)
    Result[10] = len(Same_genome_set)
    try:
        # number of hits / total number of combination of random 2 genomes
        Result[12] = "%.3f"%(Result[10] /(Result[11]*(Result[11]-1)/2))
    except ZeroDivisionError:
        Result[12] = 0
    # output
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
Function_Set_dna=function_load(os.path.join(args.s, args.t + '.all.traits.dna.txt'))
Function_Set_aa=function_load(os.path.join(args.s, args.t + '.all.traits.aa.txt'))
# cluster 16S, by 3% similarity
cluster_16S = run_16S(f16s, Cutoff_16S)
# filter usearch output into same 16S cluster and diff 16S clusters
print('compare %s trait identity to 16S' %('dna'))
script_i = compare_traits_16S(Function_Set_dna, cluster_16S[0],'dna',result_dir +
                              '/sub_sequences',os.path.split(fdna)[-1] + '*.usearch.txt',Cutoff_HGT,script_i)
print('compare %s trait identity to 16S' %('aa'))
script_i = compare_traits_16S(Function_Set_aa, cluster_16S[0],'aa',result_dir +
                             '/sub_sequences',os.path.split(faa)[-1] + '*.usearch.txt',Cutoff_aa,script_i)
print('compare %s trait identity to 16S' %('extended dna'))
script_i = compare_traits_16S(Function_Set_dna, cluster_16S[0],'dna_extended',result_dir +
                             '/sub_sequences',os.path.split(fdna_500)[-1] + '*.usearch.txt',Cutoff_extended,script_i)

# calculate index for HGT
# load each function
all_function_diff_same = glob.glob(os.path.join(result_dir + '/sub_fun','*aa.diff.cluster'))+\
glob.glob(os.path.join(result_dir + '/sub_fun','*aa.same.cluster'))
all_output = open(os.path.join(result_dir,'HGT.summary.dna.%s.aa.%s.16S.%s.txt'
                                   % (Cutoff_HGT,Cutoff_aa,Cutoff_16S)),'w')
all_function = []
all_output.write('function_name\ttype\tcutoff\tcluster_num_diff\t16S_range_diff\thit_diff\tgenome_num_diff\t'+
                      'percentage_diff_pair\tcluster_num_same\t16S_range_same\thit_same\tgenome_num_same\tpercentage_same_pair\t'+
                      'genome_to_mge\tmge_to_mge\n')
for function_diff in all_function_diff_same:
    # for each function
    function_diff = function_diff.replace('.same.cluster', '.diff.cluster')
    function_same = function_diff.replace('.diff.cluster','.same.cluster')
    function_mge = function_diff.replace('.diff.cluster','.mge.cluster')
    function_name = os.path.split(function_diff)[-1].split('.')[0]
    if function_name not in all_function:
        print('summarize potential HGT of %s %s trait with cutoff of %s' % ('dna',function_name,Cutoff_HGT))
        HGT_finder_sum(cluster_16S,function_name,'dna',
                       Cutoff_HGT,
                       function_diff.replace('.aa.diff.cluster', '.dna.diff.cluster'),
                       function_same.replace('.aa.same.cluster', '.dna.same.cluster'),
                       function_diff.replace('.aa.diff.cluster', '.dna.mge.cluster'),
                       all_output)
        print('summarize potential HGT of %s %s trait with cutoff of %s' % ('aa', function_name, Cutoff_aa))
        HGT_finder_sum(cluster_16S, function_name, 'aa',
                       Cutoff_aa,
                       function_diff, function_same, function_mge,
                       all_output)
        print('summarize potential HGT of %s %s trait with cutoff of %s' % ('extended dna', function_name, Cutoff_extended))
        HGT_finder_sum(cluster_16S, function_name, 'dna_extended',
                       Cutoff_extended,
                       function_diff.replace('.aa.diff.cluster', '.dna_extended.diff.cluster'),
                       function_same.replace('.aa.same.cluster', '.dna_extended.same.cluster'),
                       function_diff.replace('.aa.diff.cluster', '.dna_extended.mge.cluster'),
                       all_output)
        all_function.append(function_name)
all_output.close()
# collect all bash files
list_of_files = glob.glob('HGT_subscripts/HGTalign.*.sh')
f1 = open("HGTalign.sh", 'w')
f1.write("#!/bin/bash\nsource ~/.bashrc\n")
for file_name in list_of_files:
    f1.write("jobmit %s HGTalign big\n" % (file_name))
f1.close()