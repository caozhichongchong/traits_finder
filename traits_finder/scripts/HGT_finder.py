import os
from Bio import SeqIO
import argparse
import glob

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-t",
                    help="trait name", type=str, default='ARG',metavar='trait name')
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
    os.mkdir(result_dir)
except OSError:
    pass
try:
    os.mkdir(result_dir + '/sub_sequences')
except OSError:
    pass
workingdir=os.path.abspath(os.path.dirname(__file__))
try:
    os.mkdir('HGT_subscripts')
except OSError:
    pass

################################################### Function ########################################################
def checkfile(filename, i):
    try:
        f1 = open(filename, 'r')
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


def loci_seq(record_name):
    try:
        loci1 = int(record_name.split('_')[-2])
        loci2 = int(record_name.split('_')[-1])
        return [loci1,loci2]
    except ValueError:
        print(record_name)

def compare_loci(loci_new,loci_ref):
    if loci_new[0] <= loci_ref[0] and loci_new[1] >= loci_ref[1]:
        return 1
    else:
        return 0

def cutoff_set(filesize,cutoff):
    if int(filesize) <= 10*1E+6:
        # 10Mb
        return cutoff - 0.05
    elif int(filesize) <= 100*1E+6:
        # 0.1G
        return cutoff - 0.1
    elif int(filesize) <= 1*1E+9:
        # 1G
        return cutoff - 0.2
    else:
        # larger than 1G, run self-clustering
        return 0


def self_clustering(input_usearch, output_uc, cutoff):
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
                if Gene1 not in clusters and Gene2 not in clusters:
                    cluster_num += 1
                    if ID < cutoff:
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
                    if ID < cutoff:
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
        print('file %s is %s' % (input_usearch, Checkoutput))
    output_uc.close()


def run_cluster(input_fasta, output_clusters = 1, cutoff=0.99):
    # run cutoff_set
    try:
        f1 = open("%s.uc" % (input_fasta),'r')
    except IOError:
        filesize = os.path.getsize(input_fasta)
        newcutoff = cutoff_set(filesize, cutoff)
        if newcutoff > 0:
            print('Running usearch cluster for %s, filesize %s, cutoff %s' % (input_fasta,filesize,newcutoff))
            os.system('%s -sort length -cluster_fast %s -id %s -centroids %s.cluster.aa -uc %s.uc -threads %s'
                      %(args.u, input_fasta,str(newcutoff),input_fasta,input_fasta,args.th))
        else:
            print ('Please install hs-blastn because the input file is larger than 2GB\n')
            print ('Running self-clustering using hs-blastn\n')
            newcutoff = cutoff - 0.2
            print('Running self-clustering for %s, filesize %s, cutoff %s' % (input_fasta, filesize, newcutoff))
            try:
                f1 = open("%s.%s.usearch.txt" % (input_fasta,newcutoff), 'r')
            except IOError:
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
                       input_fasta, input_fasta, input_fasta, newcutoff,
                       newcutoff, str(args.th)))
            self_clustering('%s.%s.usearch.txt' % (input_fasta, newcutoff), input_fasta + '.uc', newcutoff)
    print('Finish running usearch cluster for %s' %
          (input_fasta))
    # read cluster results
    Clusters = dict()
    Clusters_seqs = dict()
    Clusters_extend_seqs = dict()
    Total_seq = 0
    Cluster_num = 0
    Checkoutput = checkfile(input_fasta+'.uc', 8)
    if Checkoutput == 'not empty':
        for lines in open(input_fasta+'.uc','r'):
                cluster = lines.split('\t')[1]
                record_name = lines.split('\t')[8].split(' ')[0]
                if cluster not in Clusters:
                    Clusters.setdefault(cluster,[record_name])
                else:
                    Clusters[cluster].append(record_name)
                Total_seq += 1
                Cluster_num = max(Cluster_num,int(cluster))
    else:
        print('file %s is %s' % (input_fasta+'.uc', Checkoutput))
    if Cluster_num > 0:
        cutoff_temp = max(3, int(Total_seq / Cluster_num))
        print('A total of %s clusters for %s sequences\nSetting cutoff as %s'
              % (Cluster_num, Total_seq, cutoff_temp))
        for cluster in Clusters:
            # at least 3 sequences in a cluster
            if len(Clusters[cluster]) >= cutoff_temp:
                for record_name in Clusters[cluster]:
                    Clusters_seqs.setdefault(record_name, str(cluster))
                    record_subname = '_'.join(record_name.split('_')[0:-2])
                    if output_clusters == 1:
                        if record_subname not in Clusters_extend_seqs:
                            Clusters_extend_seqs.setdefault(record_subname, [[cluster, record_name,
                                                                              loci_seq(record_name)]])
                        else:
                            Clusters_extend_seqs[record_subname].append([cluster, record_name,
                                                                         loci_seq(record_name)])
        # clustering sequences
        print('Clustering sequences for %s' % (input_fasta))
        for record in SeqIO.parse(input_fasta, 'fasta'):
            if str(record.id) in Clusters_seqs:
                output_sequences = open(os.path.join(result_dir + '/sub_sequences', '%s.%s.fasta' % (
                    os.path.split(input_fasta)[-1],
                    Clusters_seqs[str(record.id)])), 'a')
                output_sequences.write('>%s\n%s\n' % (str(record.id), str(record.seq)))
                output_sequences.close()
    return Clusters_extend_seqs


def extract_extend(Clusters_extend_seqs,input_extend_fasta):
    # check redundancy!!!
    print('Extracting extended sequences from %s' % (input_extend_fasta))
    output_map_file = open(os.path.join(result_dir, args.t + '.all.traits.dna-dna_extended.mapping.txt'),
                       'w')
    for record in SeqIO.parse(input_extend_fasta, 'fasta'):
        record_subname = '_'.join(str(record.id).split('_')[0:-2])
        if record_subname in Clusters_extend_seqs:
            for clusters in Clusters_extend_seqs[record_subname]:
                if compare_loci(loci_seq(str(record.id)),clusters[-1]):
                    # the right extended sequences
                    output_map_file.write('%s\t%s\t%s\t%s\n'%(record_subname,str(clusters[0]),
                                                      clusters[-2],str(record.id)))
                    output_sequences = open(os.path.join(result_dir + '/sub_sequences', '%s.%s.fasta' % (
                        os.path.split(input_extend_fasta)[-1],
                        clusters[0])), 'a')
                    output_sequences.write('>%s\n%s\n' % (str(record.id), str(record.seq)))
                    output_sequences.close()
                    break
    output_map_file.close()


def run_alignment(input_folder,type_fasta,cutoff,script_i):
    print('Running alignment for %s' % (input_folder))
    all_clusters = glob.glob(os.path.join(result_dir + '/sub_sequences', '%s.*.fasta' % (
        os.path.split(input_folder)[-1])))
    for input_fasta in all_clusters:
        if "usearch" in args.u:
            try:
                f1 = open("%s.sorted" % (input_fasta), 'r')
            except IOError:
                os.system('%s -sortbylength %s -fastaout %s.sorted' % (args.u, input_fasta,input_fasta))
            input_fasta=input_fasta + '.sorted'
            try:
                f1 = open("%s.%s.usearch.txt" % (input_fasta,str(cutoff)), 'r')
            except IOError:
                output_file = open(('HGT_subscripts/HGTfinder.%s.sh'%(int(script_i % script_i_max))), 'a')
                script_i += 1
                output_file.write("#!/bin/bash\n")
                output_file.write("%s -makeudb_usearch %s -output %s.udb\n"
                                  % (args.u, input_fasta,input_fasta))
                if  'dna' in type_fasta :
                    output_file.write("%s -usearch_global %s -db %s.udb  -strand both -id %s -maxaccepts 0 -maxrejects 32 -blast6out %s.%s.usearch.txt  -threads %s\n"
                                      %(args.u,input_fasta,input_fasta,str(cutoff),input_fasta,str(cutoff),str(args.th)))
                else:
                    output_file.write(
                        "%s  -usearch_global %s -db %s.udb  -id %s  -evalue 1e-2 -maxaccepts 0 -maxrejects 32 -blast6out %s.%s.usearch.txt  -threads %s\n"
                        % (args.u, input_fasta, input_fasta,str(cutoff),input_fasta,str(cutoff), str(args.th)))
                output_file.close()
    return script_i


################################################### Programme #######################################################
faa = os.path.join(args.s, args.t + '.all.traits.aa.fasta')
fdna = os.path.join(args.s, args.t + '.all.traits.dna.fasta')
fdna_500 = glob.glob(os.path.join(args.s, args.t + '.all.traits.dna.extra*.fasta'))[0]

# cluster dna sequences
# identity cutoff is 99%, but run cutoff - x
Clusters_extend_seqs = run_cluster(fdna,1,Cutoff_HGT)
# extract 500bp extended sequences for dna clusters and cluster extended sequences
if Clusters_extend_seqs != dict():
    extract_extend(Clusters_extend_seqs, fdna_500)
    all_extend_clusters = glob.glob(os.path.join(result_dir + '/sub_sequences', '%s.*.fasta' % (
            os.path.split(fdna_500)[-1])))
# cluster aa sequences
# identity cutoff is 90%
run_cluster(faa,0,Cutoff_aa)
# run alignment for dna
script_i = run_alignment(fdna,'dna',Cutoff_HGT, script_i)
script_i = run_alignment(fdna_500,'dna',Cutoff_extended,script_i)
# run alignment for aa
script_i = run_alignment(faa,'aa',Cutoff_aa,script_i)
# collect all bash files
list_of_files = glob.glob('HGT_subscripts/HGTfinder.*.sh')
f1 = open("HGTfinder.sh", 'w')
f1.write("#!/bin/bash\nsource ~/.bashrc\n")
for file_name in list_of_files:
    f1.write("jobmit %s HGTfinder big\n" % (file_name))
f1.close()