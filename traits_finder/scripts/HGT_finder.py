import os
from Bio import SeqIO
import argparse
import glob

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-t",
                    help="trait name", type=str, default='ARG',metavar='trait name')
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
workingdir=os.path.abspath(os.path.dirname(__file__))

try:
    os.mkdir('HGTfinder_subscripts')
except OSError:
    pass

################################################### Function ########################################################
def loci_seq(record_name):
    loci1 = int(record_name.split('_')[-2])
    loci2 = int(record_name.split('_')[-1])
    return [loci1,loci2]

def compare_loci(loci_new,loci_ref):
    if loci_new[0] <= loci_ref[0] and loci_new[1] >= loci_ref[1]:
        return 1
    else:
        return 0


def run_cluster(input_fasta,cutoff=0.99):
    print('Running usearch cluster for %s' % (input_fasta))
    os.system('usearch -sort length -cluster_fast %s -id %s -centroids %s.cluster.aa -uc %s.uc'
              %(input_fasta,cutoff,input_fasta,input_fasta))
    print('finish running usearch cluster for %s' % (input_fasta))
    # read cluster results
    Clusters = dict()
    Clusters_seqs = dict()
    Clusters_extend_seqs = dict()
    for lines in open(input_fasta+'.rc','r'):
        cluster = lines.split('\t')[1]
        record_name = lines.split('\t')[-2].split(' ')[0]
        if cluster not in Clusters:
            Clusters.setdefault(cluster,[record_name])
        else:
            Clusters[cluster].append(record_name)
    for cluster in Clusters:
        # at least 3 sequences in a cluster
        if len(Clusters[cluster] > 2):
            for record_name in Clusters[cluster]:
                Clusters_seqs.setdefault(record_name,str(cluster))
                record_subname = '_'.join(record_name.split('_')[0:-2])
                if record_subname not in Clusters_extend_seqs:
                    Clusters_extend_seqs.setdefault(record_subname,[[Clusters_seqs,loci_seq(record_name)]])
                else:
                    Clusters_extend_seqs[record_subname].append([Clusters_seqs,loci_seq(record_name)])
    # clustering sequences
    for record in SeqIO.parse(input_fasta, 'fasta'):
        if str(record.id) in Clusters_seqs:
            output_sequences=open(os.path.join(result_dir,'%s.%s.fasta'%(
                os.path.split(input_fasta)[-1],
                Clusters_seqs[str(record.id)])),'a')
            output_sequences.write('>%s\n%s\n'%(str(record.id),str(record.seq)))
            output_sequences.close()
    return Clusters_extend_seqs


def extract_extend(Clusters_extend_seqs,input_extend_fasta):
    for record in SeqIO.parse(input_extend_fasta, 'fasta'):
        record_subname = '_'.join(str(record.id).split('_')[0:-2])
        if record_subname in Clusters_extend_seqs:
            for clusters in Clusters_extend_seqs[record_subname]:
                if compare_loci(loci_seq(str(record.id)),clusters[-1]):
                    # the right extended sequences
                    output_sequences = open(os.path.join(result_dir, '%s.%s.fasta' % (
                        os.path.split(input_extend_fasta)[-1],
                        clusters[0])), 'a')
                    output_sequences.write('>%s\n%s\n' % (str(record.id), str(record.seq)))
                    output_sequences.close()


def run_alignment(input_folder,type_fasta,output_file):
    all_clusters = glob.glob(os.path.join(result_dir, '%s.*.fasta' % (
        os.path.split(input_folder)[-1])))
    for input_fasta in all_clusters:
        if "usearch" in args.u:
            os.system('usearch -sortbylength %s -fastaout %s.sorted' % (input_fasta,input_fasta))
            input_fasta=input_fasta + '.sorted'
        if args.mf != 'None' and args.ft != 'None':
            if type_fasta == 'dna':
                output_file.write("%s --nuc --adjustdirection --quiet --nofft --maxiterate 0 --retree 1 --thread %s %s > %s.align\n"
                 % (args.mf, str(args.th), input_fasta, input_fasta))
                output_file.write("%s -nt -quiet -fastest -nosupport %s.align > %s.align.nwk\n"
                 % (args.ft, input_fasta, input_fasta))
            else:
                output_file.write(
                    "%s --amino --quiet --retree 1 --maxiterate 0 --nofft --thread %s %s > %s.align\n"
                    % (args.mf, str(args.th), input_fasta, input_fasta))
                output_file.write("%s -quiet -fastest -nosupport %s.align > %s.align.nwk\n"
                                  % (args.ft, input_fasta, input_fasta))


################################################### Programme #######################################################
faa = os.path.join(args.s, args.t + '.all.traits.aa.*.fasta')
fdna = os.path.join(args.s, args.t + '.all.traits.dna.*.fasta')
fdna_500 = os.path.join(args.s, args.t + '.all.traits.dna.extra500.fasta')
output_file = open(os.path.join('HGTfinder_subscripts/HGTfinder_subscripts.sh'), 'w')
output_file.write("#!/bin/bash\n")
# cluster dna sequences
# identity cutoff is 99%
Clusters_extend_seqs = run_cluster(fdna)
# extract 500bp extended sequences for dna clusters and cluster extended sequences
# identity cutoff is 60%
extract_extend(Clusters_extend_seqs, fdna_500)
all_extend_clusters = glob.glob(os.path.join(result_dir, '%s.*.fasta' % (
        os.path.split(fdna_500)[-1])))
for input_fasta in all_extend_clusters:
    run_cluster(input_fasta,0.6)
    # run alignment for dna
    run_alignment(input_fasta, 'dna', output_file)
# cluster aa sequences
# identity cutoff is 99%
run_cluster(faa)
# run alignment for dna
run_alignment(fdna,'dna',output_file)
# run alignment for aa
run_alignment(faa,'aa',output_file)
output_file.close()
