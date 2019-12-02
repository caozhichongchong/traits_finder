################################################### Function ########################################################
def run_16S(input_fasta,cutoff=0.971):
    print('Running usearch cluster for %s' % (input_fasta))
    os.system('usearch -sort length -cluster_fast %s -id %s -centroids %s.cluster.aa -uc %s.uc'
              % (input_fasta, cutoff, input_fasta, input_fasta))
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
    return Clusters_seqs


def function_check():
    # no reference
    # aa do alignment
    # dna do extract alignment

def sum_alignment():


def compare_traits_16S():

################################################### Programme #######################################################
f16s = os.path.join(args.s, args.t + '.all.16S.fasta')
# cluster 16S
# by 3% dissimilarity
cluster_16S = run_16S(f16s, 0.971)
# compare usearch output of dna directly
# compare usearch output of extracted 500bp directly