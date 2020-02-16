import os
from Bio import SeqIO
import argparse
import glob
from datetime import datetime
import statistics
import random
import subprocess


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
                        type=str, default='_final.scaffolds.fasta'
                    , metavar='.fasta, .fna, .fastq or .fa')
# optional input setup
parser.add_argument("--r",
                    help="input directory or folder of your previous results by Traits_WGD.py",
                    type=str, default='Result',metavar='Result')
parser.add_argument("--s",
                    help="input directory or folder of your previous results by traits summary",
                    type=str, default='None',metavar='summary')
# optional search parameters
parser.add_argument('--g',
                        help="Optional: gene-level HGT finding; --g F (default function-level); --g T for gene-level",
                        metavar=['T', 'F'], action='store', default='F', type=str)
parser.add_argument('--th',
                        help="Optional: set the thread number assigned for running XXX (default 1)",
                        metavar="1 or more", action='store', default=5, type=int)
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
                      action='store', default='fastANI', type=str)

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
                except IndexError:
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
    os.system('rm -rf Filelist.temp.*.output')
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
    os.system('cat Filelist.temp.*.output > Filelist.temp.all.output')
    os.system('sed -i -e \"s/%s//g\"  -e \"s/%s//g\" Filelist.temp.all.output' % (args.fa, '-'))
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
                            except IndexError:
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
                list_fasta1 = glob.glob(
                    os.path.join('/scratch/users/anniz44/genomes/GHM/newgenomes/searched/',
                                 '*_final.scaffolds.fasta')) + \
                              glob.glob(os.path.join('/scratch/users/anniz44/genomes/GHM/newgenomes/filtered/',
                                                     '*_final.scaffolds.fasta'))
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

run_fastani('GMCBN10.fastani.out','.','_final.scaffolds.fasta')