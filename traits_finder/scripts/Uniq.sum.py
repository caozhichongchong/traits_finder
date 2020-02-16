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
# optional input setup
parser.add_argument("--taxa",
                    help="mapping file of taxa", type=str,
                    default='/scratch/users/anniz44/scripts/maffttree/GTDB_taxon_CG_GMC.brief.txt',
                    metavar='GTDB_taxon_CG_GMC.brief.txt')
parser.add_argument("--s",
                    help="input directory or folder of your previous results by traits summary",
                    type=str, default='/scratch/users/anniz44/genomes/plasmid/ARG_gene_subset/merge/summary',
                    metavar='summary')
parser.add_argument('--donor',
                        help="Optional: donor_level count",
                        metavar=['T', 'F'], action='store', default='F', type=str)


################################################## Definition ########################################################
args = parser.parse_args()

################################################### new class #########################################################
__metaclass__ = type

class Uniq_gene:
    # create a class to store HGT_function
    'a class to store HGT_function'

    def init(self, gene, function):
        self.function = function
        self.name = gene
        self.species = set()
        self.genome_set = set()
        self.mge = set()
        self.geneset = set()

    def addspecies(self, species):
        self.species.add(str(species))

    def addmge(self, mge):
        self.mge.add(str(mge))

    def addgenome(self, genome):
        self.genome_set.add(genome)

    def addgene(self, gene):
        self.geneset.add(str(gene))


################################################### Function #######################################################
def taxonomy_read(input,column_num1,column_num2):
    print('%s loading taxonomy file %s' % (datetime.now(),input))
    Taxonomy_set = dict()
    for lines in open(input,'r'):
        taxonomy = str(lines).split('\t')[column_num2].replace('\r','').replace('\n','')
        if taxonomy == 'NA':
            taxonomy = 'Other'
        Taxonomy_set.setdefault(str(lines).split('\t')[column_num1],taxonomy)
    return Taxonomy_set

def split_string_last(input_string,substring):
    last_loci = input_string.rfind(substring)
    if last_loci > -1:
        return input_string[0 : last_loci]
    else:
        return input_string

def mapping(recordID):
    if recordID.startswith('mge'):
        return ['_'.join(recordID.split('_', maxsplit=3)[0:2]).split('.')[0],'mge','None']
    elif recordID.startswith('GCA'):
        genomename = '_'.join(recordID.split('_', maxsplit=3)[0:2])
        donor = 'None'
    else:
        genomename = '_'.join(recordID.split('_', maxsplit=5)[0:4])
        donor = genomename.split('_', maxsplit=2)[0]
    return [genomename,Taxonomy_Set.get(genomename,'None'),donor]


def count_uniq(list_file):
    all_output = open(list_file + '.count', 'w')
    all_output.write('function_name\tgene_name\tuniq_gene_num\tspecie_num\tgenome_num\twithmge\tuniq_gene_list\tspecie_list\n')
    Function_count = dict()
    print('%s process unique list file %s' % (datetime.now(), list_file))
    line_num = 0
    for lines in open(list_file):
        line_set = lines.split('\t')
        gene_name = line_set[0]
        function = split_string_last(gene_name,'-')
        Genome = line_set[1]
        genomename, taxonomy, donor = mapping(Genome)
        line_num += 1
        if line_num%1000 == 0:
            print('%s process %s genes' % (datetime.now(), line_num))
        if function not in Function_count:
            Uniq_function_temp = Uniq_gene()
            Uniq_function_temp.init(function, function)
            Function_count.setdefault(function,Uniq_function_temp)
        if gene_name not in Function_count:
            Uniq_gene_temp = Uniq_gene()
            Uniq_gene_temp.init(gene_name,function)
            Function_count.setdefault(gene_name, Uniq_gene_temp)
        if args.donor == 'T':
            gene_name_donor = '%s:%s'%(gene_name,donor)
            if gene_name_donor not in Function_count:
                Uniq_gene_donor_temp = Uniq_gene()
                Uniq_gene_donor_temp.init(gene_name_donor, function)
                Function_count.setdefault(gene_name_donor, Uniq_gene_donor_temp)
            function_donor = '%s:%s' % (function, donor)
            if function_donor not in Function_count:
                Uniq_function_donor_temp = Uniq_gene()
                Uniq_function_donor_temp.init(function_donor, function)
                Function_count.setdefault(function_donor, Uniq_function_donor_temp)
            Uniq_gene_donor_temp = Function_count[gene_name_donor]
            Uniq_function_donor_temp = Function_count[function_donor]
        else:
            Uniq_gene_donor_temp = Function_count[gene_name]
            Uniq_function_donor_temp = Function_count[function]
        Uniq_function_temp = Function_count[function]
        Uniq_function_temp.addgene(gene_name)
        Uniq_function_donor_temp.addgene(gene_name)
        Uniq_gene_temp = Function_count[gene_name]
        Uniq_gene_temp.addgene(gene_name)
        Uniq_gene_donor_temp.addgene(gene_name)
        if taxonomy == 'mge':
            Uniq_function_temp.addmge(genomename)
            Uniq_function_donor_temp.addmge(genomename)
            Uniq_gene_temp.addmge(genomename)
            Uniq_gene_donor_temp.addmge(genomename)
        else:
            Uniq_function_temp.addgenome(genomename)
            Uniq_function_donor_temp.addgenome(genomename)
            Uniq_gene_temp.addgenome(genomename)
            Uniq_gene_donor_temp.addgenome(genomename)
            if taxonomy != 'None':
                Uniq_function_temp.addspecies(taxonomy)
                Uniq_function_donor_temp.addspecies(taxonomy)
                Uniq_gene_temp.addspecies(taxonomy)
                Uniq_gene_donor_temp.addspecies(taxonomy)
    Result_list = []
    print('%s output unique list summary %s' % (datetime.now(), list_file + '.count'))
    for gene_name in Function_count:
        Uniq_gene_temp = Function_count[gene_name]
        Result_list.append('\t'.join([
            Uniq_gene_temp.function,Uniq_gene_temp.name,
            str(len(Uniq_gene_temp.geneset)),
            str(len(Uniq_gene_temp.species)),
            str(len(Uniq_gene_temp.genome_set)),
            ';'.join(Uniq_gene_temp.mge),
            ';'.join(Uniq_gene_temp.geneset),
            ';'.join(Uniq_gene_temp.species)
        ]))
    all_output.write('\n'.join(Result_list)+'\n')
    all_output.close()



################################################### Programme #######################################################
Taxonomy_Set = taxonomy_read(args.taxa,0,-1)
faa = glob.glob(os.path.join(args.s,'*.aa.fasta.unique_list'))[0]
fdna = glob.glob(os.path.join(args.s,'*.dna.fasta.unique_list'))[0]
count_uniq(fdna)
count_uniq(faa)
