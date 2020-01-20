import os
from Bio import SeqIO
import argparse
import glob
import statistics
import traits_finder
import sys
#from cyvcf2 import VCF


################################################### Decalration #######################################################
print ("\
------------------------------------------------------------------------\n\
traits_finder searches and summarizes traits in genomes and metagenomes\n\
input: reference database and folder of genomes/metagenomes\n\
requirement: blast \n\n\
optional: diamond, bwa, hs-blastn, usearch, mafft, fasttree \n\n\
Copyright:An Ni Zhang, Prof. Eric Alm, MIT\n\n\
Citation:\n\
Contact anniz44@mit.edu\n\
------------------------------------------------------------------------\n\
")

def main():
    usage = ("usage: traits_finder -t your.otu.table -s your.otu.seqs")
    version_string = 'traits_finder {v}, on Python {pyv[0]}.{pyv[1]}.{pyv[2]}'.format(
        v=traits_finder.__version__,
        pyv=sys.version_info,
    )
    ############################################ Arguments and declarations ##############################################
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('command',
                        help="traits_finder genome for analyzing genomes, \
                        traits_finder mge for analyzing mobile genetic elements (MGEs), \
                        traits_finder meta for analyzing metagenomes, \
                        traits_finder sum_genome for summarizing genome results,\
                        traits_finder sum_genome for summarizing MGE results,\
                        traits_finder sum_meta for summarizing metagenome results,\
                             traits_finder HGT for finding candidate HGT of traits across species",
                        type=str,
                        default='genome',
                        choices=['genome','mge', 'meta',
                        'sum_genome','sum_mge','sum_meta',
                        'merge','HGT'],
                        action='store',
                        metavar='traits_finder command')
    required.add_argument("-db",
                        help="file name of your input database",
                        type=str, default='Butyrate.pro.aa',
                        metavar='database.aa')
    optional.add_argument("-dbf",
                        help="sequence format of your input database\
                        (1: nucleotide; 2: protein), \
                        (default \'2\' for protein)",
                        metavar="1 or 2",
                        choices=[1, 2],
                        action='store', default=2, type=int)
    required.add_argument("-m",
                        help="mapping file of traits to function", type=str,
                        default='Butyrate.pro.mapping.txt',
                        metavar='database.mapping.txt')
    required.add_argument("-i",
                        help="input folder of metagenomes or genomes", type=str,
                        default='.',metavar='current dir (.)')
    required.add_argument("-fa",
                        help="input format of genome/metagenome sequence",
                        type=str, default='.fa', metavar='.fasta, .fna, .fastq or .fa')
    optional.add_argument('-s',
                        help="set the method to search the your database \
                        (1: blast; 2: hmm; 3: alignment), \
                        (default \'1\' for blast search)",
                        metavar="1 or 2",
                        choices=[1, 2, 3],
                        action='store', default=1, type=int)
    # optional parameters
    optional.add_argument("--meta",
                        help="metadata  of metagenomes", type=str,
                        default='None',
                        metavar='metadata.metagenomes.txt')
    optional.add_argument("--orf",
                        help="input format of genomes orfs", type=str,
                          default='.genes.faa',metavar='.faa')
    optional.add_argument("--l",
                          help="input list of a subset of metagenomes", type=str,
                          default='None', metavar='list.txt')
    optional.add_argument("--g",
                        help="Optional: gene-level HGT finding; --g T (default: function-level; --g F)",
                        metavar=['T', 'F'], action='store', default='F', type=str)
    # optional output setup
    optional.add_argument("--r",
                        help="output directory or folder of your results",
                        type=str, default='Result',metavar='Result', nargs='*')
    optional.add_argument("--r16",
                        help="output directory or folder of your 16S sequences",
                        type=str, default='Result',metavar='Result')
    # optional search parameters
    optional.add_argument('--t',
                        help="Optional: set the thread number assigned for running XXX (default 1)",
                        metavar="1 or more", action='store', default=1, type=int)
    optional.add_argument('--id','--identity',
                        default=75.0, action='store', type=float, metavar='60.0',
                        help='Optional: set the amno acid based identity cutoff for blast (default is 80.0)\n'
                             'Leave it alone if hmm is used')
    optional.add_argument('--ht','--hitlength',
                        default=75.0, action='store', type=float, metavar='60.0',
                        help='Optional: set the amno acid based hit-length cutoff for blast (default is 80.0)\n'
                             'Leave it alone if hmm is used')
    optional.add_argument('--e','--evalue',
                        default=1e-2, action='store', type=float, metavar='1e-5',
                        help='Optional: set the evalue cutoff for blast or hmm (default is 1e-5)')
    # requirement for software calling
    optional.add_argument('--u','--usearch',
                        help="Optional: use two-step method for blast search,"+
                             " \'None\' for using one step, \'usearch\' for using two-step \
                             (complete path to usearch if not in PATH), (default: \'None\')",
                        metavar="None or usearch",
                        action='store', default='None', type=str)
    optional.add_argument('--dm', '--diamond',
                          help="Optional: use two-step method for blast search," +
                               " \'None\' for using one step, \'diamond\' for using two-step \
                               (complete path to diamond if not in PATH), (default: \'None\')",
                          metavar="None or diamond",
                          action='store', default='None', type=str)
    optional.add_argument('--hs',
                          help="Optional: use two-step method for blast search," +
                               " \'None\' for using one step, \'hs-blastn\' for using two-step \
                               (complete path to hs-blastn if not in PATH), (default: \'None\')",
                          metavar="None or hs-blastn",
                          action='store', default='None', type=str)
    optional.add_argument('--hmm',
                        help="Optional: complete path to hmmscan if not in PATH,",
                        metavar="/usr/local/bin/hmmscan",
                        action='store', default='hmmscan', type=str)
    optional.add_argument('--bp','--blast',
                        help="Optional: complete path to blast if not in PATH, \'None\' for no blast search",
                        metavar="/usr/local/bin/blast",
                        action='store', default='blast', type=str)
    optional.add_argument('--bwa',
                        help="Optional: complete path to bwa if not in PATH,",
                        metavar="/usr/local/bin/bwa",
                        action='store', default='None', type=str)
    optional.add_argument('--mf','--mafft',
                          help="Optional: complete path to mafft if not in PATH,",
                          metavar="/usr/local/bin/mafft",
                          action='store', default='None', type=str)
    optional.add_argument('--ft','--fasttree',
                          help="Optional: complete path to fasttree if not in PATH,",
                          metavar="/usr/local/bin/fasttree",
                          action='store', default='None', type=str)
    ################################################## Definition ########################################################
    args = parser.parse_args()
    workingdir=os.path.abspath(os.path.dirname(__file__))
    ################################################### Function ########################################################
    def split_string_last(input_string, substring):
        return input_string[0: input_string.rfind(substring)]


    def makedatabase(search_method,db_file,db_type):
        # for bwa database
        if 'bwa' in search_method:
            try:
                f1=open(db_file+'.bwt','r')
            except IOError:
                os.system('%s index %s \n' % (search_method, db_file))
        # for blast database
        if db_type == 1:
            #dna database
            if 'blast' in search_method:
                try:
                    f1 = open("%s.nhr" % (db_file), 'r')
                except IOError:
                    os.system('%s -in %s -input_type fasta -dbtype nucl' %
                              (os.path.join(os.path.split(search_method)[0], 'makeblastdb'), db_file))
            if 'hs-blastn' in search_method:
                try:
                    f1 = open("%s.counts.obinary" % (db_file), 'r')
                except IOError:
                    os.system('windowmasker -in %s -infmt blastdb -mk_counts -out %s.counts' %
                              (db_file, db_file))
                    os.system('windowmasker -in %s.counts -sformat obinary -out %s.counts.obinary -convert' %
                              (db_file, db_file))
                    os.system('%s index %s' % (os.path.join(os.path.split(search_method)[0], 'makeblastdb'), db_file))
            if 'hmm' in search_method:
                if '.hmm' not in db_file:
                    try:
                        f1 = open("%s.hmm" % (db_file), 'r')
                    except IOError:
                        os.system('%s %s %s.hmm' %
                              (os.path.join(os.path.split(search_method)[0], 'makehmmerdb'), db_file,db_file))
        else:
            # aa database
            if 'blast' in search_method:
                try:
                    f1 = open("%s.phr" % (db_file), 'r')
                except IOError:
                    os.system('%s -in %s -input_type fasta -dbtype prot' %
                              (os.path.join(os.path.split(search_method)[0], 'makeblastdb'), db_file))
            if 'diamond' in search_method:
                if '.dmnd' not in db_file:
                    try:
                        f1 = open("%s.dmnd" % (db_file), 'r')
                    except IOError:
                        os.system('%sdiamond makedb --in %s -d %s.dmnd' %
                                  (split_string_last(args.dm, 'diamond'), db_file,db_file))
            if 'hmm' in search_method:
                if '.hmm' not in db_file:
                    try:
                        f1 = open("%s.hmm" % (db_file), 'r')
                    except IOError:
                        os.system('%s %s.hmm %s' %
                                  (os.path.join(os.path.split(search_method)[0], 'hmmbuild'), db_file,db_file))
        if 'usearch' in search_method:
            if '.udb' not in db_file:
                try:
                    f1 = open("%s.udb" % (db_file), 'r')
                except IOError:
                    os.system('%s -makeudb_usearch %s -output %s.udb' %
                              (search_method, db_file, db_file))


    ################################################### Programme #######################################################
    f1 = open ('traits_finder.log','w')
    thread = int(args.t)
    # make database
    if args.db == 'but':
        args.db = workingdir + "/database/Butyrate.pro.fasta"
        args.dbf = 2
        args.s = 1
        args.m = workingdir + "/database/Butyrate.pro.fasta.mapping.txt"
    if args.db == 'ARG':
        args.db = workingdir + "/database/SARG.db.fasta"
        args.dbf = 2
        args.s = 1
        args.m = workingdir + "/database/SARG.db.fasta.mapping.txt"
    if args.s == 1:
        if args.bp != 'None':
            makedatabase(args.bp, args.db, args.dbf)
        if args.u != 'None':
            makedatabase(args.u, args.db, args.dbf)
            makedatabase(args.u, workingdir + "/database/85_otus.fasta", 1)
            makedatabase(args.u, workingdir + "/database/85_otus.fasta.all.V4_V5.fasta", 1)
        if args.dm !='None':
            makedatabase(args.dm, args.db, args.dbf)
            makedatabase(args.dm, workingdir + "/database/all_KO30.pro.fasta", 2)
        if args.hs !='None':
            makedatabase(args.hs, args.db, args.dbf)
            makedatabase(args.hs, workingdir + "/database/85_otus.fasta", 1)
            makedatabase(args.hs, workingdir + "/database/85_otus.fasta.all.V4_V5.fasta", 1)
    elif args.s == 2:
        makedatabase(args.hmm, args.db, args.dbf)
    if args.bwa != 'None':
        makedatabase(args.bwa, args.db, args.dbf)


    # run traits finding
    if args.command in ['genome','mge'] :
        cmd = ('python '+workingdir+'/Traits_WG.py -db %s -dbf %s -i %s -s %s --fa %s --orf %s --r %s --r16 %s --t %s --id %s --ht %s --e %s --u %s --dm %s --hs %s --hmm %s --bp %s --bwa %s --mf %s\n'
        % (str(args.db),str(args.dbf),str(args.i),str(args.s),str(args.fa),str(args.orf),str(args.r[0]),str(args.r16),str(thread),
           str(args.id),str(args.ht),str(args.e),str(args.u),str(args.dm),str(args.hs),str(args.hmm),str(args.bp),str(args.bwa),str(args.mf)))
        f1.write(cmd)
        os.system(cmd)
    elif args.command == 'meta':
        cmd = ('python '+workingdir+'/Traits_MG.py -db %s -dbf %s -i %s -s %s --fa %s -l %s --r %s --r16 %s --t %s --id %s --ht %s --e %s --u %s --dm %s --hs %s --hmm %s --bp %s --bwa %s\n'
        % (str(args.db),str(args.dbf),str(args.i),str(args.s),str(args.fa),str(args.l),str(args.r[0]),str(args.r16),str(thread),str(args.id),
           str(args.ht),str(args.e),str(args.u),str(args.dm),str(args.hs),str(args.hmm),str(args.bp),str(args.bwa)))
        f1.write(cmd)
        os.system(cmd)
    elif args.command == 'sum_genome':
        cmd = ('python '+workingdir+'/scripts/Traits_summary_WG.py -t %s -db %s --fa %s --orf %s -i %s -m %s --r %s --r16 %s --s %s -c %s -dbf %s \n'
        %(str(os.path.split(args.db)[1]),str(args.db),str(args.fa),str(args.orf),str(args.i),str(args.m),str(args.r[0]),
          str(args.r16),str(os.path.join(args.r[0],'summary')),str(args.id),str(args.dbf)))
        f1.write(cmd)
        os.system(cmd)
    elif args.command == 'sum_mge':
        cmd = ('python '+workingdir+'/scripts/Traits_summary_WG.py --mge %s -t %s -db %s --fa %s --orf %s -i %s -m %s --r %s --r16 %s --s %s -c %s -dbf %s \n'
        %("2",str(os.path.split(args.db)[1]),str(args.db),str(args.fa),str(args.orf),str(args.i),str(args.m),str(args.r[0]),
          str(args.r16),str(os.path.join(args.r[0],'summary')),str(args.id),str(args.dbf)))
        f1.write(cmd)
        os.system(cmd)
    elif args.command == 'sum_meta':
        cmd = ('python ' + workingdir + '/scripts/Traits_summary_MG.py -db %s -dbf %s -t %s -s %s -fa %s -m %s --r %s --r16 %s --meta %s\n'
                    % (str(args.db), str(args.dbf),str(os.path.split(args.db)[1]), str(args.s), str(args.fa), str(args.m),
                       str(args.r[0]),str(args.r16), str(args.meta)))
        f1.write(cmd)
        os.system(cmd)
    elif args.command == 'merge':
        cmd = ('python '+workingdir+'/scripts/Merge.dataset.py -t %s --r %s \n'
        %(str(os.path.split(args.db)[1]),' '.join(args.r)))
        f1.write(cmd)
        os.system(cmd)
    elif args.command == 'HGT':
        if args.u == 'None' and args.hs == 'None':
            print('please install usearch and/or hs-blastn and/or diamond')
        else:
            cmd = ('python ' + workingdir + '/scripts/HGT_finder_sum.py -db %s -dbf %s -t %s -m %s --r %s --s %s --u %s --dm %s --hs %s --mf %s --ft %s --th %s --bp %s --g %s\n'
            %(str(args.db), str(args.dbf), str(os.path.split(args.db)[1]),str(args.m),str(args.r[0]),
              str(os.path.join(args.r[0],'summary')),str(args.u),str(args.dm),str(args.hs),str(args.mf),str(args.ft),str(thread),str(args.bp),str(args.g)))
            f1.write(cmd)
            os.system(cmd)

################################################## Function ########################################################

if __name__ == '__main__':
    main()
