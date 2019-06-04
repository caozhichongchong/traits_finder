import os
from Bio import SeqIO
import argparse
import glob


############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-db",
                    help="file name of your input database", type=str, default='Butyrate.pro.aa',metavar='database.aa')
parser.add_argument("-i",
                    help="input dir of WGD", type=str, default='.',metavar='current dir (.)')
parser.add_argument('-s',
                    help="set the method to search the your database \
                    (1: blast; 2: hmm), \
                    (default \'1\' for blast search)",
                    metavar="1 or 2",
                    choices=[1, 2],
                    action='store', default=1, type=int)
# optional parameters
parser.add_argument("--fa",
                    help="input format of genome sequence", type=str, default='.fna.add',metavar='.fasta, .fna or .fa')
parser.add_argument("--orf",
                    help="input format of genome orfs", type=str, default='.genes.faa',metavar='.faa')
# optional output setup
parser.add_argument("--r",
                    help="output directory or folder of your results",
                    type=str, default='Result',metavar='Result')
parser.add_argument("--r16",
                    help="output directory or folder of your 16S sequences",
                    type=str, default='Result',metavar='Result')
# optional search parameters
parser.add_argument('--t',
                    help="Optional: set the thread number assigned for running (default 1)",
                    metavar="1 or more", action='store', default=1, type=int)
parser.add_argument('--id',
                    default=50.0, action='store', type=float, metavar='80.0',
                    help='Optional: set the amno acid based identity cutoff for blast (default is 80.0)\n'
                         'Leave it alone if hmm is used')
parser.add_argument('--ht',
                    default=50.0, action='store', type=float, metavar='80.0',
                    help='Optional: set the amno acid based hit-length cutoff for blast (default is 80.0)\n'
                         'Leave it alone if hmm is used')
parser.add_argument('--e',
                    default=1e-5, action='store', type=float, metavar='1e-5',
                    help='Optional: set the evalue cutoff for blast or hmm (default is 1e-5)')
# requirement for software calling
parser.add_argument('--u',
                    help="Optional: use two-step method for blast search,"+
                         " \'None\' for using one step, \'usearch\' or \'diamond\' for using two-step \
                         (complete path to usearch or diamond if not in PATH, \
                         please make sure the search tools can be directly called), (default: \'None\')",
                    metavar="None or usearch",
                    action='store', default='None', type=str)
parser.add_argument('--hmm',
                    help="Optional: complete path to hmmscan if not in PATH,",
                    metavar="/usr/local/bin/hmmscan",
                    action='store', default='hmmscan', type=str)
parser.add_argument('--bp',
                    help="Optional: complete path to blastp if not in PATH,",
                    metavar="/usr/local/bin/blastp",
                    action='store', default='blastp', type=str)

################################################## Definition ########################################################
args = parser.parse_args()
fasta_format = args.fa
orfs_format = args.orf
try:
    os.mkdir(args.r)
except OSError:
    pass


################################################### Function ########################################################
def addname(filedir, file_name):
    Fasta_name = open(os.path.join(filedir,file_name), 'r')
    f = open(os.path.join(filedir,file_name + '.add'), 'w')
    in_dir, input_file = os.path.split(file_name)
    for record in SeqIO.parse(Fasta_name, 'fasta'):
        if len(str(record.seq).replace(' ',''))>0:
            # remove empty ORF sequences, otherwise there could be a problem for blast, usearch and diamond
            f.write('>'+str(input_file)+'_'+str(record.id) + '\t' + str(record.description).replace('\t', ' ') + '\n' + str(
                str(record.seq)) + '\n')
    f.close()


def search(roottemp,filename):
    # search the database for each file
    cmds = ''
    cmds += "#!/bin/bash \nmodule add c3ddb/blast+/2.7.1 \nsource activate py37\n"
    # using blastp
    if args.s == 1:
        if args.u != 'None':
            # two-step search
            Usearch=0
            for root, dirs, files in os.walk(args.r + '/usearch'):
                try:
                    ftry = open(os.path.join(root, filename + '.usearch.txt.aa'), 'r')
                    searchfile = os.path.join(root, filename + '.usearch.txt.aa')
                    Usearch=1
                except IOError:
                    pass
            if Usearch == 0:
                if 'usearch' in args.u:
                    # Start search target genes by usearch
                    cmds += args.u + " -ublast " + os.path.join(roottemp, filename) + \
                            " -db " + args.db + ".udb -evalue 1e-2 -accel 0.5 -blast6out " \
                            + os.path.join(args.r + '/usearch/' + str(int(i/10000)), filename + '.usearch.txt') + \
                            " -threads " + str(int(i_max)) + " \n"
                elif "diamond" in args.u:
                    # Start search target genes by diamond!
                    cmds += args.u + " blastp --query " + os.path.join(roottemp, filename) + \
                            " --db " + args.db + ".dmnd --out " + os.path.join(args.r + '/usearch/' + str(int(i/10000)),
                                                                          filename + '.usearch.txt') + \
                            " --outfmt 6 --max-target-seqs 1 --evalue " + str(args.e) + " --threads " + str(
                        int(i_max)) + " \n"
                cmds += 'python scripts/Extract.py -i ' + roottemp + ' -f ' + filename + ' -n .usearch.txt -r ' + args.r + '/usearch/' + str(
                    int(i / 10000)) + ' \n'
                searchfile = os.path.join(args.r + '/usearch/' + str(int(i/10000)), filename + '.usearch.txt.aa')
        else:
            # one-step search
            searchfile = os.path.join(roottemp, filename)
        # blast search
        Blastsearch = 0
        for root, dirs, files in os.walk(args.r + '/search_output'):
            try:
                ftry = open(os.path.join(root, filename + '.blast.txt'), 'r')
                Blastsearch = 1
            except IOError:
                pass
        if Blastsearch == 0:
            cmds += str(args.bp) +" -query " + str(searchfile) + " -db " + args.db + " -out " + args.r + '/search_output/'+str(int(i/10000))+ \
                     "/"+filename+".blast.txt  -outfmt 6 -max_target_seqs 1 -evalue "+str(args.e)+" -num_threads " + \
                    str(int(i_max)) + " \n"
            # fiter blast result
            cmds += 'python scripts/Filter.py -i ' + args.r + '/search_output/'+str(int(i/10000))+' -f ' + filename + '.blast.txt ' +\
                    '-db ' + args.db + ' -s ' + str(args.s) + ' --ht ' + str(args.ht) + ' --id ' + str(args.id) + \
                    ' --e ' + str(args.e) + ' \n'
            cmds += 'python scripts/Extract.py -i ' + roottemp + ' -f ' + filename + ' -n .blast.txt.filter -r ' + args.r + '/search_output/'+str(int(i/10000))+' \n'
    else:
        # hmmsearch
        Blastsearch = 0
        for root, dirs, files in os.walk(args.r + '/search_output'):
            try:
                ftry = open(os.path.join(root, filename + '.hmm'), 'r')
                Blastsearch = 1
            except IOError:
                pass
        if Blastsearch == 0:
            cmds = args.hmm + ' --tblout ' + os.path.join(args.r + '/search_output/'+str(int(i/10000)), str(
                filename) + '.hmm') +  ' --cpu ' + str(int(i_max)) + ' -E ' \
                  + str(args.e) + ' ' +args.db + ' '+ os.path.join(roottemp, filename) + ' \n'
            cmds += 'python scripts/Format.py -i ' + args.r + '/search_output/'+str(int(i/10000)) + ' -f ' + str(
                    filename) + '.hmm \n'
    # 16S extraction
    Search16s = 0
    for root, dirs, files in os.walk(args.r16):
        try:
            ftry = open(os.path.join(root, filename.replace(orfs_format+'.add', fasta_format) + '.16S.txt'), 'r')
            Search16s = 1
        except IOError:
            pass
    if Search16s == 0:
        pass
        # Start search 16S by usearch
        #cmds += 'python scripts/undone.py -i '+ os.path.join(roottemp.replace('_faa','_fasta'), filename.replace(orfs_format+'.add', fasta_format)+' \n')
        # with usearch
        #cmds += args.u + " -usearch_global " + os.path.join(roottemp.replace('_faa','_fasta'), filename.replace(orfs_format+'.add', fasta_format))+ \
        #        " -db database/85_otus.fasta.udb -strand plus -id 0.7 -evalue 1e-1 -blast6out " \
        #        + os.path.join(args.r16+'/' + str(int(i/10000)), filename.replace(orfs_format+'.add', fasta_format) + '.16S.txt')  + \
        #        " -threads " + str(int(i_max)) + " \n"
        #cmds += str(args.bp).replace('blastp','blastn') +" -query " + os.path.join(roottemp.replace('_faa','_fasta'), filename.replace(orfs_format+'.add', fasta_format))\
        #        + " -db database/85_otus.fasta -out " + os.path.join(args.r16+'/' + str(int(i/10000)), filename.replace(orfs_format+'.add', fasta_format) + '.16S.txt') +\
        #"  -outfmt 6  -max_target_seqs 1 -evalue 1e-5 -num_threads " + \
        #            str(int(i_max)) + " \n"
        #cmds += 'python scripts/Extract.16S.py -i ' + roottemp.replace('_faa','_fasta') + ' -f ' + \
        #        filename.replace(orfs_format+'.add',fasta_format) + ' -n .16S.txt -r ' + args.r16 + '/' + str(
        #    int(i / 10000)) + ' \n'
    return cmds


################################################### Programme #######################################################
# load all WGD
flist=open('Filelist.txt','w')
in_dir=args.i
Targetroot=dict()
for root,dirs,files in os.walk(in_dir):
    list_fasta1 = glob.glob(os.path.join(root, '*'+orfs_format))
    if list_fasta1!=[]:
        for files in list_fasta1:
            Targetroot.setdefault(files, orfs_format)


# search the database in all genomes
i=0
#os.system("rm -rf *.sh \n")
#os.system("rm -rf nohup.sh \n")
#i_max=max(int(args.t)/len(Targetroot),1)
i_max = 40
for files in Targetroot:
    if Targetroot[files]!='None':
        f1 = open(str(i%int(args.t)) + '.sh', 'a')
        i +=1
        roottemp, filename = os.path.split(files)
        flist.write(str(files) + '\n')
        try:
            os.mkdir(args.r + '/search_output')
        except OSError:
            pass
        try:
            os.mkdir(args.r + '/usearch')
        except OSError:
            pass
        try:
            os.mkdir(args.r16)
        except OSError:
            pass
        try:
            os.mkdir(args.r + '/search_output/'+str(int(i/10000)))
        except OSError:
            pass
        try:
            os.mkdir(args.r + '/usearch/'+str(int(i/10000)))
        except OSError:
            pass
        try:
            os.mkdir(args.r16+'/'+str(int(i/10000)))
        except OSError:
            pass
        # add filename to orf
        try:
            ftry = open(os.path.join(roottemp, str(filename) + ".add"), 'r')
        except IOError:
            pass
            #addname(roottemp, str(filename))
        #filename = filename + ".add"
        # search the database in WGD
        cmds = search(roottemp, filename)
        f1.write(cmds.replace('IMG Data','IMG\ Data'))
        f1.close()
flist.close()
