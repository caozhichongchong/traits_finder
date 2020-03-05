import os
import os
from Bio import SeqIO
import argparse
import glob
#from cyvcf2 import VCF


############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
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
                    help="input dir of MG", type=str, default='.',metavar='current dir (.)')
parser.add_argument("-l",
                    help="input list of MG", type=str, default='None',metavar='rep_metagenomes.txt')
parser.add_argument('-m',
                    help="set the model to predict the results \
                    (1: simple; 2: phylogenetic; 3: advanced), \
                    (default \'1\' for simple)",
                    metavar="1 or 2 or 3",
                    choices=[1, 2, 3],
                    action='store', default=1, type=int)
parser.add_argument('-s',
                    help="set the method to search the your database \
                    (1: blast; 2: hmm; 3: alignment), \
                    (default \'1\' for blast search)",
                    metavar="1 or 2",
                    choices=[1, 2, 3],
                    action='store', default=1, type=int)
# optional parameters
parser.add_argument("--fa",
                    help="input format of metagenomes sequence", type=str,
                    default='.fasta',metavar='.fasta or .fastq')
#parser.add_argument("--orf",
#                    help="input format of genome orfs", type=str, default='.genes.faa',metavar='.faa')
# optional output setup
parser.add_argument("--r",
                    help="output directory or folder of your results",
                    type=str, default='Result_traits',metavar='Result_traits')
parser.add_argument("--r16",
                    help="output directory or folder of your 16S sequences",
                    type=str, default='Result_16S',metavar='Result_16S')
# optional search parameters
parser.add_argument('--t',
                    help="Optional: set the thread number assigned for running XXX (default 1)",
                    metavar="1 or more", action='store', default=1, type=int)
parser.add_argument('--id',
                    default=75.0, action='store', type=float, metavar='60.0',
                    help='Optional: set the amno acid based identity cutoff for blast (default is 80.0)\n'
                         'Leave it alone if hmm is used')
parser.add_argument('--ht',
                    default=75.0, action='store', type=float, metavar='60.0',
                    help='Optional: set the amno acid based hit-length cutoff for blast (default is 80.0)\n'
                         'Leave it alone if hmm is used')
parser.add_argument('--e',
                    default=1e-5, action='store', type=float, metavar='1e-5',
                    help='Optional: set the evalue cutoff for blast or hmm (default is 1e-5)')
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
parser.add_argument('--hmm',
                    help="Optional: complete path to hmmscan if not in PATH,",
                    metavar="/usr/local/bin/hmmscan",
                    action='store', default='hmmscan', type=str)
parser.add_argument('--bp',
                    help="Optional: complete path to blast if not in PATH,",
                    metavar="/usr/local/bin/blast",
                    action='store', default='blast', type=str)
parser.add_argument('--bwa',
                    help="Optional: complete path to bwa if not in PATH,",
                    metavar="/usr/local/bin/bwa",
                    action='store', default='None', type=str)
parser.add_argument('--mini',
                          help="Optional: complete path to minimap2 if not in PATH,",
                          metavar="/usr/local/bin/minimap2",
                          action='store', default='None', type=str)
parser.add_argument('--bcf',
                    help="Optional: complete path to bcftools if not in PATH,",
                    metavar="/usr/local/bin/bcftools",
                    action='store', default='bcftools', type=str)
parser.add_argument('--sam',
                    help="Optional: complete path to bwa if not in PATH,",
                    metavar="/usr/local/bin/samtools",
                    action='store', default='samtools', type=str)
parser.add_argument('--vcf',
                    help="Optional: complete path to bwa if not in PATH,",
                    metavar="/usr/local/bin/vcftools",
                    action='store', default='vcftools', type=str)
parser.add_argument('--vcfstats',
                    help="Optional: complete path to vcfstats if not in PATH,",
                    metavar="/usr/local/bin/vcfstats",
                    action='store', default='vcfstats', type=str)
parser.add_argument('--strainfinder',
                    help="Optional: complete path to strainfinder",
                    metavar="/scratch/users/anniz44/bin/miniconda3/bin/strainfinder",
                    action='store',
                    default='/scratch/users/anniz44/bin/miniconda3/bin/strainfinder',
                    type=str)


################################################## Definition ########################################################
args = parser.parse_args()
fasta_format = args.fa
#orfs_format = args.orf
try:
    os.mkdir(args.r)
except OSError:
    pass
workingdir=os.path.abspath(os.path.dirname(__file__))


################################################### Function ########################################################
def split_string_last(input_string,substring):
    last_loci = input_string.rfind(substring)
    if last_loci > -1:
        return input_string[0 : last_loci]
    else:
        return input_string


def bowtie(database, metagenomes, tempbamoutput):
    # bowtie alignment
    # run alignment
    cmds = ''
    try:
        ftest = open('%s.sorted.bam' %(tempbamoutput),'r')
    except IOError:
        if args.bwa!='None':
            cmds += args.bwa + ' mem -a -p -t %s %s %s |%s view -@ %s -S -b -F 4 >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
                min(args.t,40), database, metagenomes, args.sam, min(args.t,40),
                tempbamoutput, args.sam, min(args.t,40), tempbamoutput, tempbamoutput, args.sam, min(args.t,40), tempbamoutput)
        if args.mini !='None':
            cmds += args.mini + ' -ax sr -N 1000 -p 0.9 -t %s %s.mmi %s |%s view -@ %s -S -b -F 4 >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
                min(args.t, 40), database, metagenomes, args.sam, min(args.t, 40),
                tempbamoutput, args.sam, min(args.t, 40), tempbamoutput, tempbamoutput, args.sam, min(args.t, 40),
                tempbamoutput)
    try:
        ftest = open('%s.flt.vcf' %(tempbamoutput),'r')
    except IOError:
        cmds += '%s mpileup --threads %s -a FMT/AD -q30 -B -Ou -d3000 -f %s %s.sorted.bam  | %s call --ploidy 1 -A --threads %s -m > %s.raw.vcf' % (
            args.bcf, min(args.t, 40), database, tempbamoutput, args.bcf, min(args.t, 40), tempbamoutput)
        cmds += '\n%s filter --threads %s -s LowQual -e \'DP>100\' %s.raw.vcf > %s.flt.vcf \n' % (
            args.bcf,min(args.t,40), tempbamoutput, tempbamoutput)
        cmds += '\n%s view -H -v snps %s.flt.vcf > %s.flt.snp.vcf \n' % (
            args.bcf, tempbamoutput, tempbamoutput)
    # calculate coverage
    try:
        ftest = open('%s.sorted.bam.cov' % (tempbamoutput), 'r')
    except IOError:
        cmds += '%s depth -a -Q 10 %s.sorted.bam > %s.sorted.bam.cov\n' % (
            args.sam, tempbamoutput, tempbamoutput)
        cmds += 'echo -e "Ref_ID\\tCov_length\\tAverage\\tStdev" > %s.sorted.bam.avgcov\n' % (tempbamoutput)
        cmds += '%s depth %s.sorted.bam |  awk \'{sum[$1]+=$3; sumsq[$1]+=$3*$3; count[$1]++} END { for (id in sum) { print id,"\t",count[id],"\t",sum[id]/count[id],"\t",sqrt(sumsq[id]/count[id] - (sum[id]/count[id])**2)}}\' >> %s.sorted.bam.avgcov\n' % (
            args.sam, tempbamoutput, tempbamoutput)
    # optional cleanup
    cmds += 'rm -rf %s.bam\n' % (tempbamoutput)
    return cmds


def search(roottemp,filename):
    # search the database for each file
    cmds = "#!/bin/bash \nmodule add c3ddb/blast+/2.7.1 \n"
    # using blastp
    if args.s == 1:
        if args.u != 'None' or args.hs != 'None' or args.dm != 'None':
            # two-step search
            Usearch = 0
            for root in glob.glob(os.path.join(args.r + '/usearch', '*')):
                try:
                    ftry = open(os.path.join(root, filename + '.usearch.txt.aa'), 'r')
                    searchfile = os.path.join(root, filename + '.usearch.txt.aa')
                    Usearch = 1
                    break
                except (IOError,FileNotFoundError):
                    pass
            if Usearch == 0:
                if args.dm != 'None' and args.dbf == 2:
                    # Start search target genes by diamond
                    cmds += split_string_last(args.dm, 'diamond') + "diamond blastx --query " + os.path.join(roottemp, filename) + \
                            " --db " + split_string_last(args.db, '.dmnd') + ".dmnd --out " + os.path.join(args.r + '/usearch/' + str(int(i/10000)),
                                                                          filename + '.usearch.txt') + \
                            " --outfmt 6 --max-target-seqs 1 --evalue " + str(args.e) + " --threads " + str(
                        min(int(i_max),40)) + " \n"
                    cmds += 'python ' + workingdir + '/Extract.MG.py -p 1 -i ' + roottemp + ' -f ' + filename + ' -n .usearch.txt -r ' + args.r + '/usearch/' + str(
                        min(int(i_max), 40)) + ' \n'
                    searchfile = os.path.join(args.r + '/usearch/' + str(int(i / 10000)), filename + '.usearch.txt.aa')
                elif args.hs != 'None' and args.dbf == 1:
                    # Start search target genes by hs-blastn
                    cmds += "%s align -db %s -window_masker_db %s.counts.obinary -query %s -out %s -outfmt 6 -evalue %s -num_threads %s\n" \
                            % (args.hs, args.db, args.db, os.path.join(roottemp, filename),
                               os.path.join(args.r + '/usearch/' + str(int(i / 10000)), filename + '.usearch.txt'),
                               str(args.e), str(min(int(i_max), 40)))
                    cmds += 'python ' + workingdir + '/Extract.MG.py -p 1 -i ' + roottemp + ' -f ' + filename + ' -n .usearch.txt -r ' + args.r + '/usearch/' + str(
                        int(i / 10000)) + ' \n'
                    searchfile = os.path.join(args.r + '/usearch/' + str(int(i / 10000)), filename + '.usearch.txt.aa')
                elif args.u != 'None':
                    # Start search target genes by usearch
                    if args.dbf == 1:
                        cmds += args.u + " -ublast " + os.path.join(roottemp, filename) + \
                                " -db " + split_string_last(args.db, '.udb') + ".udb -strand both -evalue 1e-2 -accel 0.5 -blast6out " \
                                + os.path.join(args.r + '/usearch/' + str(int(i/10000)), filename + '.usearch.txt') + \
                                " -threads " + str(int(i_max)) + " \n"
                    else:
                        cmds += args.u + " -ublast " + os.path.join(roottemp, filename) + \
                                " -db " + split_string_last(args.db, '.udb') + ".udb -evalue 1e-2 -accel 0.5 -blast6out " \
                                + os.path.join(args.r + '/usearch/' + str(int(i/10000)), filename + '.usearch.txt') + \
                                " -threads " + str(int(i_max)) + " \n"
                    cmds += 'python ' + workingdir + '/Extract.MG.py -p 1 -i ' + roottemp + ' -f ' + filename + ' -n .usearch.txt -r ' + args.r + '/usearch/' + str(
                        int(i / 10000)) + ' \n'
                    searchfile = os.path.join(args.r + '/usearch/' + str(int(i / 10000)), filename + '.usearch.txt.aa')
                else:
                    print('wrong search method for this database\nrun blat directly')
                    searchfile = os.path.join(roottemp, filename)
        else:
            # one-step search
            searchfile = os.path.join(roottemp, filename)
        if args.bp != 'None':
            # blast search
            Blastsearch = 0
            for root in glob.glob(os.path.join(args.r + '/search_output', '*')):
                try:
                    ftry_blast = open(os.path.join(root, filename + '.blast.txt'), 'r')
                    ftry_blast_file = os.path.join(root, filename + '.blast.txt')
                    Blastsearch = 1
                    break
                except (IOError,FileNotFoundError):
                    pass
            if Blastsearch == 0:
                # for short metagenomic reads
                if args.dbf == 1:
                    cmds += split_string_last(args.bp, 'blast') + "blastn -query " + str(searchfile) + " -db " + args.db + " -out " + args.r + '/search_output/'+str(int(i/10000))+ \
                             "/"+filename+".blast.txt  -outfmt 6 -evalue "+str(args.e)+" -num_threads " + \
                            str(min(int(i_max), 40)) + " \n"
                else:
                    cmds += split_string_last(args.bp, 'blast') + "blastx -query " + str(
                        searchfile) + " -db " + args.db + " -out " + args.r + '/search_output/' + str(int(i / 10000)) + \
                            "/" + filename + ".blast.txt  -outfmt 6 -evalue " + str(args.e) + " -num_threads " + \
                            str(min(int(i_max), 40)) + " \n"
                # fiter blast result
            Blastsearchfilter = 0
            for root in glob.glob(os.path.join(args.r + '/search_output', '*')):
                try:
                    ftry = open(os.path.join(root, filename + '.blast.txt.filter'), 'r')
                    Blastsearchfilter = 1
                    break
                except (IOError,FileNotFoundError):
                    pass
            if Blastsearchfilter == 0:
                if Blastsearch == 0:
                    # no blast result
                    tempbamoutput_filter = args.r + '/search_output/' + str(
                        int(i / 10000))
                else:
                    # blast completed
                    tempbamoutput_filter = os.path.split(ftry_blast_file)[0]
            if Blastsearchfilter == 0:
                cmds += 'python '+ workingdir +'/Filter.MG.py --g T -i ' + tempbamoutput_filter +' -f ' + filename + '.blast.txt ' +\
                        '-db ' + args.db + ' -dbf ' + str(args.dbf) + ' -s ' + str(args.s) + ' --ht ' + str(args.ht) + ' --id ' + str(args.id) + \
                        ' --e ' + str(args.e) + ' \n'
                cmds += 'python '+ workingdir +'/Extract.MG.py -p 1  -ni .usearch.txt.aa -i ' +\
                        os.path.split(searchfile)[0] + ' -f ' + os.path.split(searchfile)[1] + ' -n .blast.txt.filter -r ' + tempbamoutput_filter +' \n'
        # bowtie alignment
        if args.bwa != 'None' or args.mini != 'None':
            tempinput = os.path.join(args.r + '/search_output/' + str(int(i / 10000)),
                                     filename + '.blast.txt.filter.aa')
            tempbamoutput = os.path.join(args.r + '/bwa/' + str(int(i / 10000)), str(
                filename) + '.blast.txt.filter.aa')
            try:
                f1 = open('%s.sorted.bam' % (tempbamoutput))
            except (IOError,FileNotFoundError):
                cmds += bowtie(args.db, tempinput, tempbamoutput)
    elif args.s == 2:
        # hmmsearch
        Blastsearch = 0
        for root in glob.glob(os.path.join(args.r + '/search_output', '*')):
            try:
                ftry = open(os.path.join(root, filename + '.hmm'), 'r')
                Blastsearch = 1
                break
            except (IOError,FileNotFoundError):
                pass
        if Blastsearch == 0:
            cmds = args.hmm + ' --tblout ' + os.path.join(args.r + '/search_output/'+str(int(i/10000)), str(
                filename) + '.hmm') +  ' --cpu ' + str(int(i_max)) + ' -E ' \
                  + str(args.e) + ' ' +split_string_last(args.db, '.hmm') + '.hmm '+ os.path.join(roottemp, filename) + ' \n'
            cmds += 'python '+ workingdir +'/Format.MG.py -i ' + args.r + '/search_output/'+str(int(i/10000)) + ' -f ' + str(
                    filename) + '.hmm \n'
    elif args.s == 3 and (args.bwa != 'None' or args.mini != 'None'):
        tempinput = os.path.join(roottemp, filename)
        if '_2' + fasta_format not in tempinput:
            tempbamoutput = os.path.join(args.r + '/bwa/' + str(int(i / 10000)), str(filename))
            Bamfile = 0
            Covfile = 0
            Avgcovfile = 0
            tempinput2 = tempinput.replace('_1' + fasta_format, '_2' + fasta_format)
            tempbamoutput2 = os.path.join(args.r + '/bwa/' + str(int(i / 10000)), str(
                filename.replace('_1' + fasta_format, '_2' + fasta_format)))
            if tempinput2 in Targetroot and tempinput2!=tempinput:
                #tempinput = '%s %s' % (tempinput, tempinput2)
                cmds += bowtie(args.db, tempinput2, tempbamoutput2)
            cmds += bowtie(args.db, tempinput, tempbamoutput)
    else:
        print('please provide --bwa for alignment (--s 3)')
    # 16S extraction
    Search16s = 0
    for root in glob.glob(os.path.join(args.r16, '*')):
        try:
            ftry = open(os.path.join(root, filename + '.16S.txt'), 'r')
            Search16s = 1
            break
        except (IOError,FileNotFoundError):
            pass
    if Search16s == 0:
        if args.hs != 'None':
            # with hs-blastn
            # Start search 16S
            cmds += 'python '+ workingdir +'/undone.MG.py -i '+ os.path.join(roottemp, filename) +' \n'
            cmds += "%s align -db %s -window_masker_db %s.counts.obinary -query %s -out %s -outfmt 6 -evalue %s -num_threads %s\n" \
                    % (args.hs, workingdir + "/../database/85_otus.fasta.all.V4_V5.fasta",
                    workingdir + "/../database/85_otus.fasta.all.V4_V5.fasta",
                       os.path.join(args.r16+'/' + str(int(i/10000)), filename+ '.16S.txt'),
                       os.path.join(args.r16+'/' + str(int(i/10000)), filename+ '.16S.txt'),
                       str(args.e), str(min(int(i_max),40)))
            cmds += 'python '+ workingdir +'/Extract.16S.MG.py -i ' + roottemp + ' -f ' + \
                    filename + ' -n .16S.txt -r ' + args.r16 + '/' + str(
                int(i / 10000)) + ' \n'
            # with usearch
        elif args.u != 'None':
                # Start search 16S
                cmds += 'python ' + workingdir + '/undone.MG.py -i ' + os.path.join(roottemp, filename) + ' \n'
                cmds += args.u + " -usearch_global " + os.path.join(roottemp, filename) + \
                        " -db " + workingdir + "/../database/85_otus.fasta.all.V4_V5.fasta.udb -strand plus -id 0.7 -evalue 1e-1 -blast6out " \
                        + os.path.join(args.r16 + '/' + str(int(i / 10000)), filename + '.16S.txt') + \
                        " -threads " + str(int(i_max)) + " \n"
                cmds += 'python ' + workingdir + '/Extract.16S.MG.py -i ' + roottemp + ' -f ' + \
                        filename + ' -n .16S.txt -r ' + args.r16 + '/' + str(
                    int(i / 10000)) + ' \n'
    # cell number calculation
    Search16s = 0
    for root in glob.glob(os.path.join(args.r16, '*')):
        try:
            ftry = open(os.path.join(root, filename + '.uscmg.blastx.txt'), 'r')
            Search16s = 1
            break
        except (IOError,FileNotFoundError):
            pass
    if Search16s == 0:
        # with diamond
        if args.dm != 'None':
            # Start search essential single copy genes
            cmds += 'python ' + workingdir + '/undone.MG.py -i ' + os.path.join(roottemp,
                                                                                filename ) + ' \n'
            cmds += split_string_last(args.dm, 'diamond') + "diamond blastx -q " + os.path.join(roottemp, filename) +\
                     " -d " + workingdir + "/../database/all_KO30.pro.fa.dmnd -o " + \
                    os.path.join(args.r16 + '/' + str(int(i / 10000)), filename + '.uscmg.blastx.txt') + \
                     " -f tab  -p 1 -e 3  --id 0.45 --max-target-seqs 1 --threads " + str(int(i_max)) + " \n"
    return cmds


################################################### Programme #######################################################
# load all WGD
flist=open('Filelist.txt','w')
in_dir=args.i
Targetroot=dict()

Targetlist=[]
if args.l != 'None':
    for lines in open(args.l,'r'):
        Targetlist.append(split_string_last(split_string_last(lines, '\r'),'\n'))

for root in glob.glob(os.path.join(in_dir, '*')) + glob.glob(in_dir):
    list_fasta1 = glob.glob(os.path.join(root, '*'+fasta_format))
    if list_fasta1!=[]:
        for files in list_fasta1:
            if args.l == 'None' or any(targets in files for targets in Targetlist):
                Targetroot.setdefault(files, fasta_format)


# set up min coverage depth
MIN_COVERAGE_DEPTH=2


# search the database in all metagenomes
i=0
try:
    os.mkdir('subscripts')
except OSError:
    pass
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
if args.bwa != 'None' or args.mini !='None':
    try:
        os.mkdir(args.r + '/bwa/')
    except OSError:
        pass
i_max=int(args.t)
for files in Targetroot:
    if Targetroot[files]!='None':
        f1 = open(os.path.join('subscripts',str(i%int(args.t)) + '.sh'), 'a')
        i +=1
        roottemp, filename = os.path.split(files)
        flist.write(str(files) + '\n')
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
        if args.bwa != 'None' or args.mini !='None':
            try:
                os.mkdir(args.r + '/bwa/' + str(int(i / 10000)))
            except OSError:
                pass
        # search the database in WGD
        cmds = search(roottemp, filename)
        f1.write(cmds)
        f1.close()
flist.close()
