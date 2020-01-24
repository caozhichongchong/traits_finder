
import os
from Bio import SeqIO
import argparse
import glob


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
                    help="input dir of WGD", type=str, default='.',metavar='current dir (.)')
parser.add_argument('-s',
                    help="set the method to search the your database \
                    (1: blast; 2: hmm; 3: alignment), \
                    (default \'1\' for blast search)",
                    metavar="1 or 2",
                    choices=[1, 2, 3],
                    action='store', default=1, type=int)
# optional parameters
parser.add_argument("--fa",
                    help="input format of genome sequence", type=str, default='.fna.add',metavar='.fasta, .fna or .fa')
parser.add_argument("--orf",
                    help="input format of genome orfs", type=str, default='.faa',metavar='.faa')
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
parser.add_argument('--mf','--mafft',
                          help="Optional: complete path to mafft if not in PATH,",
                          metavar="/usr/local/bin/mafft",
                          action='store', default='None', type=str)
parser.add_argument('--pro','--prodigal',
                        help="Optional: complete path to prodigal if not in PATH,",
                        metavar="/usr/local/bin/prodigal",
                        action='store', default='None', type=str)


################################################## Definition ########################################################
args = parser.parse_args()
fasta_format = args.fa
orfs_format = args.orf
try:
    os.mkdir(args.r)
except OSError:
    pass
workingdir=os.path.abspath(os.path.dirname(__file__))


################################################### Function ########################################################
def addname(file_name):
    Fasta_name = open(file_name, 'r')
    f = open(file_name + '.add', 'w')
    in_dir, input_file = os.path.split(file_name)
    for record in SeqIO.parse(Fasta_name, 'fasta'):
        if len(str(record.seq).replace(' ',''))>0:
            # remove empty ORF sequences, otherwise there could be a problem for blast, usearch and diamond
            f.write('>'+str(input_file)+'_'+str(record.id) + '\t' + str(record.description).replace('\t', ' ') + '\n' + str(
                str(record.seq)) + '\n')
    f.close()


def split_string_last(input_string,substring):
    last_loci = input_string.rfind(substring)
    if last_loci > -1:
        return input_string[0 : last_loci]
    else:
        return input_string


def search(roottemp,genome_output,orf_output):
    # search the database for each file
    cmds = ''
    cmds += "#!/bin/bash \nmodule add c3ddb/blast+/2.7.1 \n"
    genome_file =  os.path.join(roottemp, genome_output)
    orf_file = os.path.join(roottemp, orf_output)
    # using blastp
    if args.s == 1:
        if args.u != 'None' or args.hs != 'None' or args.dm != 'None':
            # two-step search
            Usearch=0
            for root, dirs, files in os.walk(args.r + '/usearch'):
                try:
                    ftry = open(os.path.join(root, genome_output + '.usearch.txt.aa'), 'r')
                    searchfile_genome = os.path.join(root, genome_output + '.usearch.txt.aa')
                    searchfile_orf = os.path.join(args.r + '/usearch/' + str(folder_id),
                                                  orf_output + '.usearch.txt.aa')
                    Usearch=1
                    break
                except IOError:
                    pass
            if Usearch == 0:
                searchfile_orf = os.path.join(args.r + '/usearch/' + str(folder_id),
                                              orf_output + '.usearch.txt.aa')
                searchfile_genome = os.path.join(args.r + '/usearch/' + str(folder_id),
                                                 genome_output + '.usearch.txt.aa')
                if args.dm != 'None' and args.dbf == 2:
                    # Start search target genes by diamond
                    # AA file
                    cmds += split_string_last(args.dm, 'diamond') + "diamond blastp --query " + orf_file + \
                            " --db " + split_string_last(args.db, '.dmnd') + ".dmnd --out " + os.path.join(
                        args.r + '/usearch/' + str(folder_id),
                        orf_output + '.usearch.txt') + \
                            " --outfmt 6 --max-target-seqs 1 --evalue " + str(args.e) + " --threads " + str(
                        int(i_max)) + " \n"
                    # genome file
                    cmds += split_string_last(args.dm, 'diamond') + "diamond blastx --query " + genome_file + \
                            " --db " + split_string_last(args.db, '.dmnd') + ".dmnd --out " + os.path.join(
                        args.r + '/usearch/' + str(folder_id),
                        genome_output + '.usearch.txt') + \
                            " --outfmt 6 --max-target-seqs 1 --evalue " + str(args.e) + " --threads " + str(
                        int(i_max)) + " \n"
                    # AA file
                    cmds += 'python ' + workingdir + '/Extract.WG.py -i ' \
                            + roottemp + ' -f ' + orf_output + ' -n .usearch.txt -r ' + args.r + '/usearch/' + str(
                        folder_id) + ' \n'
                    # genome file
                    cmds += 'python ' + workingdir + '/Extract.MG.py  -p 1 -i ' + roottemp + ' -f ' + \
                            genome_output + ' -n .usearch.txt -r ' + args.r + '/usearch/' + str(
                        folder_id) + ' \n'
                elif args.hs != 'None' and args.dbf == 1:
                    # Start search target genes by hs-blastn
                    if args.dbf == 1:
                        # genome file
                        cmds += "%s align -db %s -window_masker_db %s.counts.obinary -query %s -out %s -outfmt 6 -evalue %s -num_threads %s\n"\
                        %(args.hs, args.db, args.db, genome_file, os.path.join(
                            args.r + '/usearch/' + str(folder_id),
                            genome_output + '.usearch.txt'),
                          str(args.e),str(min(int(i_max),40)))
                        cmds += 'python ' + workingdir + '/Extract.MG.py -p 1 -i ' + roottemp + ' -f ' + genome_output + ' -n .usearch.txt -r ' + args.r + '/usearch/' + str(
                            folder_id) + ' \n'
                        searchfile_orf = orf_file
                elif args.u != 'None':
                    # Start search target genes by usearch
                    usearch_cmd = ".udb  -evalue 1e-2 -accel 0.5 -blast6out "
                    if args.dbf == 1:
                        usearch_cmd += " -strand both "
                    # genome file
                    cmds += args.u + " -ublast " + genome_file + \
                            " -db " + split_string_last(args.db, '.udb') + usearch_cmd \
                            + os.path.join(args.r + '/usearch/' + str(folder_id),
                                           genome_output + '.usearch.txt') + \
                            " -threads " + str(int(i_max)) + " \n"
                    # AA file
                    cmds += args.u + " -ublast " + orf_file + \
                            " -db " + split_string_last(args.db, '.udb') + usearch_cmd \
                            + os.path.join(args.r + '/usearch/' + str(folder_id), orf_output + '.usearch.txt') + \
                            " -threads " + str(int(i_max)) + " \n"
                    # AA file
                    cmds += 'python ' + workingdir + '/Extract.WG.py -i ' + roottemp + ' -f ' + orf_output + ' -n .usearch.txt -r ' + args.r + '/usearch/' + str(
                        folder_id) + ' \n'
                    # genome file
                    cmds += 'python ' + workingdir + '/Extract.MG.py  -p 1 -i ' + roottemp + ' -f ' + \
                            genome_output + ' -n .usearch.txt -r ' + args.r + '/usearch/' + str(
                        folder_id) + ' \n'
                else:
                    print('wrong search method for this database\nrun blast directly')
                    searchfile_genome = genome_file
                    searchfile_orf = orf_file
        else:
            # one-step search
            searchfile_genome = genome_file
            searchfile_orf = orf_file
        # blast search
        if args.bp != 'None':
            Blastsearch = 0
            for root, dirs, files in os.walk(args.r + '/search_output'):
                try:
                    ftry_blast_file = os.path.join(root, genome_output + '.blast.txt')
                    ftry_blast = open(ftry_blast_file, 'r')
                    Blastsearch = 1
                    break
                except IOError:
                    pass
            if Blastsearch == 0:
                if args.dbf == 2:
                    # protein database
                    cmds += split_string_last(args.bp, 'blast') + "blastp -query " + searchfile_orf + " -db " + args.db + " -out " + args.r + '/search_output/'+str(folder_id)+ \
                             "/"+orf_output+".blast.txt  -outfmt 6  -evalue "+str(args.e)+" -num_threads " + \
                            str(min(int(i_max),40)) + " \n"
                    cmds += split_string_last(args.bp, 'blast') + "blastx -query " + searchfile_genome  \
                            + " -db " + args.db + " -out " + args.r + '/search_output/' + str(folder_id) + \
                            "/" + genome_output \
                            + ".blast.txt  -outfmt 6 -evalue " + str(args.e) + " -num_threads " + \
                            str(min(int(i_max), 40)) + " \n"
                else:
                    # DNA database
                    cmds += split_string_last(args.bp, 'blast') + "tblastn -query " + searchfile_orf + " -db " + args.db + " -out " + args.r + '/search_output/' + str(folder_id) + \
                            "/" + orf_output + ".blast.txt  -outfmt 6  -evalue " + str(args.e) + " -num_threads " + \
                            str(min(int(i_max), 40)) + " \n"
                    cmds += split_string_last(args.bp, 'blast') + "blastn -query " + searchfile_genome  \
                            + " -db " + args.db + " -out " + args.r + '/search_output/' + str(folder_id) + \
                            "/" + genome_output \
                            + ".blast.txt  -outfmt 6 -evalue " + str(args.e) + " -num_threads " + \
                            str(min(int(i_max), 40)) + " \n"
            # fiter blast result
            Blastsearchfilter = 0
            for root, dirs, files in os.walk(args.r + '/search_output'):
                try:
                    ftry = open(os.path.join(root, genome_output + '.blast.txt.filter'), 'r')
                    Blastsearchfilter = 1
                    break
                except IOError:
                    pass
            if Blastsearchfilter == 0:
                if Blastsearch == 0:
                    # no blast result
                    tempbamoutput_filter = args.r + '/search_output/' + str(
                        folder_id)
                else:
                    # blast completed
                    tempbamoutput_filter = os.path.split(ftry_blast_file)[0]
                if args.dbf == 1:
                    # dna database
                    # genome file
                    cmds += 'python '+ workingdir +'/Filter.WG.py -i ' + tempbamoutput_filter + ' -f ' + genome_output + '.blast.txt ' + \
                            '-db ' + args.db + ' -s ' + str(args.s) + ' --ht ' + str(args.ht) + ' --id ' + str(args.id) + \
                            ' --e ' + str(args.e) + ' \n'
                    cmds += 'python '+ workingdir +'/Extract.MG.py -p 2 -d 500 -ni .usearch.txt.aa -i ' + os.path.split(searchfile_genome)[0] + ' -f ' + \
                            os.path.split(searchfile_genome)[1] + ' -n .blast.txt.filter -r ' + tempbamoutput_filter + ' \n'
                else:
                    # protein database
                    # AA file
                    cmds += 'python '+ workingdir +'/Filter.WG.py -i ' + tempbamoutput_filter + ' -f ' + orf_output + '.blast.txt ' + \
                            '-db ' + args.db + ' -s ' + str(args.s) + ' --ht ' + str(args.ht) + ' --id ' + str(args.id) + \
                            ' --e ' + str(args.e) + ' \n'
                    cmds += 'python '+ workingdir +'/Extract.WG.py -ni .usearch.txt.aa -i ' + os.path.split(searchfile_orf)[0]  \
                            + ' -f ' + os.path.split(searchfile_orf)[1] + ' -n .blast.txt.filter -r ' + tempbamoutput_filter + ' \n'
                    # genome file
                    cmds += 'python '+ workingdir +'/Filter.WG.py -i ' + tempbamoutput_filter + ' -f ' + genome_output + '.blast.txt ' + \
                            '-db ' + args.db + ' -s ' + str(args.s) + ' --ht ' + str(args.ht) + ' --id ' + str(args.id) + \
                            ' --e ' + str(args.e) + ' \n'
                    cmds += 'python '+ workingdir +'/Extract.MG.py  -p 2 -d 500 -ni .usearch.txt.aa  -i ' + os.path.split(searchfile_genome)[0] + ' -f ' + \
                            os.path.split(searchfile_genome)[1]\
                     + ' -n .blast.txt.filter -r ' + tempbamoutput_filter + ' \n'
                Blastsearchfilter = 1
            # bowtie alignment
            if args.bwa != 'None' or args.mf != 'None' and Blastsearchfilter == 1:
                tempinput = os.path.join(args.r + '/search_output/' + str(folder_id),
                                         genome_output + '.blast.txt.filter.aa')
                tempbamoutput = os.path.join(args.r + '/bwa/' + str(folder_id), str(
                    genome_output) + '.blast.txt.filter.aa')
                try:
                    f1 = open('%s.sorted.bam' % (tempbamoutput))
                except IOError:
                    if args.bwa != 'None':
                        cmds += args.bwa + ' mem %s %s |samtools view -S -b >%s.bam \nsamtools sort %s.bam -o %s.sorted.bam\n samtools index %s.sorted.bam\n' % (
                            args.db, tempinput,
                            tempbamoutput, tempbamoutput, tempbamoutput, tempbamoutput)
                        cmds += 'bcftools mpileup -Ou -f %s %s.sorted.bam  | bcftools call -mv > %s.vcf\n' % (
                            args.db, tempbamoutput, tempbamoutput)
                    elif args.mf != 'None':
                        # mafft for multiple alignment
                        cmds += args.bwa + ' --nuc --adjustdirection --quiet --retree 2 --maxiterate 100 --thread %s %s > %s.align \n' % (
                            str(int(i_max)), tempinput, tempbamoutput)
                        # transfer multiple alignment to vcf
                        cmds += 'snp-sites -v -o %s.align %s.vcf \n' % (
                            tempbamoutput, tempbamoutput)
                    # cmds += 'python '+ workingdir +'/VCF.reader.py -i %s.vcf -o %s.vcf.out\n' %(
                    #    tempbamoutput,tempbamoutput)
    elif args.s == 2:
        # hmmsearch
        Blastsearch = 0
        for root, dirs, files in os.walk(args.r + '/search_output'):
            try:
                ftry = open(os.path.join(root, orf_file + '.hmm'), 'r')
                Blastsearch = 1
                break
            except IOError:
                pass
        if Blastsearch == 0:
            cmds = args.hmm + ' --tblout ' + os.path.join(args.r + '/search_output/'+str(folder_id), str(
                orf_output) + '.hmm') +  ' --cpu ' + str(int(i_max)) + ' -E ' \
                  + str(args.e) + ' ' + split_string_last(args.db, '.hmm') + '.hmm '+ orf_file + ' \n'
            cmds += 'python '+ workingdir +'/Format.WG.py -i ' + args.r + '/search_output/'+str(folder_id) + ' -f ' + str(
                orf_output) + '.hmm \n'
    elif args.s == 3:
        # bowtie alignment
        tempinput = os.path.join(roottemp, str(
            genome_output))
        tempbamoutput = os.path.join(args.r + '/bwa/' + str(folder_id), str(
            genome_output))
        try:
            f1 = open('%s.sorted.bam' % (tempbamoutput))
        except IOError:
            if args.bwa != 'None':
                    cmds += args.bwa + ' mem %s %s |samtools view -S -b >%s.bam \nsamtools sort %s.bam -o %s.sorted.bam\n samtools index %s.sorted.bam\n' % (
                        args.db, tempinput,
                        tempbamoutput, tempbamoutput, tempbamoutput, tempbamoutput)
                    cmds += 'bcftools mpileup -Ou -f %s %s.sorted.bam  | bcftools call -mv > %s.vcf\n' % (
                        args.db, tempbamoutput, tempbamoutput)
            elif args.mf != 'None':
                    # mafft for multiple alignment
                    cmds += args.bwa + ' --nuc --adjustdirection --quiet --retree 2 --maxiterate 100 --thread %s %s > %s.align \n' % (
                        str(int(i_max)), tempinput, tempbamoutput)
                    # transfer multiple alignment to vcf
                    cmds += 'snp-sites -v -o %s.align %s.vcf \n' % (
                        tempbamoutput, tempbamoutput)
            else:
                print('please provide --bwa or --mafft for alignment (--s 3)')
    # 16S extraction
    Search16s = 0
    for root, dirs, files in os.walk(args.r16):
        try:
            ftry = open(os.path.join(root, genome_output + '.16S.txt'), 'r')
            Search16s = 1
            break
        except IOError:
            pass
    if Search16s == 0:
        if args.hs != 'None':
            # with hs-blastn
            # genome file
            # Start search 16S
            cmds += 'python '+ workingdir +'/undone.WG.py -i '+ os.path.join(roottemp, genome_output+' \n')
            cmds += "%s align -db %s -window_masker_db %s.counts.obinary -query %s -out %s -outfmt 6 -evalue %s -num_threads %s\n" \
                    % (args.hs, workingdir +"/../database/85_otus.fasta",
                    workingdir +"/../database/85_otus.fasta", genome_file, os.path.join(
                args.r16+'/' + str(folder_id),
                genome_output + '.16S.txt'),
                       str(args.e), str(min(int(i_max),40)))
            cmds += 'python '+ workingdir +'/Extract.16S.WG.py -i ' + roottemp + ' -f ' + \
                    genome_output + ' -n .16S.txt -r ' + args.r16 + '/' + str(
                folder_id) + ' \n'
        # with usearch
        elif args.u != 'None':
            # Start search 16S
            cmds += 'python '+ workingdir +'/undone.WG.py -i '+ os.path.join(roottemp, genome_output +' \n')
            cmds += args.u + " -usearch_global " + genome_file + \
                    " -db "+ workingdir +"/../database/85_otus.fasta.udb -strand plus -id 0.7 -evalue 1e-1 -blast6out " \
                    + os.path.join(args.r16+'/' + str(folder_id), genome_output + '.16S.txt')  + \
                    " -threads " + str(int(i_max)) + " \n"
            cmds += 'python '+ workingdir +'/Extract.16S.WG.py -i ' + roottemp + ' -f ' + \
                    genome_output + ' -n .16S.txt -r ' + args.r16 + '/' + str(
                folder_id) + ' \n'
    return cmds


################################################### Programme #######################################################
# search the database in all genomes
i=0
i_max = int(args.t)
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
if args.bwa != 'None':
    try:
        os.mkdir(args.r + '/bwa/')
    except OSError:
        pass
folder_id = 0
try:
    os.mkdir(args.r + '/search_output/' + str(folder_id))
except OSError:
    pass
try:
    os.mkdir(args.r + '/usearch/' + str(folder_id))
except OSError:
    pass
try:
    os.mkdir(args.r16 + '/' + str(folder_id))
except OSError:
    pass
if args.bwa != 'None':
    try:
        os.mkdir(args.r + '/bwa/' + str(folder_id))
    except OSError:
        pass
# load all WGD
flist=open('Filelist.txt','w')
flist_list = []
in_dir = args.i
for root,dirs,files in os.walk(in_dir):
    if fasta_format != 'None':
        list_fasta1 = glob.glob(os.path.join(root, '*' + fasta_format))
        if list_fasta1!=[]:
            for genomefile in list_fasta1:
                # set output scripts
                f1 = open(os.path.join('subscripts', str(i % int(args.t)) + '.sh'), 'a')
                i += 1
                roottemp, genomefilename = os.path.split(genomefile)
                flist_list.append(str(genomefile))
                # mkdir dir
                folder_id = int(i / 10000)
                if folder_id > int(i - 1 / 10000):
                    try:
                        os.mkdir(args.r + '/search_output/' + str(folder_id))
                    except OSError:
                        pass
                    try:
                        os.mkdir(args.r + '/usearch/' + str(folder_id))
                    except OSError:
                        pass
                    try:
                        os.mkdir(args.r16 + '/' + str(folder_id))
                    except OSError:
                        pass
                    if args.bwa != 'None':
                        try:
                            os.mkdir(args.r + '/bwa/' + str(folder_id))
                        except OSError:
                            pass
                # set ORF file
                orffile = split_string_last(genomefile, fasta_format) + orfs_format
                orffilename = split_string_last(genomefilename, fasta_format) + orfs_format
                try:
                    ftry_orf = open(orffile, 'r')
                except IOError:
                    if args.pro != 'None':
                        # predict AA
                        os.system('%s -q -a %s -i %s' % (args.pro, orffile, genomefile))
                # add genomefilename to orf
                if ".add" not in orffile:
                    try:
                        # already added filename to orf
                        ftry = open(orffile + ".add", 'r')
                    except IOError:
                        # add filename to orf
                        addname(orffile)
                    orffile = orffile + ".add"
                    orffilename = orffilename + ".add"
                # search the database in WGD
                cmds = search(roottemp, genomefilename, orffilename)
                f1.write(cmds)
                f1.close()
flist.write('\n'.join(flist_list))
flist.close()
