import glob
import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                    help="input dir of traits", type=str,
                    default='/scratch/users/anniz44/Metagenomes/NR/search_output/0',
                    metavar='current dir (.)')
parser.add_argument("-i16",
                    help="input dir of traits", type=str,
                    default='/scratch/users/anniz44/Metagenomes/16S/0',
                    metavar='current dir (.)')
parser.add_argument("-o",
                    help="output directory",
                    type=str, default='/scratch/users/anniz44/Metagenomes/NR_summary', metavar='output')
parser.add_argument("-f",
                    help="input format of traits",
                    type=str, default='.fasta.hmm2.txt', metavar='.fasta.hmm2.txt')
parser.add_argument("-f16",
                    help="input format of 16S",
                    type=str, default='.fasta.16S.txt', metavar='.fasta.16S.txt')
parser.add_argument("-l",
                    help="traits length",
                    type=str, default='/scratch/users/anniz44/scripts/database/NR.dna.new.hmm.length',
                    metavar='NR.dna.new.hmm.length')
parser.add_argument("-tf",
                    help="traits function",
                    type=str, default='/scratch/users/anniz44/scripts/database/NR.dna.new.hmm.mapping',
                    metavar='NR.dna.new.hmm.mapping')
parser.add_argument("-m",
                    help="habitat metadata",
                    type=str, default='/scratch/users/anniz44/scripts/traits_search_MG/MGD_meta_new.txt',
                    metavar='MGD_meta_new.txt')
parser.add_argument("-b",
                    help="blastresults",
                    type=str, default='/scratch/users/anniz44/scripts/traits_search_MG/alltraits.to.alltraits.90.0.dedupli.txt',
                    metavar='alltraits.to.alltraits.90.0.dedupli.txt.usearch')
parser.add_argument("-c",
                    help="cutoff for traits (evalue)",
                    type=float, default=1e-2,
                    metavar='1e-2')


################################################## Definition ########################################################
args = parser.parse_args()
try:
    os.mkdir(args.o)
except OSError:
    pass
try:
    os.mkdir(args.o+'/similaritynew')
except OSError:
    pass

################################################## Function ########################################################
def copycal(traitsfile,file16S,i):
    copy1 = copytraits(traitsfile,i)
    copy2 = copy16S(file16S)
    temp = []
    temp.append(str(copy2))
    for traitscopy in copy1:
        temp.append(str(traitscopy))
        if copy2 > 0:
            temp.append(str(traitscopy / copy2))
        else:
            print(file16S, copy2,traitsfile)
            temp.append(str(0))
    return temp


def findfunction(gene):
    return Functionlist[Function[gene]]


def copytraits(traitsfile,i):
    #totaltraits=[0,0,0,0,0]
    totaltraits=[]
    while i > 0:
        totaltraits.append(0)
        i = i-1
    for lines in open(traitsfile,'r'):
        if float(lines.split('\t')[-4]) <= args.c:
            gene = lines.split('\t')[2]
            totaltraits[findfunction(gene)] += float(lines.split('\t')[-6])/Length[gene]
    return totaltraits


def copy16S(file16S):
    total16S=0
    try:
        for lines in open(file16S,'r'):
            total16S += float(lines.split('\t')[3])/1522.0
    except FileNotFoundError:
        pass
    return total16S


def traits_filter(blastfile, aafile):
    f1 = open(aafile+'.'+str(args.c)+'.filter','w')
    filtered = []
    for lines in open(blastfile,'r'):
        if float(lines.split('\t')[-4]) <= args.c:
            filtered.append(lines.split('\t')[0])
    for record in SeqIO.parse(aafile, "fasta"):
        if record.id in filtered:
            f1.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
    f1.close()


def traits_dedupli(aallseqlist,blastfile):
    tempseq=[]
    for lines in open(blastfile,'r'):
        aallseqlist[lines.split('\t')[0]] = Function[lines.split('\t')[1]]


def traits_dedupli2(aseqlist, aallseqlist,blastfile,aafile):
    tempseq=[]
    f1 = open(aafile + '.dedupli', 'w')
    for record in SeqIO.parse(aafile, "fasta"):
        #aallseqlist.setdefault(str(record.id),'None')
        if str(record.seq) not in tempseq:
            #aseqlist.setdefault(str(record.id),'None')
            tempseq.append(str(record.seq))
            f1.write('>' + str(record.id) + '\n' + str(record.seq) + '\n')
    f1.close()


def habitat_tag(habitattype):
    if habitattype in ['fecal_animal','fecal_human']:
        return 'Fecal'
    elif habitattype in ['wwtp']:
        return 'WWTP'
    else:
        return 'Nature'

def habitat_pair(habitattype1,habitattype2,tagpair):
    tag1=habitat_tag(habitattype1)
    tag2 = habitat_tag(habitattype2)
    if tag1 == tag2:
        return tag1
    else:
        if tag1 + '__' + tag2 in tagpair:
            return tag1 + '__' + tag2
        elif tag2 + '__' + tag1 in tagpair:
            return tag2 + '__' + tag1
        else:
            tagpair.append(tag1 + '__' + tag2)
            return tag1 + '__' + tag2

def swap(key):
    # for subtype
    if 'wwtp_' in key:
        return 'wwtp_'+key.split('wwtp_')[1].split('_')[0]+'_'+key.replace('wwtp_'+key.split('wwtp_')[1].split('_')[0],'')
    elif 'natural_' in key:
        try:
            return 'natural_'+key.split('natural_')[1].split('_')[0]+'_'+key.split('natural_')[1].split('_')[1]+'_'+key.replace('natural_'+key.split('natural_')[1].split('_')[0]+'_'+key.split('natural_')[1].split('_')[1],'')
        except IndexError:
            return 'natural_'+key.split('natural_')[1].split('_')[0]+'_'+key.replace('natural_'+key.split('natural_')[1].split('_')[0],'')
    elif 'fecal_animal_' in key:
        return 'fecal_animal_'+key.split('fecal_animal_')[1].split('_')[0]+'_'+key.replace('fecal_animal_'+key.split('fecal_animal_')[1].split('_')[0],'')
    else:
        return key


def swap2(key):
    # for type
    if 'wwtp' in key:
        return 'wwtp'+'_'+key.replace('wwtp','')
    elif 'natural_' in key:
        return 'natural_'+key.split('natural_')[1].split('_')[0]+'_'+key.replace('natural_'+key.split('natural_')[1].split('_')[0],'')
    elif 'fecal_animal' in key:
        return 'fecal_animal'+'_'+key.replace('fecal_animal','')
    else:
        return key


def swap3(key):
    return key.replace('__','_').replace('fecal_animal','Animal_feces').replace('fecal_human','Human_feces').replace('natural_','Nature_').replace('wwtp','WWTP')


def traits_habitat(blastresult,allseqlist):
    habitatset=[]
    tagpair = []
    for lines in open(blastresult):
            if allseqlist[lines.split('\t')[0]] == allseqlist[lines.split('\t')[1]]:
                MG1 = lines.split('\t')[0].split('.')[0]
                MG2 = lines.split('\t')[1].split('.')[0]
                if lines.split('\t')[0] != lines.split('\t')[1]:
                    if Meta[MG1][0] == Meta[MG2][0]:
                        #file1 = open(os.path.join(args.o+'/similaritynew',Meta[MG1][0]+'.intrasubtype.similaritynew'),'a')
                        #file1.write(str(lines.split('\t')[2])+'\n')
                        #file1.close()
                        file1 = open(os.path.join(args.o+'/similaritynew','Intrasubtype.similaritynew'),'a')
                        file1.write(Meta[MG1][1]+'\t'+Meta[MG1][0]+'\t'+str(lines.split('\t')[2]) + '\tIntrasubtype\n')
                        file1.close()
                    else:
                        #file1 = open(os.path.join(args.o+'/similaritynew', Meta[MG1][0] + '_' + Meta[MG2][0] +'.intersubtype.similaritynew'), 'a')
                        #file1.write(str(lines.split('\t')[2]) + '\n')
                        #file1.close()
                        if Meta[MG1][0] + '_' + Meta[MG2][0] in habitatset:
                            habitatpair = Meta[MG1][0] + '_' + Meta[MG2][0]
                        elif Meta[MG2][0] + '_' + Meta[MG1][0] in habitatset:
                            habitatpair = Meta[MG2][0] + '_' + Meta[MG1][0]
                        else:
                            habitatpair = Meta[MG1][0] + '_' + Meta[MG2][0]
                            habitatset.append(habitatpair)
                        file1 = open(os.path.join(args.o+'/similaritynew', 'Intersubtype.similaritynew'), 'a')
                        file1.write(habitat_pair(Meta[MG1][1],Meta[MG2][1],tagpair)+'\t'+swap3(swap(habitatpair)) + '\t' + str(lines.split('\t')[2]) + '\tIntersubtype\n')
                        file1.close()
                    if Meta[MG1][1] == Meta[MG2][1]:
                        #file1 = open(os.path.join(args.o+'/similaritynew',Meta[MG1][1]+'.intratype.similaritynew'),'a')
                        #file1.write(str(lines.split('\t')[2])+'\n')
                        #file1.close()
                        file1 = open(os.path.join(args.o+'/similaritynew','Intratype.similaritynew'),'a')
                        file1.write(Meta[MG1][1] + '\t' + Meta[MG1][1] + '\t' + str(lines.split('\t')[2])+'\tIntratype\n')
                        file1.close()
                    else:
                        #file1 = open(os.path.join(args.o+'/similaritynew', Meta[MG1][1] + '_' + Meta[MG2][1] +'.intertype.similaritynew'), 'a')
                        #file1.write(str(lines.split('\t')[2]) + '\n')
                        #file1.close()
                        if Meta[MG1][1] + '_' + Meta[MG2][1] in habitatset:
                            habitatpair = Meta[MG1][1] + '_' + Meta[MG2][1]
                        elif Meta[MG2][1] + '_' + Meta[MG1][1] in habitatset:
                            habitatpair =Meta[MG2][1] + '_' + Meta[MG1][1]
                        else:
                            habitatpair = Meta[MG1][1] + '_' + Meta[MG2][1]
                            habitatset.append(habitatpair)
                        file1 = open(os.path.join(args.o+'/similaritynew', 'Intertype.similaritynew'), 'a')
                        file1.write(habitat_pair(Meta[MG1][1],Meta[MG2][1],tagpair)+'\t'+swap3(swap2(habitatpair))+'\t'+str(lines.split('\t')[2]) + '\tIntertype\n')
                        file1.close()


################################################## Programme ########################################################
# load traits length
Length=dict()
for lines in open(args.l,'r'):
    Length.setdefault(lines.split('\t')[0],
                      float(lines.split('\t')[-1].split('\r')[0].split('\n')[0]))

# load traits function
Function=dict()
Functionlist=dict()
genenum=0
for lines in open(args.tf,'r'):
    Function.setdefault(lines.split('\t')[0],
                      (lines.split('\t')[-1].split('\r')[0].split('\n')[0]))
    if (lines.split('\t')[-1].split('\r')[0].split('\n')[0]) not in Functionlist:
        Functionlist.setdefault((lines.split('\t')[-1].split('\r')[0].split('\n')[0]),genenum)
        genenum+=1

print(Function)
print(Functionlist)

# load metadata
Meta=dict()
for lines in open(args.m,'r'):
    Meta.setdefault(lines.split('\t')[-1].split('.')[0].split('\r')[0].split('\n')[0],
                    [lines.split('\t')[2],lines.split('\t')[1]])

# load all MG results
alltraitsfile = glob.glob(os.path.join(args.i,'*'+args.f))
fout = open(os.path.join(args.o,'traits_copy_per_16S.'+str(args.c)+'.txt'),'w')
fout.write('SampleID\t16Scopy\t')
for functions in Functionlist:
    fout.write(str(functions)+'copy\t'+str(functions)+'copyper16S\t')
fout.write('Habitat_subtype\tHabitat_type\n')
for files in alltraitsfile:
    filedir, filename = os.path.split(files)
    if '.fasta.16S.txt' not in filename:
        files16S = os.path.join(args.i16,filename.replace(args.f,args.f16))
        MG = filename.split('.')[0]
        try:
            temp = MG+'\t'+'\t'.join(copycal(files, files16S,genenum))
            fout.write(temp+'\t'+Meta[MG][0]+'\t'+Meta[MG][1]+'\n')
        except KeyError:
            fout.write(temp+'\n')
fout.close()

#for files in alltraitsfile:
    #filedir, filename = os.path.split(files)
    #traits_filter(files, files+'.aa')

#seqlist=dict()
#allseqlist=dict()
#for files in alltraitsfile:
    #traits_dedupli2(seqlist, allseqlist, files, files + '.aa.90.0.filter')
    #traits_dedupli(allseqlist,files)

#traits_habitat(args.b,allseqlist)
