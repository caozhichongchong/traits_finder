import os
import argparse
import glob

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-t",
                    help="trait name", type=str, default='ARG',metavar='trait name')
parser.add_argument("--r",
                    help="input directory or folder of your previous results by Traits_WGD.py",
                    type=str, default='Result',metavar='Result', nargs='+')

################################################## Definition ########################################################
files_to_merge=dict()
args = parser.parse_args()


# set merge results dir
result_dir = args.r[0]
merge_dir = os.path.join(result_dir,'merge')
try:
    os.mkdir(merge_dir)
except OSError:
    pass
print('merge results can be found in %s' % (merge_dir))
merge_dir = os.path.join(merge_dir,'summary')
try:
    os.mkdir(merge_dir)
except OSError:
    pass

# set all files for merging
# sequences, blast summary, 16S sequences
merge_list = [args.t + '.all.traits.aa.fasta',args.t + '.all.traits.dna.fasta',\
args.t + '.all.traits.dna.txt',args.t + '.all.traits.aa.txt',\
args.t + '.all.16S.fasta']


# set all summary files
result_dir = os.path.join(result_dir,'summary')
summary_files = glob.glob(os.path.join(result_dir,args.t+'.all.traits.*.summarize.*.txt'))
for files in summary_files:
    merge_list.append(os.path.split(files)[-1])


################################################### Programme #######################################################
# collect merge files
for result_dir in args.r:
    result_dir = os.path.join(result_dir,'summary')
    # search output files
    # extra 500 sequences
    extra_file=glob.glob(os.path.join(result_dir, args.t + '.all.traits.dna.extra*.fasta'))
    if extra_file != []:
        if args.t + '.all.traits.dna.extra500.fasta' not in files_to_merge:
            files_to_merge.setdefault(args.t + '.all.traits.dna.extra500.fasta',
            [extra_file[0]])
        else:
            files_to_merge[args.t + '.all.traits.dna.extra500.fasta'].append(
                extra_file[0])
    for file_name in merge_list:
        try:
            ftest = open(os.path.join(result_dir, file_name),'r')
            if file_name not in files_to_merge:
                files_to_merge.setdefault(file_name,
                [os.path.join(result_dir, file_name)])
            else:
                files_to_merge[file_name].append(
                os.path.join(result_dir, file_name))
        except IOError:
            pass


for file_name in files_to_merge:
    os.system('cat %s > %s'%(' '.join(files_to_merge[file_name]),os.path.join(merge_dir,file_name)))
