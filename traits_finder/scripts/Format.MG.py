import argparse
import os

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                    help="input dir", type=str, default='.',metavar='current dir (.)')
parser.add_argument("-f",
                    help="input filename", type=str, default='input.faa',metavar='input.faa')


################################################## Definition ########################################################
args = parser.parse_args()


################################################### Function #######################################################
def hmmformat(hmmresult):
    # format hmm results
    f1 = open(hmmresult + '2.txt', 'a')
    for line in open(hmmresult, 'r'):
        if str(line)[0] == '#':
            pass
        else:
            line = str(line).replace(' # ', '#')
            while line != str(line).replace('  ', ' '):
                line = str(line).replace('  ', ' ')
            line = str(line).replace(' ', '\t')
            line = str(line).replace('#', ' # ')
            filedir, filename = os.path.split(hmmresult)
            filename = filename.split('.hmm')[0]
            f1.write(filename + '_' + line)
    f1.close()


################################################### Programme #######################################################
hmmformat(os.path.join(args.i,args.f))
