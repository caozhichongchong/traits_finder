import argparse
import vcfpy

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                    help="input vcf filename", type=str, default='input.vcf',metavar='input.vcf')
parser.add_argument("-o",
                    help="output filename", type=str, default='input.vcf.out',metavar='input.vcf.out')

################################################## Definition ########################################################
args = parser.parse_args()


################################################### Function #######################################################
def process_record(record,output):
    line = [record.CHROM, record.POS, record.REF]
    line += [alt.value for alt in record.ALT]
    line += [call.data.get('GT') or './.' for call in record.calls]
    output.write('\t'.join(map(str, line)))

def readvcf(input,output):
    # Build and print header
    vcf_file=vcfpy.Reader.from_path(input)
    #vcf_file2=vcfpy.Call(input.replace('.vcf','.sorted.bam'),input)
    output.write('#'+str(vcf_file.header.samples.names[0]) + '\n')
    header = ['#CHROM', 'POS', 'REF', 'ALT']
    output.write('\t'.join(header)+'\n')
    #print(vcf_file2.gt_alleles)
    for record in vcf_file:
        if record.is_snv(): #a single nucleotide variant
            process_record(record,output)
            output.write('\n')
            print(record.FORMAT)
            print(record.affected_end, record.affected_start, record.begin,record.end)
            print(record.QUAL)
            print(record.calls)
            print(record.is_snv())
            for calls in record.calls:
                print(calls.gt_bases,calls.gt_phase_char,calls.gt_type)
                print(calls.is_het,calls.is_phased,calls.is_variant)
                print(calls.plodity,calls.sample)
                print(calls.site,calls.is_filtered)
            print('\n\n')
################################################### Programme #######################################################
readvcf(args.i, open(args.o,'w'))
