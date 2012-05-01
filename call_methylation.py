#!/usr/bin/python

from optparse import OptionParser
from bs_utils.utils import *
import pysam

if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-i", "--input", type="string", dest="infilename",help="BAM output from bs_seeker2.py", metavar="INFILE")
    parser.add_option("-g", "--genome", type="string", dest="genome",help="Name of the reference genome (the same as the reference genome file in the preprocessing step) [ex. chr21_hg18.fa]")
    parser.add_option("--db", type="string", dest="dbpath",help="Path to the reference genome library (generated in preprocessing genome) [%default]" , metavar="DBPATH", default = reference_genome_path)
    parser.add_option("-o", "--output", type="string", dest="outfilename",help="The name of output file [INFILE.meth_calls]", metavar="OUTFILE")

    (options, args) = parser.parse_args()


    # if no options were given by the user, print help and exit
    if len(sys.argv) == 1:
        print parser.print_help()
        exit(0)

    if options.infilename is None:
        error('-i option is required')
    if not os.path.isfile(options.infilename):
        error('Cannot find input file: %s' % options.infilename)

    sorted_input_filename = options.infilename+'_sorted'
    pysam.sort(options.infilename, sorted_input_filename)
    sorted_input_filename += '.bam'
    pysam.index(sorted_input_filename)

    sorted_input = pysam.Samfile(sorted_input_filename)
    for column in sorted_input.pileup():
        print column