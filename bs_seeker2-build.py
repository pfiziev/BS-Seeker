import fileinput, string,os, operator, shelve, time, subprocess
import json
import re
from optparse import OptionParser, OptionGroup
from rrbs_build import rrbs_build
from utils import *
from wg_build import wg_build


if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-f", "--file", dest="filename",help="Input your reference genome file (fasta)", metavar="FILE")

    parser.set_defaults(taginfo="Y")
    parser.add_option("-t", "--tag", dest="taginfo",help="Yes for undirectional lib, no for directional [Y]", metavar="TAG")

    parser.set_defaults(bowtiepath = default_bowtie_path)
    parser.add_option("-p", "--path", dest="bowtiepath",help="Path to Bowtie [%s]" % default_bowtie_path, metavar="PATH")

    parser.set_defaults(dbpath = reference_genome_path)
    parser.add_option("-d", "--db", type="string", dest="dbpath",help="Path to the reference genome library (generated in preprocessing genome) [%s]" % reference_genome_path, metavar="DBPATH")


    # RRBS options
    rrbs_opts = OptionGroup(parser, "Reduced Representation Bisulfite Sequencing Options",
                                "Use this options with conjuction of -r [--rrbs]")

    rrbs_opts.add_option("-r", "--rrbs", action="store_true", dest="rrbs", default = False, help = 'Preprocess the genome for analysis of Reduced Representation Bisulfite Sequencing experiments')

    rrbs_opts.add_option("-l", "--low", dest="low_bound",help="lower bound", default = 75)
    rrbs_opts.add_option("-u", "--up", dest="up_bound",help="upper bound", default = 280)
    parser.add_option_group(rrbs_opts)


    (options, args) = parser.parse_args()

    # if no options were given by the user, print help and exit
    import sys
    if len(sys.argv) == 1:
        print parser.print_help()
        exit(0)



    rrbs = options.rrbs


    fasta_file=options.filename
    if fasta_file is None:
        error('Fasta file for the reference genome must be supported')

    if not os.path.isfile(fasta_file):
        error('%s cannot be found' % fasta_file)

    asktag=str(options.taginfo).upper()

    if asktag not in 'YN':
        error('-t option should be either Y or N, not %s' % asktag)


    bowtie_path=os.path.join(options.bowtiepath,'bowtie-build')


    print "Reference genome file: %s" % fasta_file
    print "Reduced Representation Bisulfite Sequencing: %s" % rrbs
    print "BS reads from undirectional/directional library: %s" % asktag
    print "Bowtie path: %s" % bowtie_path
    #---------------------------------------------------------------

    ref_path = options.dbpath

    if os.path.exists(ref_path):
        if not os.path.isdir(ref_path):
            error("%s must be a directory. Please, delete it or change the -d option." % ref_path)
    else:
        os.mkdir(ref_path)


    if rrbs: # RRBS preprocessing
        rrbs_build(fasta_file, asktag, bowtie_path, ref_path, options.low_bound, options.up_bound)
    else: # Whole genome preprocessing
        wg_build(fasta_file, asktag, bowtie_path, ref_path)