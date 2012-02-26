from optparse import OptionParser, OptionGroup
from bs_single_end import bs_single_end
from utils import *

if __name__ == '__main__':

    parser = OptionParser()

    opt_group = OptionGroup(parser, "For single end reads")
    opt_group.add_option("-i", "--input", type="string", dest="infilename",help="Input your read file name (FORMAT: sequences, illumina fastq, qseq,fasta)", metavar="INFILE")
    parser.add_option_group(opt_group)


    opt_group = OptionGroup(parser, "For pair end reads")
    opt_group.add_option("-1", "--input_1", type="string", dest="infilename_1",help="Input your read file end 1 (FORMAT: sequences, illumina fastq, qseq)", metavar="FILE")
    opt_group.add_option("-2", "--input_2", type="string", dest="infilename_2",help="Input your read file end 2 (FORMAT: sequences, illumina fastq, qseq)", metavar="FILE")

    opt_group.add_option("--minins",type = "int",dest = "min_insert_size", help="The minimum insert size for valid paired-end alignments [-1]", default = -1)
    opt_group.add_option("--maxins",type = "int",dest = "max_insert_size", help="The maximum insert size for valid paired-end alignments [400]", default = 400)

    parser.add_option_group(opt_group)

    opt_group = OptionGroup(parser, "General options")

    opt_group.add_option("-t", "--tag", type="string", dest="taginfo",help="Yes for undirectional lib, no for directional [Y]", metavar="TAG", default = 'Y')

    opt_group.add_option("-s","--start_base",type = "int",dest = "cutnumber1", help="The first base of your read to be mapped [1]", default = 1)

    opt_group.add_option("-e","--end_base",type = "int",dest = "cutnumber2", help="The last cycle number of your read to be mapped [36]", default = 36)

    opt_group.add_option("-a", "--adapter", type="string", dest="adapterfilename",help="Input text file of your adaptor sequences (to be trimed from the 3'end of the reads). Input 1 seq for dir. lib., 2 seqs for undir. lib. One line per sequence", metavar="FILE", default = '')

    opt_group.add_option("-g", "--genome", type="string", dest="genome",help="Name of the reference genome (the same as the reference genome file in the preprocessing step) [ex. chr21_hg18.fa]")

    opt_group.add_option("-m", "--mis",type = "int", dest="indexname",help="Number of mismatches (0,1,...,read length) [3]", default = 3)

    opt_group.add_option("--path", type="string", dest="bowtiepath",help="Path to Bowtie [%s]" % default_bowtie_path, metavar="PATH", default = default_bowtie_path)

    opt_group.add_option("--db", type="string", dest="dbpath",help="Path to the reference genome library (generated in preprocessing genome) [%s]" % reference_genome_path, metavar="DBPATH", default = reference_genome_path)

    opt_group.add_option("-l", "--split_line",type = "int", dest="no_split",help="Number of lines per split (the read file will be split into small files for mapping. The result will be merged. [4000000]", default = 4000000)

    opt_group.add_option("-o", "--output", type="string", dest="outfilename",help="The name of output file [INFILE.bs(se|pe|rrbs)]", metavar="OUTFILE")

    parser.add_option_group(opt_group)


    # Pair-end options




    #----------------------------------------------------------------
    (options, args) = parser.parse_args()


    # if no options were given by the user, print help and exit
    import sys
    if len(sys.argv) == 1:
        print parser.print_help()
        exit(0)



    if options.infilename and (options.infilename_1 or options.infilename_2):
        error('-i and [-1|-2] options are exclusive. You should use only one of them.')


    if not (options.infilename or (options.infilename_1 and options.infilename_2)):
        error('You should set either -i or -1 and -2 options.')




    asktag=str(options.taginfo).upper()
    if asktag not in 'YN':
        error('-t option should be either Y or N, not %s' % asktag)


    adapter_file=options.adapterfilename

    cut1=options.cutnumber1
    cut2=options.cutnumber2

    no_small_lines=options.no_split

    int_no_mismatches=min(int(options.indexname),cut2)
    indexname=str(int_no_mismatches)

    bowtie_path=os.path.join(options.bowtiepath,'bowtie')

    genome = options.genome
    if genome is None:
        error('-g is a required option')

    db_path = os.path.join(options.dbpath, genome + '_' + asktag)

    if not os.path.isfile(os.path.join(db_path,'ref.json')):
        error(genome + ' cannot be found in ' + options.dbpath +'. Please, run the Preprocessing_genome to create it.')


    if options.infilename is not None:
        # single end reads
        bs_single_end(  options.infilename,
                        asktag,
                        adapter_file,
                        cut1,
                        cut2,
                        no_small_lines,
                        int_no_mismatches,
                        indexname,
                        bowtie_path,
                        db_path,
                        options.outfilename or options.infilename+'.bsse' # this is the output file name
                        )
    else:
        # pair end reads

        main_read_file_1=options.infilename_1
        main_read_file_2=options.infilename_2



