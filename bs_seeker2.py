from optparse import OptionParser, OptionGroup
import re
from bs_pair_end import bs_pair_end
from bs_rrbs import bs_rrbs
from bs_single_end import bs_single_end
from utils import *
import sys




if __name__ == '__main__':

    parser = OptionParser()

    opt_group = OptionGroup(parser, "For single end reads")
    opt_group.add_option("-i", "--input", type="string", dest="infilename",help="Input your read file name (FORMAT: sequences, illumina fastq, qseq,fasta)", metavar="INFILE")
    parser.add_option_group(opt_group)


    opt_group = OptionGroup(parser, "For pair end reads")
    opt_group.add_option("-1", "--input_1", type="string", dest="infilename_1",help="Input your read file end 1 (FORMAT: sequences, illumina fastq, qseq)", metavar="FILE")
    opt_group.add_option("-2", "--input_2", type="string", dest="infilename_2",help="Input your read file end 2 (FORMAT: sequences, illumina fastq, qseq)", metavar="FILE")
    opt_group.add_option("--minins",type = "int",dest = "min_insert_size", help="The minimum insert size for valid paired-end alignments [%default]", default = -1)
    opt_group.add_option("--maxins",type = "int",dest = "max_insert_size", help="The maximum insert size for valid paired-end alignments [%default]", default = 400)

    parser.add_option_group(opt_group)

    opt_group = OptionGroup(parser, "Reduced Representation Bisulfite Sequencing Options")
    opt_group.add_option("-r", "--rrbs", action="store_true", dest="rrbs", default = False, help = 'Preprocess the genome for analysis of Reduced Representation Bisulfite Sequencing experiments')
    opt_group.add_option("--rrbs-tag", type="string",dest="rrbs_taginfo",help="Msp-I tag: CGG TGG CGA or CGG/TGG (both)", metavar="TAG", default = "CGG/TGG")

    opt_group.add_option("--low", dest="rrbs_low_bound",help="lower bound")
    opt_group.add_option("--up",  dest="rrbs_up_bound",help="upper bound")


    parser.add_option_group(opt_group)


    opt_group = OptionGroup(parser, "General options")

    opt_group.add_option("-t", "--tag", type="string", dest="taginfo",help="Yes for undirectional lib, no for directional [%default]", metavar="TAG", default = 'Y')

    opt_group.add_option("-s","--start_base",type = "int",dest = "cutnumber1", help="The first base of your read to be mapped [%default]", default = 1)

    opt_group.add_option("-e","--end_base",type = "int",dest = "cutnumber2", help="The last cycle number of your read to be mapped [%default]", default = 36)

    opt_group.add_option("-a", "--adapter", type="string", dest="adapter_file",help="Input text file of your adaptor sequences (to be trimed from the 3'end of the reads). Input 1 seq for dir. lib., 2 seqs for undir. lib. One line per sequence", metavar="FILE", default = '')

    opt_group.add_option("-g", "--genome", type="string", dest="genome",help="Name of the reference genome (the same as the reference genome file in the preprocessing step) [ex. chr21_hg18.fa]")

    opt_group.add_option("-m", "--mis",type = "int", dest="int_no_mismatches",help="Number of mismatches (0,1,...,read length) [%default]", default = 3)

    opt_group.add_option("--aligner", dest="aligner",help="Aligner program to perform the analisys: " + ', '.join(supported_aligners) + " [%default]", metavar="ALIGNER", default = BOWTIE)

    opt_group.add_option("-p", "--path",   dest="aligner_path",help="Path to the aligner program. Defaults: " +' '*70+ '\t'.join(('%s: %s '+' '*70) % (al, aligner_path[al]) for al in sorted(supported_aligners)),
        metavar="PATH"
    )


    opt_group.add_option("--db", type="string", dest="dbpath",help="Path to the reference genome library (generated in preprocessing genome) [%default]" , metavar="DBPATH", default = reference_genome_path)

    opt_group.add_option("-l", "--split_line",type = "int", dest="no_split",help="Number of lines per split (the read file will be split into small files for mapping. The result will be merged. [%default]", default = 4000000)

    opt_group.add_option("-o", "--output", type="string", dest="outfilename",help="The name of output file [INFILE.bs(se|pe|rrbs)]", metavar="OUTFILE")

    parser.add_option_group(opt_group)

    opt_group = OptionGroup(parser, "Aligner Options",
        "You may specify any additional options for the aligner. You just have to prefix them with " +
        ', '.join('%s for %s' % (aligner_options_prefixes[aligner], aligner) for aligner in supported_aligners)+
        ', and BS Seeker will pass them on. For example: --bt-p 4 will increase the number of threads for bowtie to 4, '
        '--bt--tryhard will instruct bowtie to try as hard as possible to find valid alignments when they exist, and so on. '
        'Be sure that you know what you are doing when using these options! Also, we don\'t do any validation on the values.')

    parser.add_option_group(opt_group)

    # Pair-end options


    #----------------------------------------------------------------
    # separate aligner options from BS Seeker options
    aligner_options = {}
    bs_seeker_options = []
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        m = re.match(r'^%s' % '|'.join('(%s)'% aligner_options_prefixes[al] for al in supported_aligners), arg)
        if m:
            a_opt = arg.replace(m.group(0),'-',1)
            aligner_options[a_opt] = []
            while i + 1 < len(sys.argv) and sys.argv[i+1][0] != '-':
                aligner_options[a_opt].append(sys.argv[i+1])
                i += 1
            if len(aligner_options[a_opt]) == 0: # if it is a key-only option
                aligner_options[a_opt] = True
        else:
            bs_seeker_options.append(arg)
        i += 1


    (options, args) = parser.parse_args(args = bs_seeker_options)


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


    if options.aligner not in supported_aligners:
        error('-a option should be: %s' % ' ,'.join(supported_aligners)+'.')

    aligner_exec = os.path.join(options.aligner_path or aligner_path[options.aligner], options.aligner)

    int_no_mismatches=min(options.int_no_mismatches, options.cutnumber2)
    indexname=str(int_no_mismatches)

    genome = options.genome
    if genome is None:
        error('-g is a required option')

    genome_subdir = genome + '_' + asktag

    # try to guess the location of the reference genome for RRBS
    if options.rrbs:
        if options.rrbs_low_bound and options.rrbs_up_bound:
            genome_subdir += '_rrbs_%d_%d'  % (options.rrbs_low_bound, options.rrbs_up_bound)
        else:
            possible_refs = filter(lambda dir: dir.startswith(genome+'_'+asktag+'_rrbs_'), os.listdir(options.dbpath))
            if len(possible_refs) == 1:
                genome_subdir = possible_refs[0]
            else:
                error('Cannot localize unambiguosly the reference genome for RRBS. '
                      'Please, specify the --low and --up options that you used at the preprocessing step.\n'
                      'Possible choices are:\n' + '\n'.join([pr.split('_rrbs_')[-1].replace('_',', ') for pr in possible_refs]))

    db_path = os.path.join(options.dbpath, genome_subdir + '_' + options.aligner)

    if not os.path.isfile(os.path.join(db_path,'ref.json')):
        error(genome + ' cannot be found in ' + options.dbpath +'. Please, run the bs_seeker2-build.py to create it.')


    # handle aligner options
    #

    # default aligner options
    aligner_options_defaults = {
                                BOWTIE  : { '-e'              : 40*int_no_mismatches,
                                            '--nomaqround'    : True,
                                            '--norc'          : True,
                                            '-k'              : 2,
                                            '--quiet'         : True,
                                            '--best'          : True,
                                            '--suppress'      : '2,5,6',
                                            '-p'              : 2
                                },
                                BOWTIE2 : {
                                            '-k'              : 2,
                                            '--norc'          : True,
                                            '--quiet'         : True,
                                            '-p'              : 2,
                                            '--sam-nohead'    : True,
                                            # run bowtie2 in local mode by default
                                            '--local' : '--end-to-end' not in aligner_options

                                },
                                SOAP    : {}
                                }


    aligner_options = dict(aligner_options_defaults[options.aligner], **aligner_options)

    aligner_options_string = lambda : ' %s ' % (' '.join(opt_key +
                                                         (' ' + ' '.join(map(str,opt_val)) # join all values if the value is an array
                                                          if type(opt_val) is list else
                                                                ('' if type(opt_val) is bool and opt_val # output an empty string if it is a key-only option
                                                                 else ' ' +str(opt_val)) # output the value if it is a single value
                                                         )
                                                        for opt_key, opt_val in aligner_options.iteritems() if opt_val not in [None, False]))


    tmp_path = (options.infilename or options.infilename_1) +'-'+ options.aligner+ '-TMP'
    clear_dir(tmp_path)


    print 'Reduced Representation Bisulfite Sequencing:', options.rrbs
    if options.infilename is not None:
        print 'Single end'


#        aligner_options = dict({BOWTIE  : {},
#                                BOWTIE2 : {},
#                                SOAP    : {}}[options.aligner],
#            **aligner_options)


        aligner_command = 'nohup ' + aligner_exec  + aligner_options_string() + { BOWTIE   : ' %(reference_genome)s  -f %(input_file)s %(output_file)s',
                                                                                  BOWTIE2  : ' -x %(reference_genome)s -f -U %(input_file)s -S %(output_file)s',
                                                                                  SOAP     : ' -D %(reference_genome)s -o %(output_file)s -a %(input_file)s'}[options.aligner]
        print 'Aligner command:', aligner_command
        # single end reads
        if options.rrbs: # RRBS scan
            bs_rrbs(options.infilename,
                    options.rrbs_taginfo,
                    options.adapter_file,
                    options.cutnumber1,
                    options.cutnumber2,
                    options.no_split,
                    indexname,
                    aligner_command,
                    db_path,
                    tmp_path,
                    options.outfilename or options.infilename+'.rrbsse')
        else: # Normal single end scan
            bs_single_end(  options.infilename,
                        asktag,
                        options.adapter_file,
                        options.cutnumber1,
                        options.cutnumber2,
                        options.no_split,
                        indexname,
                        aligner_command,
                        db_path,
                        tmp_path,
                        options.outfilename or options.infilename+'.bsse' # this is the output file name
                        )
    else:
        print 'Pair end'
        # pair end specific default options
        aligner_options = dict({BOWTIE: {'--ff'  : asktag == 'N',
                                         '--fr'  : asktag == 'Y',
                                         '-X'    : options.max_insert_size,
                                         '-I'    : options.min_insert_size if options.min_insert_size > 0 else None
                                },
                                BOWTIE2 : {
                                         '--ff'  : asktag == 'N',
                                         '--fr'  : asktag == 'Y',
                                         '-X'    : options.max_insert_size,
                                         '-I'    : options.min_insert_size if options.min_insert_size > 0 else None
                                },
                                SOAP: {}}[options.aligner],
                                **aligner_options)



        aligner_command = 'nohup ' + aligner_exec + aligner_options_string() + { BOWTIE   : ' %(reference_genome)s  -f -1 %(input_file_1)s -2 %(input_file_2)s %(output_file)s',
                                                                                 BOWTIE2  : ' -x %(reference_genome)s  -f -1 %(input_file_1)s -2 %(input_file_2)s -S %(output_file)s',
                                                                                 SOAP     : ' '}[options.aligner]

        print 'Aligner command:', aligner_command

        bs_pair_end(options.infilename_1,
                    options.infilename_2,
                    asktag,
                    options.adapter_file,
                    options.cutnumber1,
                    options.cutnumber2,
                    options.no_split,
                    indexname,
                    aligner_command,
                    db_path,
                    tmp_path,
                    options.outfilename or options.infilename_1+'.bspe' # this is the output file name
             )

