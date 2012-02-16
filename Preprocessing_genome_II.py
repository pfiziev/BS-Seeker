import fileinput, string,os, operator, shelve, time, subprocess
import json
import re
from subprocess import Popen
from optparse import OptionParser
from utils import *


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",help="Input your reference genome file (fasta)", metavar="FILE")

    parser.set_defaults(taginfo="N")
    parser.add_option("-t", "--tag", dest="taginfo",help="Yes for undirectional lib, no for directional [N]", metavar="TAG")

    parser.set_defaults(bowtiepath = default_bowtie_path)
    parser.add_option("-p", "--path", dest="bowtiepath",help="Path to Bowtie [%s]" % default_bowtie_path, metavar="PATH")

    (options, args) = parser.parse_args()

    # if no options were given by the user, print help and exit
    import sys
    if len(sys.argv) == 1:
        print parser.print_help()
        exit(0)



    fasta_file=options.filename
    if fasta_file is None:
        error('Fasta file for the reference genome must be supported')

    if not os.path.isfile(fasta_file):
        error('%s cannot be found' % fasta_file)

    asktag=str(options.taginfo).upper()

    if asktag not in 'YN':
        error('-t option should be either Y or N, not %s' % asktag)


    bowtie_path = options.bowtiepath

    if bowtie_path[-1] != "/":
        bowtie_path += "/"

    print "Reference genome file: %s" % fasta_file
    print "BS reads from undirectional/directional library: %s" % asktag
    print "Bowtie path: %s" % bowtie_path
    #---------------------------------------------------------------

    ref_path = reference_genome_path

    if os.path.exists(ref_path):
        if not os.path.isdir(ref_path):
            error("%s must be a directory" % ref_path)
    else:
        os.mkdir(ref_path)


    # ref_path is a string that containts the directory where the reference genomes are stored with
    # the input fasta filename appended
    ref_path = os.path.join(ref_path,
                            os.path.split(fasta_file)[1] + '_' + asktag + '_')

    #---------------------------------------------------------------
    # 1. First get the complementary genome (also do the reverse)
    # 2. Then do CT and GA conversions
    #---------------------------------------------------------------
    FW_genome={}
    header=""
    g=''
    n=0

    ref_log=open(ref_path[:-1] + ".log","w")

#    refd = shelve.open(ref_path + "refname.shelve",'n')

    refd = {}


    for line in fileinput.input(fasta_file):
        l=line.split()
        if line[0]!=">":
            g=g+line[:-1]
        elif line[0]==">":
            if header=="":
                n+=1
                header=l[0][1:]
                short_header=str(n).zfill(4)
            else:
                g=g.upper()
                print "reference seq: %s (renamed as %s ) %d bp"%(header,short_header,len(g))
                ref_log.write("reference seq: %s (renamed as %s ) %d bp"%(header,short_header,len(g))+"\n")
                refd[short_header]=[header,len(g)]
                FW_genome[short_header]=g

                g=""
                header=l[0][1:]
                n+=1
                short_header=str(n).zfill(4)

    g=g.upper()
    short_header=str(n).zfill(4)
    print "reference seq: %s (renamed as %s) %d bp"%(header,short_header,len(g))
    ref_log.write("reference seq: %s (renamed as %s ) %d bp"%(header,short_header,len(g))+"\n")
    refd[short_header]=[header,len(g)]
    FW_genome[short_header]=g
    g=""

    json.dump(refd, open(ref_path + 'refname.json', 'w'))

#    refd.close()
    ref_log.close()

    FW_lst=FW_genome.keys()
    FW_lst.sort()

    #---------------- Python shelve -----------------------------------------------
    json.dump(FW_genome, open(ref_path + 'ref.json', 'w'))

#    d = shelve.open(ref_path + "ref.shelve",'n')
#    for chr_id in FW_genome:
#        d[chr_id]=FW_genome[chr_id]
#    d.close()

    #---------------- Reverse complement (Crick strand) ----------------------------
    header=""
    RC_genome={}
    for header in FW_lst:
        g=FW_genome[header]
        g=reverse_compl_seq(g)
        RC_genome[header]=g
    RC_lst=RC_genome.keys()
    RC_lst.sort()



    path_dict = {'bowtie_path' : bowtie_path, 'ref_path' : ref_path}

    if asktag=="Y":
        #---------------- 4 converted fasta -------------------------------------------

        outf=open(ref_path + 'W_C2T.fa','w')
        for header in FW_lst:
            outf.write('>%s\n' % header)
            g=FW_genome[header]
            g=g.replace("c","t")
            g=g.replace("C","T")
            outf.write('%s\n' % g)
        outf.close()
        print 'end 4-1'

        outf=open(ref_path + 'C_C2T.fa','w')
        for header in RC_lst:
            outf.write('>%s\n'% header)
            g=RC_genome[header]
            g=g.replace("c","t")
            g=g.replace("C","T")
            outf.write('%s\n'% g)
        outf.close()
        print 'end 4-2'

        outf=open(ref_path + 'W_G2A.fa','w')
        for header in FW_lst:
            outf.write('>%s\n'% header)
            g=FW_genome[header]
            g=g.replace("g","a")
            g=g.replace("G","A")
            outf.write('%s\n'% g)
        outf.close()
        print 'end 4-3'

        outf=open(ref_path + 'C_G2A.fa','w')
        for header in RC_lst:
            outf.write('>%s\n'% header)
            g=RC_genome[header]
            g=g.replace("g","a")
            g=g.replace("G","A")
            outf.write('%s\n' % g)
        outf.close()
        print 'end 4-4'
        #---------------- bowtie libraries -------------------------------------------
        to_bowtie = ['W_C2T', 'W_G2A', 'C_C2T', 'C_G2A']

    else: # asktag=="N"
        #---------------- 2 converted fasta -------------------------------------------

        outf=open(ref_path+'W_C2T.fa','w')
        for header in FW_lst:
            outf.write('>%s\n' % header)
            g=FW_genome[header]
            g=g.replace("c","t")
            g=g.replace("C","T")
            outf.write('%s\n' % g)
        outf.close()
        print 'end 2-1'
        FW_lst={}

        outf=open(ref_path+'C_C2T.fa','w')
        for header in RC_lst:
            outf.write('>%s\n'% header)
            g=RC_genome[header]
            g=g.replace("c","t")
            g=g.replace("C","T")
            outf.write('%s\n'% g)
        outf.close()
        print 'end 2-2'
        to_bowtie = ['W_C2T', 'C_C2T']


    # start bowtie-build for all converted genomes and wait for the processes to finish
    for proc in [Popen('nohup %(bowtie_path)sbowtie-build -f %(fname)s.fa %(fname)s > %(fname)s.log'% {'bowtie_path' : bowtie_path,
                                                                                                       'fname' : ref_path + fname} ,
                       shell=True) for fname in to_bowtie]:
        proc.wait()

    # delete fasta files of converted genomes
    os.system('rm -rf ' + ' '.join(ref_path + fname + '.fa' for fname in to_bowtie))


    elapsed('Done')

