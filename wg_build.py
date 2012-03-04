import fileinput
import json
from subprocess import Popen
from utils import *


def wg_build(fasta_file, asktag, build_command, ref_path, aligner):

    # ref_path is a string that containts the directory where the reference genomes are stored with
    # the input fasta filename appended
    ref_path = os.path.join(ref_path,
                            os.path.split(fasta_file)[1] + '_' + asktag+'_'+aligner)

    clear_dir(ref_path)
    #---------------------------------------------------------------
    # 1. First get the complementary genome (also do the reverse)
    # 2. Then do CT and GA conversions
    #---------------------------------------------------------------
    FW_genome={}
    header=""
    g=''
    n=0

    ref_log=open(os.path.join(ref_path, "log"),"w")

    refd = {}


    for line in fileinput.input(fasta_file):
        if line[0]==">":
            l = line.split()
            if header == "":
                n += 1
                header = l[0][1:]
                short_header = str(n).zfill(4)
            else:
                print "reference seq: %s (renamed as %s ) %d bp" % (header, short_header, len(g))
                ref_log.write("reference seq: %s (renamed as %s ) %d bp\n" % (header, short_header, len(g)))
                refd[short_header] = [header,len(g)]
                FW_genome[short_header] = g

                g=""
                header = l[0][1:]
                n += 1
                short_header = str(n).zfill(4)

        else:
            g += line.strip().upper()


    short_header=str(n).zfill(4)

    print "reference seq: %s (renamed as %s) %d bp"%(header,short_header,len(g))
    ref_log.write("reference seq: %s (renamed as %s ) %d bp\n"%(header,short_header,len(g)))
    refd[short_header] = [header, len(g)]
    FW_genome[short_header] = g

    json.dump(refd, open(os.path.join(ref_path, 'refname.json'), 'w'))

    ref_log.close()

    FW_lst=sorted(FW_genome.keys())

    json.dump(FW_genome, open(os.path.join(ref_path, 'ref.json'), 'w'))

    #---------------- Reverse complement (Crick strand) ----------------------------

    RC_genome = dict((header, reverse_compl_seq(FW_genome[header])) for header in FW_lst)
    RC_lst=sorted(RC_genome.iterkeys())


    if asktag=="Y":
        #---------------- 4 converted fasta -------------------------------------------

        outf=open(os.path.join(ref_path, 'W_C2T.fa'),'w')
        for header in FW_lst:
            outf.write('>%s\n' % header)
            g=FW_genome[header]
            g=g.replace("C","T")
            outf.write('%s\n' % g)
        outf.close()
        print 'end 4-1'

        outf=open(os.path.join(ref_path, 'C_C2T.fa'),'w')
        for header in RC_lst:
            outf.write('>%s\n'% header)
            g=RC_genome[header]
            g=g.replace("C","T")
            outf.write('%s\n'% g)
        outf.close()
        print 'end 4-2'

        outf=open(os.path.join(ref_path, 'W_G2A.fa'),'w')
        for header in FW_lst:
            outf.write('>%s\n'% header)
            g=FW_genome[header]
            g=g.replace("G","A")
            outf.write('%s\n'% g)
        outf.close()
        print 'end 4-3'

        outf=open(os.path.join(ref_path, 'C_G2A.fa'),'w')
        for header in RC_lst:
            outf.write('>%s\n'% header)
            g=RC_genome[header]
            g=g.replace("G","A")
            outf.write('%s\n' % g)
        outf.close()
        print 'end 4-4'
        #---------------- bowtie libraries -------------------------------------------
        to_bowtie = ['W_C2T', 'W_G2A', 'C_C2T', 'C_G2A']

    else: # asktag=="N"
        #---------------- 2 converted fasta -------------------------------------------

        outf=open(os.path.join(ref_path,'W_C2T.fa'),'w')
        for header in FW_lst:
            outf.write('>%s\n' % header)
            g=FW_genome[header]
            g=g.replace("C","T")
            outf.write('%s\n' % g)
        outf.close()
        print 'end 2-1'

        outf=open(os.path.join(ref_path,'C_C2T.fa'),'w')
        for header in RC_lst:
            outf.write('>%s\n'% header)
            g=RC_genome[header]
            g=g.replace("C","T")
            outf.write('%s\n'% g)
        outf.close()
        print 'end 2-2'
        to_bowtie = ['W_C2T', 'C_C2T']


    # append ref_path to all elements of to_bowtie
    to_bowtie = map(lambda f: os.path.join(ref_path, f), to_bowtie)

    # start bowtie-build for all converted genomes and wait for the processes to finish
    for proc in [Popen(build_command % { 'fname' : fname }, shell=True) for fname in to_bowtie]:
        proc.wait()

    # delete fasta files of converted genomes
    delete_files(f+'.fa' for f in to_bowtie)


    elapsed('Done')

