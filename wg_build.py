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

    open_log(os.path.join(ref_path, 'log'))

    FW_genome, refd = fasta2dict(fasta_file)
    serialize(refd, os.path.join(ref_path,"refname.data"))
    serialize(FW_genome, os.path.join(ref_path, 'ref.data'))
    elapsed('Storing genome')

    FW_lst = sorted(FW_genome.iterkeys())

    #---------------- Reverse complement (Crick strand) ----------------------------

    outf=open(os.path.join(ref_path, 'W_C2T.fa'),'w')
    for header in FW_lst:
        outf.write('>%s\n%s\n' % (header, FW_genome[header].replace("C","T")))
    outf.close()
    elapsed('Watson strand C->T')

    if asktag == 'Y':
        outf=open(os.path.join(ref_path, 'W_G2A.fa'),'w')
        for header in FW_lst:
            outf.write('>%s\n%s\n' % (header, FW_genome[header].replace("G","A")))
        outf.close()
        elapsed('Watson strand G->A')


    # reverse complement in place to save memory
    RC_genome = FW_genome
    for key in RC_genome:
        RC_genome[key] = reverse_compl_seq(RC_genome[key])
    RC_lst = FW_lst

#    RC_genome = dict((header, reverse_compl_seq(FW_genome[header])) for header in FW_lst)
#    RC_lst = sorted(RC_genome.iterkeys())


    outf=open(os.path.join(ref_path, 'C_C2T.fa'),'w')
    for header in RC_lst:
        outf.write('>%s\n%s\n' % (header, RC_genome[header].replace("C","T")))
    outf.close()
    elapsed('Crick strand C->T')


    if asktag == 'Y':
        outf=open(os.path.join(ref_path, 'C_G2A.fa'),'w')
        for header in RC_lst:
            outf.write('>%s\n%s\n' % (header, RC_genome[header].replace("G","A")))
        outf.close()
        elapsed('Crick strand G->A')




#
#
#    if asktag=="Y":
#        #---------------- 4 converted fasta -------------------------------------------
#
#        outf=open(os.path.join(ref_path, 'W_C2T.fa'),'w')
#        for header in FW_lst:
#            outf.write('>%s\n' % header)
#            g=FW_genome[header]
#            g=g.replace("C","T")
#            outf.write('%s\n' % g)
#        outf.close()
#        print 'end 4-1'
#
#        outf=open(os.path.join(ref_path, 'C_C2T.fa'),'w')
#        for header in RC_lst:
#            outf.write('>%s\n'% header)
#            g=RC_genome[header]
#            g=g.replace("C","T")
#            outf.write('%s\n'% g)
#        outf.close()
#        print 'end 4-2'
#
#        outf=open(os.path.join(ref_path, 'W_G2A.fa'),'w')
#        for header in FW_lst:
#            outf.write('>%s\n'% header)
#            g=FW_genome[header]
#            g=g.replace("G","A")
#            outf.write('%s\n'% g)
#        outf.close()
#        print 'end 4-3'
#
#        outf=open(os.path.join(ref_path, 'C_G2A.fa'),'w')
#        for header in RC_lst:
#            outf.write('>%s\n'% header)
#            g=RC_genome[header]
#            g=g.replace("G","A")
#            outf.write('%s\n' % g)
#        outf.close()
#        print 'end 4-4'
#        #---------------- bowtie libraries -------------------------------------------
#        to_bowtie = ['W_C2T', 'W_G2A', 'C_C2T', 'C_G2A']
#
#    else: # asktag=="N"
#        #---------------- 2 converted fasta -------------------------------------------
#
#        outf=open(os.path.join(ref_path,'W_C2T.fa'),'w')
#        for header in FW_lst:
#            outf.write('>%s\n' % header)
#            g=FW_genome[header]
#            g=g.replace("C","T")
#            outf.write('%s\n' % g)
#        outf.close()
#        print 'end 2-1'
#
#        outf=open(os.path.join(ref_path,'C_C2T.fa'),'w')
#        for header in RC_lst:
#            outf.write('>%s\n'% header)
#            g=RC_genome[header]
#            g=g.replace("C","T")
#            outf.write('%s\n'% g)
#        outf.close()
#        print 'end 2-2'
#        to_bowtie = ['W_C2T', 'C_C2T']
#

    # append ref_path to all elements of to_bowtie
    to_bowtie = map(lambda f: os.path.join(ref_path, f), ['W_C2T', 'C_C2T'] if asktag == 'N' else ['W_C2T', 'W_G2A', 'C_C2T', 'C_G2A'])

    # start bowtie-build for all converted genomes and wait for the processes to finish

    run_in_parallel([build_command % { 'fname' : fname } for fname in to_bowtie])

    # delete fasta files of converted genomes
    delete_files(f+'.fa' for f in to_bowtie)


    elapsed('Done')

