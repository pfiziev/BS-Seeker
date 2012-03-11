import os
import shutil
from subprocess import Popen
from utils import *

def rrbs_build(fasta_file, asktag, build_command, ref_path, low_bound, up_bound, aligner):

    # ref_path is a string that containts the directory where the reference genomes are stored with
    # the input fasta filename appended
    ref_path = os.path.join(ref_path,
        os.path.split(fasta_file)[1] + '_' + asktag + '_rrbs_%d_%d' % (low_bound, up_bound) +'_' + aligner)

    clear_dir(ref_path)
    open_log(os.path.join(ref_path, 'log'))

    genome, refd = fasta2dict(fasta_file)
    serialize(refd, os.path.join(ref_path,"refname.data"))

    all_base = sum(len(genome[key]) for key in genome)


    logm("--- In total %d reference seqs ==> %d bp"%(len(genome),all_base) )


    #--- positive strand -------------------------------------------------------------
#    outf = open(mappable_genome_fn,"w")
    FW_genome = {}
    d3   = {}
    f3   = open(os.path.join(ref_path, "RRBS_mapable_regions.txt"),"w")

    all_L=0
    all_mappable_length=0
    all_unmappable_length=0

    no_mapable_region=0
    for chr in sorted(genome.keys()):
        chr_regions = []
        seq = genome[chr]
        L = len(seq)
        CCGG_sites = []
        CCGG_CCGG = []

        #-- collect all "CCGG sites ---------------------------------
        i=1
        while i <= L-4:
            if seq[i:i+4]=="CCGG":
                CCGG_sites.append(i)
            i+=1

        #-- find "CCGG" pairs that are within the length of fragment ---------------------------------
        for j in xrange(len(CCGG_sites)-1):
            DD = (CCGG_sites[j+1]-CCGG_sites[j])-4 # NOT including both CCGG
            if  low_bound <= DD <= up_bound:
                CCGG_CCGG.append([CCGG_sites[j], CCGG_sites[j+1]+3]) # leftmost <--> rightmost
                mapable_seq = seq[CCGG_sites[j]:CCGG_sites[j+1]+4]
                no_mapable_region += 1


                chr_regions.append([CCGG_sites[j], CCGG_sites[j+1]+3, no_mapable_region, mapable_seq])
                # start_position, end_position, serial, sequence
                f3.write("%s\t%d\t%d\t%d\t%s\n"%(chr,no_mapable_region,CCGG_sites[j],CCGG_sites[j+1]+3,mapable_seq))

        d3[chr] = chr_regions

        #-----------------------------------
#        map_seq = ""
        _map_seq = []
        mappable_length = 0
        unmappable_length = 0
        m = 0
#        mark = "close"
        mark = False
        while m < L:
            if len(CCGG_CCGG) > 0:
                pair = CCGG_CCGG[0]
                p1 = pair[0]
                p2 = pair[1]
                if p1 <= m < p2 + 1:
#                    map_seq += seq[m]
                    _map_seq.append(seq[m])
                    mappable_length+=1
#                    mark="open"
                    mark = True
                else:
                    if not mark: #mark=="close":
#                        map_seq+="-"
                        _map_seq.append("-")
                        unmappable_length+=1
                    else: # mark=="open": # the last eligible base
                        _ = CCGG_CCGG.pop(0)
                        if len(CCGG_CCGG)>0:
                            pair = CCGG_CCGG[0]
                            p1 = pair[0]
                            p2 = pair[1]
                            if  p1 <= m < p2 + 1:
#                                map_seq += seq[m]
                                _map_seq.append(seq[m])
                                mappable_length += 1
#                                mark="open"
                                mark = True
                            else:
#                                map_seq+="-"
                                _map_seq.append("-")
                                unmappable_length += 1
#                                mark = "close"
                                mark = False
                        else:
#                            map_seq+="-"
                            _map_seq.append("-")
                            unmappable_length+=1
#                            mark = "close"
                            mark = False
            else:
                if not mark: # =="close":
#                    map_seq+="-"
                    _map_seq.append("-")
                    unmappable_length+=1
                else: #if mark=="open":
#                    map_seq+="-"
                    _map_seq.append("-")
                    unmappable_length+=1
#                    mark="close"
                    mark = False

            m+=1

        del genome[chr]
        #-----------------------------------
#        map_seq = ''.join(_map_seq)

        FW_genome[chr] = ''.join(_map_seq)

#        outf.write(">%s\n" % chr)
#        for ii in xrange(0, L, 50):
#            y = min(ii+50, L)
#            outf.write("%s\n" % map_seq[ii:y])

        #-----------------------------------
        logm("# %s : all %d : %d (unmapable -) %d (mapable) (%1.5f)"%(chr,
                                                                               L,
                                                                               unmappable_length,
                                                                               mappable_length,
                                                                               float(mappable_length)/L) )
        all_L += L
        all_mappable_length += mappable_length
        all_unmappable_length += unmappable_length

        elapsed(chr)

    logm("# total %d chromosome seqs ==> %d : %d (unmapable -) %d (mapable) (%1.5f)" %(len(genome.keys()),
                                                                                                all_L,
                                                                                                all_unmappable_length,
                                                                                                all_mappable_length,
                                                                                                float(all_mappable_length)/all_L) )
    logm("#       %d eligible fragments" % no_mapable_region )

    elapsed('Cutting mappable regions')
    serialize(FW_genome, os.path.join(ref_path, "ref.data"))
    serialize(d3, os.path.join(ref_path, "RRBS_mapable_regions.data"))

#    outf.close()
    f3.close()
    elapsed('Store mappable regions and genome')

#    gzip = Popen('nohup gzip %s &' % mappable_genome_fn, shell=True)


    # Part 2
    #----------------------------------------------------------------
    logm("\n")
    logm("----------------         Pre-processing mapable genome         (2-2) ----------------" )


    #---------------------------------------------------------------
    # 1. First get the complementary genome (also do the reverse)
    # 2. Then do CT and GA conversions
    #---------------------------------------------------------------



#    FW_genome = {}
#    header = ""
#    g = ''
#    n = 0
#
#    refd = {}
#
#    for line in fileinput.input(mappable_genome_fn):
#
#        if line[0] == ">":
#            l=line.split()
#            if header=="":
#                n += 1
#                header = l[0][1:]
#                short_header = str(n).zfill(4)
#            else:
#
#                #print "reference seq: %s (renamed as %s ) %d bp"%(header,short_header,len(g));
#                logm("Pre-processing reference seq: %s ( %d bp)\n" % (short_header, len(g)))
#                refd[short_header] = [header, len(g)]
#                FW_genome[short_header] = g
#
#                g = ""
#                header = l[0][1:]
#                n += 1
#                short_header = str(n).zfill(4)
#        else:
#            g += line.strip().upper()
#
#    short_header = str(n).zfill(4)
#
#    logm("Pre-processing reference seq: %s ( %d bp)"%(short_header,len(g)) )
#    refd[short_header] = [header, len(g)]
#    FW_genome[short_header] = g


#    FW_genome, refd = fasta2dict(mappable_genome_fn, header_info = True)


    FW_lst = sorted(FW_genome.iterkeys())

#    serialize(FW_genome ,open(os.path.join(ref_path,"ref.data"),'w'))

    elapsed('Storing forward genome')


    #---------------- 4 converted fasta -------------------------------------------

    outf=open(os.path.join(ref_path,'W_C2T.fa'),'w')
    for header in FW_lst:
        outf.write('>%s\n%s\n' % (header, FW_genome[header].replace("C", "T")))
    outf.close()
    elapsed('end 4-1')

    outf=open(os.path.join(ref_path,'W_G2A.fa'),'w')
    for header in FW_lst:
        outf.write('>%s\n%s\n' % (header, FW_genome[header].replace("G", "A")))
    outf.close()
    elapsed('end 4-3')


    #---------------- Reverse complement (Crick strand) ----------------------------

    # reverse complement in place to save memory
    RC_genome = FW_genome
    for key in RC_genome:
        RC_genome[key] = reverse_compl_seq(RC_genome[header])
    RC_lst=sorted(RC_genome.iterkeys())

    outf=open(os.path.join(ref_path,'C_C2T.fa'),'w')
    for header in RC_lst:
        outf.write('>%s\n%s\n' % (header, RC_genome[header].replace("C", "T")))
    outf.close()
    elapsed('end 4-2')


    outf=open(os.path.join(ref_path,'C_G2A.fa'),'w')
    for header in RC_lst:
        outf.write('>%s\n%s\n' % (header, RC_genome[header].replace("G", "A")))
    outf.close()
    elapsed('end 4-4')

    #---------------- bowtie-build -------------------------------------------

    # append ref_path to all elements of to_bowtie
    to_bowtie = map(lambda f: os.path.join(ref_path, f), ['W_C2T', 'W_G2A', 'C_C2T', 'C_G2A'])

    for proc in [Popen( build_command % { 'fname' : fname}, shell=True) for fname in to_bowtie]:
        proc.wait()
    elapsed('Index building')
    # delete all fasta files
    delete_files(f+'.fa' for f in to_bowtie)

#    ref_log.close()
#    gzip.wait()
    elapsed('END')
