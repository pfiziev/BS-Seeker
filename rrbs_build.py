import fileinput
import json
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
    ref_log=open(os.path.join(ref_path,'log'),"w")

    mappable_genome_fn = os.path.join(ref_path, "RRBS_mapable_genome.fa")

    g=""
    header=""
    genome={}
    all_base=0
    n=0
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
                ref_log.write("reference seq: %s (renamed as %s ) %d bp"%(header,short_header,len(g))+"\n")
                genome[short_header]=g
                all_base+=len(g)

                header=l[0][1:]
                n+=1
                short_header=str(n).zfill(4)

                g=""
    fileinput.close()

    ref_log.write("reference seq: %s (renamed as %s ) %d bp"%(header,short_header,len(g))+"\n")
    genome[short_header]=g.upper()
    all_base+=len(g)
    ref_log.write("--- In total %d reference seqs ==> %d bp"%(len(genome),all_base)+"\n")


    #--- positive strand -------------------------------------------------------------
    outf = open(mappable_genome_fn,"w")
    d2   = {}
    d3   = {}
    f3   = open(os.path.join(ref_path, "RRBS_mapable_regions.txt"),"w")

    all_L=0
    all_mappable_length=0
    all_unmappable_length=0

    no_mapable_region=0
    for chr in sorted(genome.keys()):
        chr_regions=[]
        seq=genome[chr].upper()
        L=len(seq)
        CCGG_sites=[]
        CCGG_CCGG=[]

        #-- collect all "CCGG sites ---------------------------------
        i=1
        while i <= L-4:
            if seq[i:i+4]=="CCGG":
                CCGG_sites.append(i)
            i+=1

        #-- find "CCGG" pairs that are within the length of fragment ---------------------------------
        for j in range(len(CCGG_sites)-1):
            DD = (CCGG_sites[j+1]-CCGG_sites[j])-4 # NOT including both CCGG
            if  low_bound <= DD <= up_bound:
                CCGG_CCGG.append([CCGG_sites[j], CCGG_sites[j+1]+3]) # leftmost <--> rightmost
                mapable_seq=seq[CCGG_sites[j]:CCGG_sites[j+1]+4]
                no_mapable_region+=1


                chr_regions.append([CCGG_sites[j],CCGG_sites[j+1]+3,no_mapable_region,mapable_seq])
                # start_position, end_position, serial, sequence
                f3.write("%s\t%d\t%d\t%d\t%s"%(chr,no_mapable_region,CCGG_sites[j],CCGG_sites[j+1]+3,mapable_seq)+"\n")

        d3[chr]=chr_regions

        #-----------------------------------
        map_seq=""
        mappable_length=0
        unmappable_length=0
        m=0
        mark="close"
        while m < L:
            if len(CCGG_CCGG)>0:
                pair=CCGG_CCGG[0]
                p1=pair[0]
                p2=pair[1]
                if p1 <= m < p2 + 1:
                    map_seq+=seq[m]
                    mappable_length+=1
                    mark="open"
                else:
                    if mark=="close":
                        map_seq+="-"
                        unmappable_length+=1
                    elif mark=="open": # the last eligible base
                        _ = CCGG_CCGG.pop(0)
                        if len(CCGG_CCGG)>0:
                            pair=CCGG_CCGG[0]
                            p1=pair[0]
                            p2=pair[1]
                            if  p1 <= m < p2 + 1:
                                map_seq+=seq[m]
                                mappable_length += 1
                                mark="open"
                            else:
                                map_seq+="-"
                                unmappable_length += 1
                                mark="close"
                        else:
                            map_seq+="-"
                            unmappable_length+=1
                            mark="close"
            else:
                if mark=="close":
                    map_seq+="-"
                    unmappable_length+=1
                elif mark=="open":
                    map_seq+="-"
                    unmappable_length+=1
                    mark="close"

            m+=1

        del genome[chr]
        #-----------------------------------
        d2[chr]=map_seq
        outf.write(">%s"% chr +"\n")
        for ii in range(0,L,50):
            y=min(ii+50,L)
            outf.write("%s"%(map_seq[ii:y])+"\n")

        #-----------------------------------
        ref_log.write("# %s : all %d : %d (unmapable -) %d (mapable) (%1.5f)"%(chr,
                                                                               L,
                                                                               unmappable_length,
                                                                               mappable_length,
                                                                               float(mappable_length)/L)+"\n")
        all_L+=L
        all_mappable_length+=mappable_length
        all_unmappable_length+=unmappable_length

    ref_log.write("# total %d chromosome seqs ==> %d : %d (unmapable -) %d (mapable) (%1.5f)" %(len(genome.keys()),
                                                                                                all_L,
                                                                                                all_unmappable_length,
                                                                                                all_mappable_length,
                                                                                                float(all_mappable_length)/all_L)+"\n")
    ref_log.write("#       %d eligible fragments" % no_mapable_region+"\n")

    json.dump(d2, open(os.path.join(ref_path, "RRBS_mapable_genome.json"), 'w'))
    json.dump(d3, open(os.path.join(ref_path, "RRBS_mapable_regions.json"), 'w'))

    outf.close()
    f3.close()

    gzip = Popen('nohup gzip %s &' % mappable_genome_fn, shell=True)


    # Part 2
    #----------------------------------------------------------------
    ref_log.write("\n")
    ref_log.write("----------------         Pre-processing mapable genome         (2-2) ----------------"+"\n")


    #---------------------------------------------------------------
    # 1. First get the complementary genome (also do the reverse)
    # 2. Then do CT and GA conversions
    #---------------------------------------------------------------
    FW_genome={}
    header=""
    g=''
    n=0

    refd = {}

    for line in fileinput.input(mappable_genome_fn):
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
                #print "reference seq: %s (renamed as %s ) %d bp"%(header,short_header,len(g));
                ref_log.write("Pre-processing reference seq: %s ( %d bp)"%(short_header,len(g))+"\n")
                refd[short_header]=[header,len(g)]
                FW_genome[short_header]=g

                g=""
                header=l[0][1:]
                n+=1
                short_header=str(n).zfill(4)


    short_header=str(n).zfill(4)

    ref_log.write("Pre-processing reference seq: %s ( %d bp)"%(short_header,len(g))+"\n")
    refd[short_header]=[header,len(g)]
    FW_genome[short_header]=g.upper()

    json.dump(refd, open(os.path.join(ref_path,"refname.json"), 'w'))

    FW_lst=sorted(FW_genome.iterkeys())

    json.dump(FW_genome ,open(os.path.join(ref_path,"ref.json"),'w'))


    #---------------- Reverse complement (Crick strand) ----------------------------

    RC_genome = dict((header, reverse_compl_seq(FW_genome[header])) for header in FW_lst)
    RC_lst=sorted(RC_genome.iterkeys())


    #---------------- 4 converted fasta -------------------------------------------

    outf=open(os.path.join(ref_path,'W_C2T.fa'),'w')
    for header in FW_lst:
        outf.write('>%s'%header+"\n")
        g=FW_genome[header]
        g=g.replace("c","t")
        g=g.replace("C","T")
        outf.write('%s'%g+'\n')
    outf.close()
    print 'end 4-1'

    outf=open(os.path.join(ref_path,'C_C2T.fa'),'w')
    for header in RC_lst:
        outf.write('>%s' % header+"\n")
        g=RC_genome[header]
        g=g.replace("c","t")
        g=g.replace("C","T")
        outf.write('%s' % g +'\n')
    outf.close()
    print 'end 4-2'

    outf=open(os.path.join(ref_path,'W_G2A.fa'),'w')
    for header in FW_lst:
        outf.write('>%s'%header+"\n")
        g=FW_genome[header]
        g=g.replace("g","a")
        g=g.replace("G","A")
        outf.write('%s'%g+'\n')
    outf.close()
    print 'end 4-3'

    outf=open(os.path.join(ref_path,'C_G2A.fa'),'w')
    for header in RC_lst:
        outf.write('>%s'%header+"\n")
        g=RC_genome[header]
        g=g.replace("g","a")
        g=g.replace("G","A")
        outf.write('%s'%g+'\n')
    outf.close()
    print 'end 4-4'

    #---------------- bowtie-build -------------------------------------------

    # append ref_path to all elements of to_bowtie
    to_bowtie = map(lambda f: os.path.join(ref_path, f), ['W_C2T', 'W_G2A', 'C_C2T', 'C_G2A'])

    for proc in [Popen( build_command % { 'fname' : fname}, shell=True) for fname in to_bowtie]:
        proc.wait()

    # delete all fasta files
    delete_files(f+'.fa' for f in to_bowtie)

    ref_log.close()
    gzip.wait()
    elapsed('END')
