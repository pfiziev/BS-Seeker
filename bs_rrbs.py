import fileinput, copy , random, math, os.path
from subprocess import Popen
from utils import *

from bs_single_end import extract_mapping


def my_mapable_region(chr_regions, mapped_location, FR): # start_position (first C), end_position (last G), serial, sequence
    #print len(chr_regions)
    out_serial=0
    out_start=-1
    out_end=-1
    if FR == "+FW":
        my_location=str(mapped_location-2)
        if my_location in chr_regions:
            my_lst = chr_regions[my_location]
            out_start = int(my_location)
            out_end = my_lst[0]
            out_serial = my_lst[1]
    elif FR == "-FW":
        my_location = str(mapped_location-1)
        if my_location in chr_regions:
            my_lst = chr_regions[my_location]
            out_end = int(my_location)
            out_start = my_lst[0]
            out_serial = my_lst[1]
    return out_serial, out_start, out_end




#----------------------------------------------------------------

def bs_rrbs(main_read_file, mytag, adapter_file, cut1, cut2, no_small_lines, indexname, aligner_command, db_path, tmp_path, outfilename):
    
    mytag_lst = mytag.split("/")
    #----------------------------------------------------------------

    # helper method to join fname with tmp_path
    tmp_d = lambda fname: os.path.join(tmp_path, fname)

    #----------------------------------------------------------------
    adapter_seq=""
    if adapter_file !="":
        adapter_inf=open(adapter_file,"r")
        adapter_seq=adapter_inf.readline()
        adapter_inf.close()
        adapter_seq=adapter_seq.rstrip("\n")

    #----------------------------------------------------------------
    # splitting the big read file

    input_fname = os.path.split(main_read_file)[1]

    split_file(main_read_file, tmp_d(input_fname)+'-s-', no_small_lines)
    my_files = sorted(tmp_d(splitted_file) for splitted_file in os.listdir(tmp_path)
                                if splitted_file.startswith("%s-s-" % input_fname))


    #----------------------------------------------------------------
    # output files

    outfile=outfilename
    outf=open(outfile,'w')


    open_log(outfilename+'.log_RRBS_Seeker_SE')

    #----------------------------------------------------------------
    print "Read filename: %s" % main_read_file
    print "Starting Msp-1 tag: %s"% mytag
    print "The last cycle (for mapping): %d"% cut2
    print "Bowtie path: %s" % aligner_command
    print "Reference genome library path: %s" % db_path
    print "Number of mismatches allowed: %s" % indexname

    logm("I Read filename: %s" % main_read_file)
    logm("I Starting Msp-1 tag: %s" % mytag )
    logm("I The last cycle (for mapping): %d" % cut2 )
    logm("I Bowtie path: %s" % aligner_command )
    logm("I Reference genome library path: %s" % db_path )
    logm("I Number of mismatches allowed: %s" % indexname)
    logm("I adapter seq: %s" % adapter_seq)
    logm("----------------------------------------------")
    #--- Reference genome -------------------------------------------------------------

    print "== Reading reference genome =="

    genome_seqs = deserialize(os.path.join(db_path,"ref.data"))

    logm("G %d ref sequence(s)"%(len(genome_seqs)))
    logm("----------------------------------------------")


    #--- Mappable regions -------------------------------------------------------------
    FW_regions={}
    RC_regions={}
    d2 = deserialize(os.path.join(db_path,"RRBS_mapable_regions.data"))
    n_mapable_regions=0
    for chr in d2:
        FW_temp_regions={}
        RC_temp_regions={}
        for x in d2[chr]:
            FW_temp_regions[str(x[0])]=[x[1],x[2]]
            RC_temp_regions[str(x[1])]=[x[0],x[2]]
        FW_regions[chr]=FW_temp_regions
        RC_regions[chr]=RC_temp_regions
        n_mapable_regions+=len(FW_temp_regions)

    del d2
    logm("G %d mapable fragments" % n_mapable_regions + "\n")
    logm("----------------------------------------------"+"\n")

    #----------------------------------------------------------------
    all_raw_reads=0
    all_tagged=0
    all_tagged_trimed=0
    all_mapped=0
    all_mapped_passed=0
    n_mytag_lst={}
    for x in mytag_lst:
        n_mytag_lst[x]=0

    mC_lst=[0,0,0]
    uC_lst=[0,0,0]

    no_my_files=0
    original_bs_reads = {}
    #----------------------------------------------------------------
    print "== Start mapping =="
    for read_file in my_files:
        logm("Processing read file: %s\n" % read_file)

        no_my_files+=1
        random_id = ".tmp-"+str(random.randint(1000000,9999999))
        outfile2=tmp_d('Trimed_C2T.fa'+random_id)

        outf2=open(outfile2,'w')

        #--- Checking input format ------------------------------------------
        read_inf=open(read_file,"r")
        oneline=read_inf.readline()
        l=oneline.split()
        input_format=""
        if oneline[0]=="@":	# Illumina GAII FastQ (Lister et al Nature 2009)
            input_format="FastQ"
            n_fastq=0
        elif len(l)==1 and oneline[0]!=">": 	# pure sequences
            input_format="list of sequences"
        elif len(l)==11:	# Illumina GAII qseq file
            input_format="Illumina GAII qseq file"
        elif oneline[0]==">":	# fasta
            input_format="fasta"
            n_fasta=0
        read_inf.close()

        #----------------------------------------------------------------
        id=""
        seq=""
        seq_ready=0
        for line in fileinput.input(read_file):
            l=line.split()

            if input_format=="Old Solexa Seq file":
                all_raw_reads+=1
                seq_ready="N"
                id=str(all_raw_reads)
                id=id.zfill(12)
                seq=l[4]
                seq_ready="Y"
            elif input_format=="list of sequences":
                all_raw_reads+=1
                seq_ready="N"
                id=str(all_raw_reads)
                id=id.zfill(12)
                seq=l[0]
                seq_ready="Y"
            elif input_format=="FastQ":
                m_fastq=math.fmod(n_fastq,4)
                n_fastq+=1
                seq_ready="N"
                if m_fastq==0:
                    all_raw_reads+=1
                    id=str(all_raw_reads)
                    id=id.zfill(12)
                    seq=""
                elif m_fastq==1:
                    seq=l[0]
                    seq_ready="Y"
                else:
                    seq=""
            elif input_format=="Illumina GAII qseq file":
                all_raw_reads+=1
                seq_ready="N"
                id=str(all_raw_reads)
                id=id.zfill(12)
                seq=l[8]
                seq_ready="Y"
            elif input_format=="fasta":
                m_fasta=math.fmod(n_fasta,2)
                n_fasta+=1
                seq_ready="N"
                if m_fasta==0:
                    all_raw_reads+=1
                    #id=str(all_raw_reads)
                    id=l[0][1:]
                    seq=""
                elif m_fasta==1:
                    seq=l[0]
                    seq_ready="Y"
                else:
                    seq=""
            #---------------------------------------------------------------
            if seq_ready=="Y":
                seq=seq[0:cut2].upper().replace(".","N") #<----------------------selecting 0..52 from 1..72  -e 52

                #-- Selecting Reads with mytag (i.e., CGG or TGG or CGA) -----------------------
                has_tag="N"
                for i in range(cut2):
                    if seq[i:i+3] in mytag_lst and has_tag=="N":
                        all_tagged+=1
                        n_mytag_lst[seq[i:i+3]]+=1
                        has_tag="Y"
                        seq=seq[i:]

                        #-- Trimming adapter sequence ---
                        if adapter_seq !="":
                            small_adapter=adapter_seq[:6]
                            if small_adapter in seq:
                                adapter_index=seq.index(small_adapter)
                                res_seq=seq[adapter_index:]
                                if res_seq in adapter_seq or adapter_seq in res_seq:
                                    all_tagged_trimed+=1
                                    seq=seq[:adapter_index+1]
                        if len(seq)<=4:
                            seq = "N" * cut2

                        break


                if has_tag=="Y":
                    #---------  trimmed_raw_BS_read and qscore ------------------
                    original_bs_reads[id] = seq

                    #---------  FW_C2T  ------------------
                    outf2.write('>%s\n%s\n'%(id, seq.replace('C', 'T')))

        fileinput.close()

        outf2.close()

        delete_files(read_file)

        #--------------------------------------------------------------------------------
        # Bowtie mapping
        #--------------------------------------------------------------------------------
        WC2T=tmp_d("W_C2T_m"+indexname+".mapping"+random_id)
        CC2T=tmp_d("C_C2T_m"+indexname+".mapping"+random_id)

        for proc in [ Popen(aligner_command % {'reference_genome' : os.path.join(db_path,'W_C2T'),
                                               'input_file' : outfile2,
                                               'output_file' : WC2T} ,shell=True),
                      Popen(aligner_command % {'reference_genome' : os.path.join(db_path,'C_C2T'),
                                               'input_file' : outfile2,
                                               'output_file' : CC2T} ,shell=True)]:
            proc.wait()
        logm("Aligning reads is done")

#        b1='%s -e %d --nomaqround --norc --best --quiet -k 2 --suppress 2,5,6 -p 3 %s -f %s %s '%(bowtie_path,40*int_no_mismatches,os.path.join(db_path, 'W_C2T'),outfile2,WC2T)
#        b2='%s -e %d --nomaqround --norc --best --quiet -k 2 --suppress 2,5,6 -p 3 %s -f %s %s '%(bowtie_path,40*int_no_mismatches,os.path.join(db_path, 'C_C2T'),outfile2,CC2T)
#
#        bowtie_map1=Popen(b1,shell=True)
#        bowtie_map2=Popen(b2,shell=True)
#        bowtie_map1.wait()
#        bowtie_map2.wait()

        delete_files(outfile2)

        #--------------------------------------------------------------------------------
        # Post processing
        #--------------------------------------------------------------------------------

        FW_C2T_U,FW_C2T_R=extract_mapping(WC2T)
        RC_C2T_U,RC_C2T_R=extract_mapping(CC2T)
        logm("Extracting alignments is done")

        #----------------------------------------------------------------
        # get uniq-hit reads
        #----------------------------------------------------------------
        Union_set=set(FW_C2T_U.iterkeys()) | set(RC_C2T_U.iterkeys())

        Unique_FW_C2T=set() # +
        Unique_RC_C2T=set() # -


        for x in Union_set:
            list=[]
            for dx in [FW_C2T_U,RC_C2T_U]:
                mis_lst=dx.get(x,[99])
                mis=int(mis_lst[0])
                list.append(mis)
            for dx in [FW_C2T_R,RC_C2T_R]:
                mis=dx.get(x,99)
                list.append(mis)
            mini=min(list)
            if list.count(mini)==1:
                mini_index=list.index(mini)
                if mini_index==0:
                    Unique_FW_C2T.add(x)
                elif mini_index==1:
                    Unique_RC_C2T.add(x)

        del Union_set
        del FW_C2T_R
        del RC_C2T_R

        FW_uniq_lst=[[FW_C2T_U[u][1],u] for u in Unique_FW_C2T]
        RC_uniq_lst=[[RC_C2T_U[u][1],u] for u in Unique_RC_C2T]
        FW_uniq_lst.sort()
        RC_uniq_lst.sort()
        FW_uniq_lst=[x[1] for x in FW_uniq_lst]
        RC_uniq_lst=[x[1] for x in RC_uniq_lst]

        del Unique_FW_C2T
        del Unique_RC_C2T

        #----------------------------------------------------------------

        nn=0
        for ali in [(FW_uniq_lst,FW_C2T_U),(RC_uniq_lst,RC_C2T_U)]:
            nn += 1
            ali_unique_lst = ali[0]
            ali_dic = ali[1]
            mapped_chr0 = ""
            for header in ali_unique_lst:

                _, mapped_chr, mapped_location, cigar_string = ali_dic[header]

                original_BS=original_bs_reads[header]
                #-------------------------------------
                if mapped_chr != mapped_chr0:
                    FW_chr_regions=FW_regions[mapped_chr]
                    RC_chr_regions=RC_regions[mapped_chr]
                    my_gseq=genome_seqs[mapped_chr]
                    chr_length=len(my_gseq)

                    mapped_chr0=mapped_chr

                #-------------------------------------

                # cigar_string is not None only for aligners that output in SAM format.
                # BS Seeker reconstructs the alignments and handles mismatches accordingly.
                original_BS_length = len(original_BS)

                if cigar_string is not None:
                    r_start, r_end, g_len = get_read_start_end_and_genome_length(cigar_string)
                else:
                    r_start, r_end, g_len = 0, original_BS_length, original_BS_length

                all_mapped+=1

                checking_first_C = False
                if nn == 1: 							# +FW mapped to + strand:
                    FR = "+FW"
                    mapped_location += 1 # 1 based (after plus 1)
                    origin_genome_long = my_gseq[mapped_location - 2 - 1 : mapped_location + g_len + 2 - 1]
                    checking_first_C = (origin_genome_long[1:5] == "CCGG")
                    mapped_strand = "+"
                    origin_genome = origin_genome_long[2:-2]

                elif nn==2: 						# -FW mapped to - strand:
                    mapped_strand = "-"
                    FR = "-FW"
                    mapped_location = chr_length - mapped_location - g_len + 1
                    origin_genome_long = my_gseq[mapped_location - 2 - 1 : mapped_location + g_len + 2 - 1]
                    origin_genome_long = reverse_compl_seq(origin_genome_long)
                    checking_first_C = (origin_genome_long[1:5] == "CCGG")

                    origin_genome = origin_genome_long[2:-2]


                if cigar_string is not None:
                    original_BS = original_BS[r_start : r_end]
                    r_aln, g_aln = cigar_to_alignment(cigar_string, original_BS, origin_genome)
                else:
                    r_aln = original_BS
                    g_aln = origin_genome


                if len(r_aln) == len(g_aln) and checking_first_C:
                    #---------------------------------------------
                    if FR=="+FW":
                        my_region_serial, my_region_start, my_region_end = my_mapable_region(FW_chr_regions, mapped_location, "+FW")
                    elif FR=="-FW":
                        my_region_serial, my_region_start, my_region_end = my_mapable_region(RC_chr_regions, mapped_location + g_len, "-FW")
                    #---------------------------------------------

                    N_mismatch = N_MIS(r_aln, g_aln) #+ original_BS_length - (r_end - r_start) # mismatches in the alignment + soft clipped nucleotides

                    if N_mismatch <= int(indexname) and my_region_serial != 0:
                        all_mapped_passed+=1
                        #---------------------------------------------

                        mapped_location = str(mapped_location).zfill(10)

                        coordinate = mapped_chr + mapped_strand + mapped_location

                        output_genome = origin_genome_long[0:2] + "_" + origin_genome + "_" + origin_genome_long[-2:]

                        methy = methy_seq(r_aln, g_aln + origin_genome_long[-2:])

                        mC_lst, uC_lst = mcounts(methy, mC_lst, uC_lst)

                        #---STEVE FILTER----------------
                        condense_seq = methy.replace('-','')
                        STEVE=0
                        if "ZZZ" in condense_seq:
                            STEVE=1

                        outf.write('%s	%2d	%3s	%s	%s	%s	%s	%d	%d	%d	%d\n' % (header, N_mismatch, FR, coordinate, output_genome, original_BS, methy, my_region_serial, my_region_start, my_region_end, STEVE))


        print "--> %s (%d/%d) "%(read_file,no_my_files,len(my_files));

        delete_files(WC2T, CC2T)

    outf.close()
    delete_files(tmp_path)

    #----------------------------------------------------------------


    logm("O Number of raw reads: %d "% all_raw_reads)
    if all_raw_reads >0:
        logm("O Number of CGG/TGG tagged reads: %d (%1.3f)"%(all_tagged,float(all_tagged)/all_raw_reads))
        for kk in range(len(n_mytag_lst)):
            logm("O Number of raw reads with %s tag: %d (%1.3f)"%(mytag_lst[kk],n_mytag_lst[mytag_lst[kk]],float(n_mytag_lst[mytag_lst[kk]])/all_raw_reads))
        logm("O Number of CGG/TGG reads having adapter removed: %d "%all_tagged_trimed)
        logm("O Number of unique-hits reads for post-filtering: %d"%all_mapped)

        logm("O ------ %d uniqlely aligned reads, passed fragment check, with mismatches <= %s"%(all_mapped_passed, indexname))
        logm("O Mapability= %1.4f%%"%(100*float(all_mapped_passed)/all_raw_reads))

        n_CG=mC_lst[0]+uC_lst[0]
        n_CHG=mC_lst[1]+uC_lst[1]
        n_CHH=mC_lst[2]+uC_lst[2]

        logm("----------------------------------------------")
        logm("M Methylated C in mapped reads ")
        logm("M mCG %1.3f%%"%(100*float(mC_lst[0])/n_CG))
        logm("M mCHG %1.3f%%"%(100*float(mC_lst[1])/n_CHG))
        logm("M mCHH %1.3f%%"%(100*float(mC_lst[2])/n_CHH))

    logm("----------------------------------------------")
    logm("------------------- END ----------------------")
    elapsed(main_read_file)

    close_log()
