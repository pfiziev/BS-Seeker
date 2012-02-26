import fileinput, copy , random, math, os.path
import json
from subprocess import Popen
from utils import *


#----------------------------------------------------------------
def extract_mapping(ali_file):
    U={}
    R={}
    header0=""
    lst=[]
    for line in fileinput.input(ali_file):
        l=line.split()
        header=l[0]
        chr=str(l[1])
        location=int(l[2])
        #-------- mismatchs -----------
        if len(l)==4:
            no_mismatch=0
        elif len(l)==5:
            no_mismatch=l[4].count(":")
        else:
            print l
        #------------------------------
        if header != header0:
            #---------- output -----------
            if len(lst)==1:
                U[header0]=lst[0]      #[no_mismatch,chr,location]
            elif len(lst)==2:
                if lst[0][0]<lst[1][0]:
                    U[header0]=lst[0]
                else:
                    R[header0]=lst[0][0]
            header0=header
            lst=[[no_mismatch,chr,location]]
        elif header == header0:
            lst.append([no_mismatch,chr,location])
    fileinput.close()

    #logoutf.write("# %s"%(ali_file)+"\n")
    #logoutf.write("# -- %15d unique-hit reads"%(len(U))+"\n")
    #logoutf.write("# -- %15d multiple-hit reads"%(len(R))+"\n")
    return U,R



def my_mapable_region(chr_regions,mapped_location,FR): # start_position (first C), end_position (last G), serial, sequence
    #print len(chr_regions);
    out_serial=0
    out_start=-1
    out_end=-1
    if FR=="+FW":
        my_location=str(mapped_location-2)
        if my_location in chr_regions:
            my_lst=chr_regions[my_location]
            out_start=int(my_location)
            out_end=my_lst[0]
            out_serial=my_lst[1]
    elif FR=="-FW":
        my_location=str(mapped_location-1)
        if my_location in chr_regions:
            my_lst=chr_regions[my_location]
            out_end=int(my_location)
            out_start=my_lst[0]
            out_serial=my_lst[1]
    return out_serial,out_start,out_end



#----------------------------------------------------------------
def methy_seq(r,g_long):
    H=['A','C','T']
    g_long=g_long.replace("_","")
    m_seq=str()
    xx="-"
    for i in range(len(r)):
        if r[i] not in ['C','T']:
            xx="-"
        elif r[i]=="T" and g_long[i+2]=="C": #(unmethylated):
            if g_long[i+3]=="G":
                xx="x"
            elif g_long[i+3] in H :
                if g_long[i+4]=="G":
                    xx="y"
                elif g_long[i+4] in H :
                    xx="z"

        elif r[i]=="C" and g_long[i+2]=="C": #(methylated):
            if g_long[i+3]=="G":
                xx="X"
            elif g_long[i+3] in H :
                if g_long[i+4]=="G":
                    xx="Y"
                elif g_long[i+4] in H :
                    xx="Z"
        else:
            xx="-"

        m_seq += xx
    return m_seq

#----------------------------------------------------------------

def bs_rrbs(main_read_file, mytag, adapter_file, cut1, cut2, no_small_lines, int_no_mismatches, indexname, bowtie_path, db_path, outfilename):
    
    mytag_lst=mytag.split("/")
    #----------------------------------------------------------------
    tmp_path = main_read_file + '-TMP'
    clear_dir(tmp_path)
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
    Popen('split -l %d %s %s-s-'%(no_small_lines, main_read_file, tmp_d(input_fname)),shell=True).wait()
    my_files = sorted(tmp_d(splitted_file) for splitted_file in os.listdir(tmp_path)
                                if splitted_file.startswith("%s-s-" % input_fname))


    #----------------------------------------------------------------
    # output files

    outfile=outfilename
    outf=open(outfile,'w')

    logoutf=open(outfilename+'.log_RRBS_Seeker_SE', 'w')

    #----------------------------------------------------------------
    print "Read filename: %s" % main_read_file
    print "Starting Msp-1 tag: %s"% mytag
    print "The last cycle (for mapping): %d"% cut2
    print "Bowtie path: %s" % bowtie_path
    print "Reference genome library path: %s" % db_path
    print "Number of mismatches allowed: %s" % indexname

    logoutf.write("I Read filename: %s" % main_read_file+"\n")
    logoutf.write("I Starting Msp-1 tag: %s" % mytag + "\n")
    logoutf.write("I The last cycle (for mapping): %d" % cut2 + "\n")
    logoutf.write("I Bowtie path: %s" % bowtie_path + "\n")
    logoutf.write("I Reference genome library path: %s" % db_path + "\n")
    logoutf.write("I Number of mismatches allowed: %s" % indexname + "\n")
    logoutf.write("I adapter seq: %s" % adapter_seq + "\n")
    logoutf.write("----------------------------------------------"+"\n")
    #--- Reference genome -------------------------------------------------------------

    print "== Reading reference genome =="

    genome_seqs = json.load(open(os.path.join(db_path,"ref.json")))

    logoutf.write("G %d ref sequence(s)"%(len(genome_seqs))+"\n")
    logoutf.write("----------------------------------------------"+"\n")


    #--- Mapable regions -------------------------------------------------------------
    FW_regions={}
    RC_regions={}
    d2 = json.load(open(os.path.join(db_path,"RRBS_mapable_regions.json"))) #shelve.open(genome_path+regions_file,'r')
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


    logoutf.write("G %d mapable fragments" % n_mapable_regions + "\n")
    logoutf.write("----------------------------------------------"+"\n")

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

    #----------------------------------------------------------------
    print "== Start mapping =="
    for read_file in my_files:
        no_my_files+=1
        random_id = ".tmp-"+str(random.randint(1000000,9999999))
        outfile1=tmp_d('Trimed_BS.fa'+random_id)
        outfile2=tmp_d('Trimed_C2T.fa'+random_id)

        outf1=open(outfile1,'w')
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
                seq=seq[0:cut2] #<----------------------selecting 0..52 from 1..72  -e 52
                seq=seq.upper()
                seq=seq.replace(".","N")

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
                            seq=''.join(["N" for x in range(cut2)])

                        break


                if has_tag=="Y":
                    #---------  trimmed_raw_BS_read and qscore ------------------
                    outf1.write('>%s'%id+"\n")
                    outf1.write('%s'%seq+"\n")

                    #---------  FW_C2T  ------------------
                    FWseq=copy.deepcopy(seq)
                    FWseq=FWseq.replace("C","T")
                    outf2.write('>%s'%id+"\n")
                    outf2.write('%s'%FWseq+"\n")

        fileinput.close()

        outf1.close()
        outf2.close()

        delete_files(read_file)

        #--------------------------------------------------------------------------------
        # Bowtie mapping
        #--------------------------------------------------------------------------------
        WC2T=tmp_d("W_C2T_m"+str(indexname)+".mapping"+random_id)
        CC2T=tmp_d("C_C2T_m"+str(indexname)+".mapping"+random_id)

        b1='%s -e %d --nomaqround --norc --best --quiet -k 2 --suppress 2,5,6 -p 3 %s -f %s %s '%(bowtie_path,40*int_no_mismatches,os.path.join(db_path, 'W_C2T'),outfile2,WC2T)
        b2='%s -e %d --nomaqround --norc --best --quiet -k 2 --suppress 2,5,6 -p 3 %s -f %s %s '%(bowtie_path,40*int_no_mismatches,os.path.join(db_path, 'C_C2T'),outfile2,CC2T)

        bowtie_map1=Popen(b1,shell=True)
        bowtie_map2=Popen(b2,shell=True)
        bowtie_map1.wait()
        bowtie_map2.wait()

        delete_files(outfile2)

        #--------------------------------------------------------------------------------
        # Post processing
        #--------------------------------------------------------------------------------

        FW_C2T_U,FW_C2T_R=extract_mapping(WC2T)
        RC_C2T_U,RC_C2T_R=extract_mapping(CC2T)

        #----------------------------------------------------------------
        # get uniq-hit reads
        #----------------------------------------------------------------
        Union_set=set(FW_C2T_U.keys()) | set(RC_C2T_U.keys())

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

        Union_set=set()

        FW_C2T_R={}
        RC_C2T_R={}

        FW_uniq_lst=[[FW_C2T_U[u][1],u] for u in Unique_FW_C2T]
        RC_uniq_lst=[[RC_C2T_U[u][1],u] for u in Unique_RC_C2T]
        FW_uniq_lst.sort()
        RC_uniq_lst.sort()
        FW_uniq_lst=[x[1] for x in FW_uniq_lst]
        RC_uniq_lst=[x[1] for x in RC_uniq_lst]
        Unique_FW_C2T=set()
        Unique_RC_C2T=set()
        #------ Original BS reads --------------------------------------
        original_bs_reads={}

        for line in fileinput.input(outfile1):
            if line[0]==">":
                header=line[1:-1]
            else:
                l=line.split()
                original_bs_reads[header]=l[0]

        fileinput.close()

        delete_files(outfile1)

        #----------------------------------------------------------------

        #log0=[n1,n2]
        #log=[0,0]
        nn=0
        for ali in [(FW_uniq_lst,FW_C2T_U),(RC_uniq_lst,RC_C2T_U)]:
            nn+=1
            ali_unique_lst=ali[0]
            ali_dic=ali[1]
            mapped_chr0=""
            for xx in ali_unique_lst:
                l=ali_dic[xx]
                header=xx
                mapped_chr=str(l[1])
                out_chr=mapped_chr

                mapped_location=int(l[2])
                original_BS=original_bs_reads[header]
                #-------------------------------------
                if mapped_chr != mapped_chr0:
                    FW_chr_regions=FW_regions[mapped_chr]
                    RC_chr_regions=RC_regions[mapped_chr]
                    my_gseq=genome_seqs[mapped_chr]
                    chr_length=len(my_gseq)

                    mapped_chr0=mapped_chr

                #-------------------------------------

                all_mapped+=1

                checking_first_C="N"
                if nn==1: 							# +FW mapped to + strand:
                    FR="+FW"
                    mapped_location += 1 # 1 based (after plus 1)
                    origin_genome_long=my_gseq[mapped_location-2-1:mapped_location+len(original_BS)+2-1]
                    if origin_genome_long[1:5]=="CCGG":checking_first_C="Y"
                    origin_genome_long=origin_genome_long.upper()
                    mapped_strand="+"
                    origin_genome=origin_genome_long[2:-2]

                elif nn==2: 						# -FW mapped to - strand:
                    mapped_strand="-"
                    FR="-FW"
                    mapped_location=chr_length-mapped_location-len(original_BS)+1
                    origin_genome_long=my_gseq[mapped_location-2-1:mapped_location+len(original_BS)+2-1]
                    origin_genome_long=reverse_compl_seq(origin_genome_long)
                    if origin_genome_long[1:5]=="CCGG":checking_first_C="Y"
                    origin_genome=origin_genome_long[2:-2]

                if len(origin_genome)==len(original_BS) and checking_first_C=="Y":
                    #---------------------------------------------
                    if FR=="+FW":
                        my_region_serial, my_region_start, my_region_end = my_mapable_region(FW_chr_regions,mapped_location,"+FW")
                    elif FR=="-FW":
                        my_region_serial, my_region_start, my_region_end = my_mapable_region(RC_chr_regions,mapped_location+len(original_BS),"-FW")
                    #---------------------------------------------
                    N_mismatch=N_MIS(original_BS,origin_genome)
                    if N_mismatch <= int(indexname) and my_region_serial != 0:
                        all_mapped_passed+=1
                        #---------------------------------------------

                        mapped_location=str(mapped_location).zfill(10)

                        coordinate=mapped_chr+mapped_strand+mapped_location
                        output_genome=origin_genome_long[0:2]+"_"+origin_genome+"_"+origin_genome_long[-2:]
                        methy=methy_seq(original_BS,output_genome)

                        mC_lst,uC_lst=mcounts(methy,mC_lst,uC_lst)

                        #---STEVE FILTER----------------
                        condense_seq=methy.replace('-','')
                        STEVE=0
                        if "ZZZ" in condense_seq:
                            STEVE=1

                        outf.write('%s	%2d	%3s	%s	%s	%s	%s	%d	%d	%d	%d'%(header,N_mismatch,FR,coordinate,output_genome,original_BS,methy,my_region_serial, my_region_start, my_region_end,STEVE)+"\n")


        print "--> %s (%d/%d) "%(read_file,no_my_files,len(my_files));

        delete_files(WC2T, CC2T)

    outf.close()
    shutil.rmtree(tmp_path)

    #----------------------------------------------------------------


    logoutf.write("O Number of raw reads: %d "% all_raw_reads + "\n")
    if all_raw_reads >0:
        logoutf.write("O Number of CGG/TGG tagged reads: %d (%1.3f)"%(all_tagged,float(all_tagged)/all_raw_reads)+"\n")
        for kk in range(len(n_mytag_lst)):
            logoutf.write("O Number of raw reads with %s tag: %d (%1.3f)"%(mytag_lst[kk],n_mytag_lst[mytag_lst[kk]],float(n_mytag_lst[mytag_lst[kk]])/all_raw_reads)+"\n")
        logoutf.write("O Number of CGG/TGG reads having adapter removed: %d "%all_tagged_trimed+"\n")
        logoutf.write("O Number of unique-hits reads for post-filtering: %d"%all_mapped+"\n")

        logoutf.write("O ------ %d uniqlely aligned reads, passed fragment check, with mismatches <= %s"%(all_mapped_passed,str(indexname))+"\n")
        logoutf.write("O Mapability= %1.4f%%"%(100*float(all_mapped_passed)/all_raw_reads)+"\n")

        n_CG=mC_lst[0]+uC_lst[0]
        n_CHG=mC_lst[1]+uC_lst[1]
        n_CHH=mC_lst[2]+uC_lst[2]

        logoutf.write("----------------------------------------------"+"\n")
        logoutf.write("M Methylated C in mapped reads "+'\n')
        logoutf.write("M mCG %1.3f%%"%(100*float(mC_lst[0])/n_CG)+'\n')
        logoutf.write("M mCHG %1.3f%%"%(100*float(mC_lst[1])/n_CHG)+'\n')
        logoutf.write("M mCHH %1.3f%%"%(100*float(mC_lst[2])/n_CHH)+'\n')

    logoutf.write("----------------------------------------------"+"\n")
    logoutf.write("------------------- END --------------------"+"\n")
    elapsed(main_read_file)

    logoutf.close()
