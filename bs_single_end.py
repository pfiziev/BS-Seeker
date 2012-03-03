import fileinput, string,os, gzip,copy, time, subprocess, random, math
import json
from subprocess import Popen
from utils import *


#----------------------------------------------------------------
def extract_mapping(ali_file):
    unique_hits = {}
    non_unique_hits = {}

    header0 = ""
    lst = []

    for header, chr, location, no_mismatch, cigar_string in process_aligner_output(ali_file):

        #------------------------------
        if header != header0:
            #---------- output -----------
            if len(lst) == 1:
                unique_hits[header0] = lst[0]      # [no_mismatch, chr, location]
            elif len(lst) == 2:
                if lst[0][0] < lst[1][0]:
                    unique_hits[header0] = lst[0]
                else:
                    non_unique_hits[header0] = lst[0][0]
            header0 = header
            lst = [(no_mismatch, chr, location, cigar_string)]
        else:
            lst.append((no_mismatch, chr, location, cigar_string))



    if len(lst) == 1:
            unique_hits[header0] = lst[0]      # [no_mismatch, chr, location]
    elif len(lst) == 2:
        if lst[0][0] < lst[1][0]:
            unique_hits[header0] = lst[0]
        else:
            non_unique_hits[header0] = lst[0][0]

#    print "# %s" % ali_file
#    print "# -- %15d unique-hit reads"%(len(unique_hits))
#    print  "# -- %15d multiple-hit reads"%(len(non_unique_hits))
    return unique_hits, non_unique_hits


def bs_single_end(main_read_file, asktag, adapter_file, cut1, cut2, no_small_lines, indexname, aligner_command, db_path, tmp_path, outfilename):
    #----------------------------------------------------------------
    adapter=""
    adapter_fw=""
    adapter_rc=""
    if adapter_file !="":
        adapter_inf=open(adapter_file,"r")
        if asktag=="N": #<--- directional library
            adapter=adapter_inf.readline()
            adapter_inf.close()
            adapter=adapter.rstrip("\n")
        elif asktag=="Y":#<--- undirectional library
            adapter_fw=adapter_inf.readline()
            adapter_rc=adapter_inf.readline()
            adapter_inf.close()
            adapter_fw=adapter_fw.rstrip("\n")
            adapter_rc=adapter_rc.rstrip("\n")
    #----------------------------------------------------------------

    outf = open(outfilename ,'w')
    logoutf = open(outfilename + '.log_BS_Seeker_SE','w')

    #----------------------------------------------------------------
    logoutf.write("Read filename: %s"% main_read_file +"\n")
    logoutf.write("Output filename: %s"% outfilename +"\n")
    logoutf.write("Undirectional library: %s" % asktag + "\n")
    logoutf.write("The first base (for mapping): %d" % cut1 +"\n")
    logoutf.write("The last base (for mapping): %d" % cut2 + "\n")
    logoutf.write("Max. lines per mapping: %d"% no_small_lines +"\n")
    logoutf.write("Aligner: %s" % aligner_command +"\n")
    logoutf.write("Reference genome library path: %s" % db_path + "\n")
    logoutf.write("Number of mismatches allowed: %s" % indexname + "\n")
    if adapter_file !="":
        if asktag=="N":
            logoutf.write("Adapter to be removed from 3' reads: %s"%(adapter.rstrip("\n"))+"\n")
        elif asktag=="Y":
            logoutf.write("Adapter to be removed from 3' FW reads: %s"%(adapter_fw.rstrip("\n"))+"\n")
            logoutf.write("Adapter to be removed from 3' RC reads: %s"%(adapter_rc.rstrip("\n"))+"\n")
    #----------------------------------------------------------------

    # helper method to join fname with tmp_path
    tmp_d = lambda fname: os.path.join(tmp_path, fname)

    #----------------------------------------------------------------
    # splitting the big read file

    input_fname = os.path.split(main_read_file)[1]

    split_file(main_read_file, tmp_d(input_fname)+'-s-', no_small_lines)
    my_files = sorted(splitted_file for splitted_file in os.listdir(tmp_path)
                                            if splitted_file.startswith("%s-s-" % input_fname))
    #--- Reference genome -------------------------------------------------------------
    print "== Reading reference genome =="

    genome_seqs = json.load(open(os.path.join(db_path,"ref.json")))

    logoutf.write("%d ref sequence(s)"%(len(genome_seqs))+"\n")
    logoutf.write("----------------------------------------------"+"\n")

    #---- Stats ------------------------------------------------------------
    all_raw_reads=0
    all_trimed=0
    all_mapped=0
    all_mapped_passed=0

    numbers_premapped_lst=[0,0,0,0]
    numbers_mapped_lst=[0,0,0,0]

    mC_lst=[0,0,0]
    uC_lst=[0,0,0]


    no_my_files=0

    #----------------------------------------------------------------
    print "== Start mapping =="
    original_bs_reads = {}

    for read_file in my_files:
        no_my_files+=1
        random_id = ".tmp-"+str(random.randint(1000000,9999999))
        if asktag=="Y":

            #----------------------------------------------------------------
            outfile2=tmp_d('Trimed_C2T.fa'+random_id)
            outfile3=tmp_d('Trimed_G2A.fa'+random_id)

            outf2=open(outfile2,'w')
            outf3=open(outfile3,'w')

            #----------------------------------------------------------------
            read_inf=open(tmp_d(read_file),"r")
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
            seq_ready="N"
            for line in fileinput.input(tmp_d(read_file)):
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

                #----------------------------------------------------------------
                if seq_ready=="Y":
                    #seq=seq[cut1-1:cut2]  #<----------------------selecting 6..52 from 1..72  -s 6 -e 52
                    seq=seq[cut1-1:cut2] #<----------------------selecting 0..52 from 1..72  -e 52
                    seq=seq.upper()
                    seq=seq.replace(".","N")


                    #--striping BS adapter from 3' read --------------------------------------------------------------
                    if (adapter_fw !="") and (adapter_rc !=""):
                        signature=adapter_fw[:6]
                        if signature in seq:
                            signature_pos=seq.index(signature)
                            if seq[signature_pos:] in adapter_fw:
                                seq=seq[:signature_pos]#+"".join(["N" for x in range(len(seq)-len(signature_pos))])
                                all_trimed+=1
                        else:
                            signature=adapter_rc[:6]
                            if signature in seq:
                                #print id,seq,signature
                                signature_pos=seq.index(signature)
                                if seq[signature_pos:] in adapter_rc:
                                    seq=seq[:signature_pos]#+"".join(["N" for x in range(len(seq)-len(signature_pos))])
                                    all_trimed+=1

                    if len(seq)<=4:
                        seq=''.join(["N" for x in xrange(cut2-cut1+1)])


                    #---------  trimmed_raw_BS_read  ------------------
                    original_bs_reads[id] = seq

                    #---------  FW_C2T  ------------------
                    outf2.write('>%s\n' % id)
                    outf2.write('%s\n' % seq.replace("C","T"))
                    #---------  RC_G2A  ------------------
                    outf3.write('>%s\n' % id)
                    outf3.write('%s\n' % seq.replace("G","A"))


            fileinput.close()

            outf2.close()
            outf3.close()

            os.remove(tmp_d(read_file))

           #--------------------------------------------------------------------------------
            # Bowtie mapping
            #-------------------------------------------------------------------------------
            WC2T=tmp_d("W_C2T_m"+indexname+".mapping"+random_id)
            CC2T=tmp_d("C_C2T_m"+indexname+".mapping"+random_id)
            WG2A=tmp_d("W_G2A_m"+indexname+".mapping"+random_id)
            CG2A=tmp_d("C_G2A_m"+indexname+".mapping"+random_id)

#            print aligner_command % {'int_no_mismatches' : int_no_mismatches,
#                                     'reference_genome' : os.path.join(db_path,'W_C2T'),
#                                     'input_file' : outfile2,
#                                     'output_file' : WC2T}

            for proc in [ Popen(aligner_command % {'reference_genome' : os.path.join(db_path,'W_C2T'),
                                                   'input_file' : outfile2,
                                                   'output_file' : WC2T} ,shell=True),

                          Popen(aligner_command % {'reference_genome' : os.path.join(db_path,'C_C2T'),
                                                   'input_file' : outfile2,
                                                   'output_file' : CC2T} ,shell=True),

                          Popen(aligner_command % {'reference_genome' : os.path.join(db_path,'W_G2A'),
                                                   'input_file' : outfile3,
                                                   'output_file' : WG2A} ,shell=True),

                          Popen(aligner_command % {'reference_genome' : os.path.join(db_path,'C_G2A'),
                                                   'input_file' : outfile3,
                                                   'output_file' : CG2A} ,shell=True)]:
                proc.wait()



#            bowtie_map1=Popen('%s -e %d --nomaqround --norc -k 2 --quiet --best --suppress 2,5,6 -p 2 %s  -f %s %s ' % (bowtie_path,40*int_no_mismatches,os.path.join(db_path,'W_C2T'),outfile2,WC2T),shell=True)
#            bowtie_map2=Popen('%s -e %d --nomaqround --norc -k 2 --quiet --best --suppress 2,5,6 -p 2 %s  -f %s %s ' % (bowtie_path,40*int_no_mismatches,os.path.join(db_path,'C_C2T'),outfile2,CC2T),shell=True)
#            bowtie_map3=Popen('%s -e %d --nomaqround --norc -k 2 --quiet --best --suppress 2,5,6 -p 2 %s  -f %s %s ' % (bowtie_path,40*int_no_mismatches,os.path.join(db_path,'W_G2A'),outfile3,WG2A),shell=True)
#            bowtie_map4=Popen('%s -e %d --nomaqround --norc -k 2 --quiet --best --suppress 2,5,6 -p 2 %s  -f %s %s ' % (bowtie_path,40*int_no_mismatches,os.path.join(db_path,'C_G2A'),outfile3,CG2A),shell=True)
#            bowtie_map1.wait()
#            bowtie_map2.wait()
#            bowtie_map3.wait()
#            bowtie_map4.wait()
#            error('aaa')
            delete_files(outfile2, outfile3)


            #--------------------------------------------------------------------------------
            # Post processing
            #--------------------------------------------------------------------------------

            FW_C2T_U,FW_C2T_R=extract_mapping(WC2T)
            RC_G2A_U,RC_G2A_R=extract_mapping(CG2A)

            FW_G2A_U,FW_G2A_R=extract_mapping(WG2A)
            RC_C2T_U,RC_C2T_R=extract_mapping(CC2T)

            #----------------------------------------------------------------
            # get uniq-hit reads
            #----------------------------------------------------------------
            Union_set=set(FW_C2T_U.iterkeys()) | set(RC_G2A_U.iterkeys()) | set(FW_G2A_U.iterkeys()) | set(RC_C2T_U.iterkeys())

            Unique_FW_C2T=set() # +
            Unique_RC_G2A=set() # +
            Unique_FW_G2A=set() # -
            Unique_RC_C2T=set() # -


            for x in Union_set:
                list=[]
                for d in [FW_C2T_U, RC_G2A_U, FW_G2A_U, RC_C2T_U]:
                    mis_lst=d.get(x,[99])
                    mis=int(mis_lst[0])
                    list.append(mis)
                for d in [FW_C2T_R, RC_G2A_R, FW_G2A_R, RC_C2T_R]:
                    mis=d.get(x,99)
                    list.append(mis)
                mini=min(list)
                if list.count(mini) == 1:
                    mini_index=list.index(mini)
                    if mini_index == 0:
                        Unique_FW_C2T.add(x)
                    elif mini_index == 1:
                        Unique_RC_G2A.add(x)
                    elif mini_index == 2:
                        Unique_FW_G2A.add(x)
                    elif mini_index == 3:
                        Unique_RC_C2T.add(x)


            FW_C2T_uniq_lst=[[FW_C2T_U[u][1],u] for u in Unique_FW_C2T]
            FW_G2A_uniq_lst=[[FW_G2A_U[u][1],u] for u in Unique_FW_G2A]
            RC_C2T_uniq_lst=[[RC_C2T_U[u][1],u] for u in Unique_RC_C2T]
            RC_G2A_uniq_lst=[[RC_G2A_U[u][1],u] for u in Unique_RC_G2A]
            FW_C2T_uniq_lst.sort()
            RC_C2T_uniq_lst.sort()
            FW_G2A_uniq_lst.sort()
            RC_G2A_uniq_lst.sort()
            FW_C2T_uniq_lst=[x[1] for x in FW_C2T_uniq_lst]
            RC_C2T_uniq_lst=[x[1] for x in RC_C2T_uniq_lst]
            FW_G2A_uniq_lst=[x[1] for x in FW_G2A_uniq_lst]
            RC_G2A_uniq_lst=[x[1] for x in RC_G2A_uniq_lst]

            #----------------------------------------------------------------
            numbers_premapped_lst[0] += len(Unique_FW_C2T)
            numbers_premapped_lst[1] += len(Unique_RC_G2A)
            numbers_premapped_lst[2] += len(Unique_FW_G2A)
            numbers_premapped_lst[3] += len(Unique_RC_C2T)


            #----------------------------------------------------------------

            nn=0
            for ali_unique_lst, ali_dic in [(FW_C2T_uniq_lst,FW_C2T_U),(RC_G2A_uniq_lst,RC_G2A_U),(FW_G2A_uniq_lst,FW_G2A_U),(RC_C2T_uniq_lst,RC_C2T_U)]:
                nn += 1
                mapped_chr0 = ""
                for header in ali_unique_lst:
                    l=ali_dic[header]
                    mapped_chr = l[1]
                    mapped_location = l[2]
                    original_BS = original_bs_reads[header]
                    #-------------------------------------
                    if mapped_chr != mapped_chr0:
                        my_gseq=genome_seqs[mapped_chr]
                        chr_length=len(my_gseq)
                        mapped_chr0=mapped_chr
                    #-------------------------------------
                    all_mapped+=1

                    if nn==1: 							# +FW mapped to + strand:
                        FR="+FW"
                        mapped_location += 1
                        origin_genome_long=my_gseq[mapped_location-2-1:mapped_location+len(original_BS)+2-1]
                        origin_genome_long=origin_genome_long.upper()
                        mapped_strand="+"
                        origin_genome=origin_genome_long[2:-2]

                    elif nn==2:  						# +RC mapped to + strand:
                        FR="+RC" # RC reads from -RC reflecting the methylation status on Watson strand (+)
                        mapped_location=chr_length-mapped_location-len(original_BS)+1
                        origin_genome_long=my_gseq[mapped_location-2-1:mapped_location+len(original_BS)+2-1]
                        origin_genome_long=origin_genome_long.upper()
                        mapped_strand="+"
                        origin_genome=origin_genome_long[2:-2]
                        original_BS=reverse_compl_seq(original_BS)  # for RC reads

                    elif nn==3:  						# -RC mapped to - strand:
                        mapped_strand="-"
                        FR="-RC" # RC reads from +RC reflecting the methylation status on Crick strand (-)
                        mapped_location += 1
                        origin_genome_long=my_gseq[mapped_location-2-1:mapped_location+len(original_BS)+2-1]
                        origin_genome_long=reverse_compl_seq(origin_genome_long)
                        origin_genome=origin_genome_long[2:-2]
                        original_BS=reverse_compl_seq(original_BS)  # for RC reads

                    elif nn==4: 						# -FW mapped to - strand:
                        mapped_strand="-"
                        FR="-FW"
                        mapped_location=chr_length-mapped_location-len(original_BS)+1
                        origin_genome_long=my_gseq[mapped_location-2-1:mapped_location+len(original_BS)+2-1]
                        origin_genome_long=reverse_compl_seq(origin_genome_long)
                        origin_genome=origin_genome_long[2:-2]

                    if len(original_BS)==len(origin_genome):
                        N_mismatch=N_MIS(original_BS,origin_genome)
                        if N_mismatch <= int(indexname):
                            numbers_mapped_lst[nn-1]+=1
                            all_mapped_passed+=1
                            mapped_location=str(mapped_location)
                            mapped_location=mapped_location.zfill(10)
                            coordinate=mapped_chr+mapped_strand+mapped_location
                            output_genome=origin_genome_long[0:2]+"_"+origin_genome+"_"+origin_genome_long[-2:]
                            methy=methy_seq(original_BS,output_genome)
                            mC_lst,uC_lst=mcounts(methy,mC_lst,uC_lst)

                            #---STEVE FILTER----------------
                            condense_seq=methy.replace('-','')
                            STEVE=0
                            if "ZZZ" in condense_seq:
                                STEVE=1

                            outf.write('%s	%2d	%3s	%s	%s	%s	%s	%d'%(header,N_mismatch,FR,coordinate,output_genome,original_BS,methy,STEVE)+"\n")


            #----------------------------------------------------------------
            print "--> %s (%d/%d) "%(read_file,no_my_files,len(my_files))
            delete_files(WC2T, WG2A, CC2T, CG2A)

        if asktag=="N":
            #----------------------------------------------------------------
            outfile2=tmp_d('Trimed_C2T.fa'+random_id)

            outf2=open(outfile2,'w')

            n=0
            #----------------------------------------------------------------
            read_inf=open(tmp_d(read_file),"r")
            oneline=read_inf.readline()
            l=oneline.split()
            input_format=""
            if oneline[0]=="@":	# Illumina GAII FastQ (Lister et al Nature 2009)
                input_format="Illumina GAII FastQ"
                n_fastq=0
            elif len(l)==1 and oneline[0]!=">": 	# pure sequences
                input_format="list of sequences"
            elif len(l)==11:	# Illumina GAII qseq file
                input_format="Illumina GAII qseq file"
            elif oneline[0]==">":	# fasta
                input_format="fasta"
                n_fasta=0
            read_inf.close()
            #print "detected data format: %s"%(input_format)
            #----------------------------------------------------------------
            id=""
            seq=""
            seq_ready="N"
            for line in fileinput.input(tmp_d(read_file)):
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
                elif input_format=="Illumina GAII FastQ":
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

                #----------------------------------------------------------------
                if seq_ready=="Y":
                    #seq=seq[cut1-1:cut2]  #<----------------------selecting 6..52 from 1..72  -s 6 -e 52
                    seq=seq[cut1-1:cut2] #<----------------------selecting 0..52 from 1..72  -e 52
                    seq=seq.upper()
                    seq=seq.replace(".","N")

                    #--striping adapter from 3' read --------------------------------------------------------------
                    if adapter !="":
                        signature=adapter[:6]
                        if signature in seq:
                            signature_pos=seq.index(signature)
                            if seq[signature_pos:] in adapter:
                                seq=seq[:signature_pos]#+"".join(["N" for x in range(len(seq)-len(signature_pos))])
                                all_trimed+=1
                    if len(seq)<=4:
                        seq=''.join(["N" for x in range(cut2-cut1+1)])

                    #---------  trimmed_raw_BS_read  ------------------
                    original_bs_reads[id] = seq


                    #---------  FW_C2T  ------------------
                    outf2.write('>%s\n' % id )
                    outf2.write('%s\n' % seq.replace("C","T"))

            fileinput.close()

            outf2.close()
            delete_files(tmp_d(read_file))

            #--------------------------------------------------------------------------------
            # Bowtie mapping
            #--------------------------------------------------------------------------------
            WC2T=tmp_d("W_C2T_m"+indexname+".mapping"+random_id)
            CC2T=tmp_d("C_C2T_m"+indexname+".mapping"+random_id)

#            print aligner_command % {'reference_genome' : os.path.join(db_path,'W_C2T'),
#                                     'input_file' : outfile2,
#                                     'output_file' : WC2T}


            for proc in [Popen(aligner_command % {'reference_genome' : os.path.join(db_path,'W_C2T'),
                                                  'input_file' : outfile2,
                                                  'output_file' : WC2T} ,shell=True),
                         Popen(aligner_command % {'reference_genome' : os.path.join(db_path,'C_C2T'),
                                                  'input_file' : outfile2,
                                                  'output_file' : CC2T} ,shell=True)]:
                proc.wait()

#            bowtie_map1=Popen('%s -e %d --nomaqround --norc --best --quiet -k 2 --suppress 2,5,6 -p 3 %s -f %s %s '%(bowtie_path,40*int_no_mismatches,os.path.join(db_path,'W_C2T'),outfile2,WC2T),shell=True)
#            bowtie_map2=Popen('%s -e %d --nomaqround --norc --best --quiet -k 2 --suppress 2,5,6 -p 3 %s -f %s %s '%(bowtie_path,40*int_no_mismatches,os.path.join(db_path,'C_C2T'),outfile2,CC2T),shell=True)

#            bowtie_map1.wait()
#            bowtie_map2.wait()

            os.remove(outfile2)

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
                for d in [FW_C2T_U,RC_C2T_U]:
                    mis_lst=d.get(x,[99])
                    mis=int(mis_lst[0])
                    list.append(mis)
                for d in [FW_C2T_R,RC_C2T_R]:
                    mis=d.get(x,99)
                    list.append(mis)
                mini=min(list)
                if list.count(mini)==1:
                    mini_index=list.index(mini)
                    if mini_index==0:
                        Unique_FW_C2T.add(x)
                    elif mini_index==1:
                        Unique_RC_C2T.add(x)


            FW_C2T_uniq_lst=[[FW_C2T_U[u][1],u] for u in Unique_FW_C2T]
            RC_C2T_uniq_lst=[[RC_C2T_U[u][1],u] for u in Unique_RC_C2T]
            FW_C2T_uniq_lst.sort()
            RC_C2T_uniq_lst.sort()
            FW_C2T_uniq_lst=[x[1] for x in FW_C2T_uniq_lst]
            RC_C2T_uniq_lst=[x[1] for x in RC_C2T_uniq_lst]

            #----------------------------------------------------------------

            n1=len(Unique_FW_C2T)
            n2=len(Unique_RC_C2T)


            numbers_premapped_lst[0]+=n1
            numbers_premapped_lst[1]+=n2


            #----------------------------------------------------------------

            nn=0
            for ali_unique_lst, ali_dic in [(FW_C2T_uniq_lst,FW_C2T_U),(RC_C2T_uniq_lst,RC_C2T_U)]:
                nn+=1
                mapped_chr0=""
                for header in ali_unique_lst:
                    l = ali_dic[header]
                    mapped_chr = l[1]
                    mapped_location = l[2]
                    original_BS = original_bs_reads[header]
                    #-------------------------------------
                    if mapped_chr != mapped_chr0:
                        my_gseq = genome_seqs[mapped_chr]
                        chr_length = len(my_gseq)
                        mapped_chr0 = mapped_chr
                    #-------------------------------------
                    all_mapped+=1
                    if nn==1: 							# +FW mapped to + strand:
                        FR="+FW"
                        mapped_location += 1
                        origin_genome_long=my_gseq[mapped_location-2-1:mapped_location+len(original_BS)+2-1]
                        origin_genome_long=origin_genome_long.upper()
                        mapped_strand="+"
                        origin_genome=origin_genome_long[2:-2]


                    elif nn==2: 						# -FW mapped to - strand:
                        mapped_strand="-"
                        FR="-FW"
                        mapped_location=chr_length-mapped_location-len(original_BS)+1
                        origin_genome_long=my_gseq[mapped_location-2-1:mapped_location+len(original_BS)+2-1]
                        origin_genome_long=reverse_compl_seq(origin_genome_long)
                        origin_genome=origin_genome_long[2:-2]

                    if len(origin_genome)==len(original_BS):
                        N_mismatch=N_MIS(original_BS,origin_genome)
                        if N_mismatch<= int(indexname):
                            numbers_mapped_lst[nn-1]+=1
                            all_mapped_passed+=1
                            mapped_location=str(mapped_location)
                            mapped_location=mapped_location.zfill(10)
                            coordinate=mapped_chr+mapped_strand+mapped_location
                            output_genome=origin_genome_long[0:2]+"_"+origin_genome+"_"+origin_genome_long[-2:]
                            methy=methy_seq(original_BS,output_genome)
                            mC_lst,uC_lst=mcounts(methy,mC_lst,uC_lst)

                            #---STEVE FILTER----------------
                            condense_seq=methy.replace('-','')
                            STEVE=0
                            if "ZZZ" in condense_seq:
                                STEVE=1

                            outf.write('%s	%2d	%3s	%s	%s	%s	%s	%d'%(header,N_mismatch,FR,coordinate,output_genome,original_BS,methy,STEVE)+"\n")


            #----------------------------------------------------------------
            print "--> %s (%d/%d) "%(read_file,no_my_files,len(my_files))
            os.remove(WC2T)
            os.remove(CC2T)


    #----------------------------------------------------------------

    outf.close()
    shutil.rmtree(tmp_path)

    logoutf.write("Number of raw reads: %d \n"% all_raw_reads)
    if all_raw_reads >0:
        logoutf.write("Number of reads having adapter removed: %d \n" % all_trimed )
        logoutf.write("Number of unique-hits reads for post-filtering: %d\n" % all_mapped)
        if asktag=="Y":
            logoutf.write(" ---- %7d FW reads mapped to Watson strand (before post-filtering)"%(numbers_premapped_lst[0])+"\n")
            logoutf.write(" ---- %7d RC reads mapped to Watson strand (before post-filtering)"%(numbers_premapped_lst[1])+"\n")
            logoutf.write(" ---- %7d FW reads mapped to Crick strand (before post-filtering)"%(numbers_premapped_lst[2])+"\n")
            logoutf.write(" ---- %7d RC hreads mapped to Crick strand (before post-filtering)"%(numbers_premapped_lst[3])+"\n")
        elif asktag=="N":
            logoutf.write(" ---- %7d FW reads mapped to Watson strand (before post-filtering)"%(numbers_premapped_lst[0])+"\n")
            logoutf.write(" ---- %7d FW reads mapped to Crick strand (before post-filtering)"%(numbers_premapped_lst[1])+"\n")

        logoutf.write("Post-filtering %d uniqlely aligned reads with mismatches <= %s"%(all_mapped_passed, indexname)+"\n")
        if asktag=="Y":
            logoutf.write(" ---- %7d FW reads mapped to Watson strand"%(numbers_mapped_lst[0])+"\n")
            logoutf.write(" ---- %7d RC reads mapped to Watson strand"%(numbers_mapped_lst[1])+"\n")
            logoutf.write(" ---- %7d FW reads mapped to Crick strand"%(numbers_mapped_lst[2])+"\n")
            logoutf.write(" ---- %7d RC reads mapped to Crick strand"%(numbers_mapped_lst[3])+"\n")
        elif asktag=="N":
            logoutf.write(" ---- %7d FW reads mapped to Watson strand"%(numbers_mapped_lst[0])+"\n")
            logoutf.write(" ---- %7d FW reads mapped to Crick strand"%(numbers_mapped_lst[1])+"\n")
        logoutf.write("Mapability= %1.4f%%"%(100*float(all_mapped_passed)/all_raw_reads)+"\n")

        n_CG=mC_lst[0]+uC_lst[0]
        n_CHG=mC_lst[1]+uC_lst[1]
        n_CHH=mC_lst[2]+uC_lst[2]

        logoutf.write("----------------------------------------------"+"\n")
        logoutf.write("Methylated C in mapped reads "+'\n')

        logoutf.write(" mCG %1.3f%%"%((100*float(mC_lst[0])/n_CG) if n_CG != 0 else 0)+'\n')
        logoutf.write(" mCHG %1.3f%%"%((100*float(mC_lst[1])/n_CHG) if n_CHG != 0 else 0)+'\n')
        logoutf.write(" mCHH %1.3f%%"%((100*float(mC_lst[2])/n_CHH) if n_CHH != 0 else 0)+'\n')

        
    logoutf.write("----------------------------------------------"+"\n")
    logoutf.write("------------------- END --------------------"+"\n")
    elapsed("=== END %s ===" % main_read_file)


    logoutf.close()
