import fileinput, string,os, gzip,copy, time, subprocess, random, math
import json
from subprocess import Popen
from utils import *


#----------------------------------------------------------------
def extract_mapping(ali_file):
    U={}
    R={}
    header0=""
    lst=[]
    for line in fileinput.input(tmp_path+ali_file):
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


from optparse import OptionParser


#noinspection PyUnboundLocalVariable,PyUnboundLocalVariable
if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-i", "--input", type="string", dest="infilename",help="Input your read file name (FORMAT: sequences, illumina fastq, qseq,fasta)", metavar="FILE")

    parser.set_defaults(taginfo="N")
    parser.add_option("-t", "--tag", type="string", dest="taginfo",help="Yes for undirectional lib, no for directional [N]", metavar="TAG")

    parser.set_defaults(cutnumber1=1)
    parser.add_option("-s","--start_base",type = "int",dest = "cutnumber1", help="The first base of your read to be mapped [1]")

    parser.set_defaults(cutnumber2=36)
    parser.add_option("-e","--end_base",type = "int",dest = "cutnumber2", help="The last cycle number of your read to be mapped [36]")

    parser.set_defaults(adapterfilename="")
    parser.add_option("-a", "--adapter", type="string", dest="adapterfilename",help="Input text file of your adaptor sequences (to be trimed from the 3'end of the reads). Input 1 seq for dir. lib., 2 seqs for undir. lib. One line per sequence", metavar="FILE")


    parser.set_defaults(bowtiepath=default_bowtie_path)
    parser.add_option("-p", "--path", type="string", dest="bowtiepath",help="Path to Bowtie [%s]" % default_bowtie_path, metavar="PATH")

    parser.set_defaults(dbpath = reference_genome_path)
    parser.add_option("-d", "--db", type="string", dest="dbpath",help="Path to the reference genome library (generated in preprocessing genome) [%s]" % reference_genome_path, metavar="DBPATH")

    parser.add_option("-g", "--genome", type="string", dest="genome",help="Name of the reference genome (the same as the reference genome file in the preprocessing step) [ex. chr21_hg18]")

    parser.set_defaults(indexname=3)
    parser.add_option("-m", "--mis",type = "int", dest="indexname",help="Number of mismatches (0,1,...,read length) [3]")

    parser.set_defaults(no_split=4000000)
    parser.add_option("-l", "--split_line",type = "int", dest="no_split",help="Number of lines per split (the read file will be split into small files for mapping. The result will be merged. [4000000]")

    parser.set_defaults(outfilename="BS_SEEKER_SE_OUTPUT.txt")
    parser.add_option("-o", "--output", type="string", dest="outfilename",help="The name of output file [BS_SEEKER_SE_OUTPUT.txt]", metavar="OUTFILE")
    #----------------------------------------------------------------
    (options, args) = parser.parse_args()


    # if no options were given by the user, print help and exit
    import sys
    if len(sys.argv) == 1:
        print parser.print_help()
        exit(0)


    main_read_file=options.infilename

    asktag=str(options.taginfo).upper()
    if asktag not in 'YN':
        error('-t option should be either Y or N, not %s' % asktag)


    adapter_file=options.adapterfilename

    cut1=options.cutnumber1
    cut2=options.cutnumber2

    no_small_lines=options.no_split

    indexname=options.indexname
    int_no_mismatches=min(int(indexname),cut2)
    indexname=str(int_no_mismatches)

    bowtie_path=options.bowtiepath
    if bowtie_path[-1] !="/":
        bowtie_path += "/"

    genome = options.genome
    if genome is None:
        error('-g is a required option')

    db_path=options.dbpath
    if db_path[-1] !="/":
        db_path += "/"

    db_path += genome +  '_' + asktag + '_'

    if not os.path.isfile(db_path+'ref.json'):
        print db_path+'ref.json'
        error(genome + ' cannot be found in ' + options.dbpath +'. Please, run the Preprocessing_genome to create it.')

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

    outf = open(options.outfilename,'w')
    logoutf = open(options.outfilename + '.log_BS_Seeker_SE','w')

    #----------------------------------------------------------------
    logoutf.write("Read filename: %s"% main_read_file +"\n")
    logoutf.write("Undirectional library: %s" % asktag + "\n")
    logoutf.write("The first base (for mapping): %d" % cut1 +"\n")
    logoutf.write("The last base (for mapping): %d" % cut2 + "\n")
    logoutf.write("Max. lines per mapping: %d"% no_small_lines +"\n")
    logoutf.write("Bowtie path: %s" % bowtie_path +"\n")
    logoutf.write("Reference genome library path: %s" % db_path + "\n")
    logoutf.write("Number of mismatches allowed: %s" % indexname + "\n")
    if adapter_file !="":
        if asktag=="N":
            logoutf.write("Adapter to be removed from 3'read: %s"%(adapter.rstrip("\n"))+"\n")
        elif asktag=="Y":
            logoutf.write("Adapter to be removed from 3' FW reads: %s"%(adapter_fw.rstrip("\n"))+"\n")
            logoutf.write("Adapter to be removed from 3' RC reads: %s"%(adapter_rc.rstrip("\n"))+"\n")
    #----------------------------------------------------------------

    tmp_path = main_read_file + '-TMP'
    if not os.path.exists(tmp_path):
        os.mkdir(tmp_path)
    tmp_path += '/'

#    tmp_path='./%s-TMP/' % main_read_file
    #----------------------------------------------------------------
    # splitting the big read file

    input_fname = os.path.split(main_read_file)[1]
    Popen('split -l %d %s %s%s-s-' % (no_small_lines, main_read_file, tmp_path, input_fname),shell=True).wait()

    my_files = sorted(splitted_file for splitted_file in os.listdir(tmp_path)
                                            if splitted_file.startswith("%s-s-" % input_fname))

    #--- Reference genome -------------------------------------------------------------
    print "== Reading reference genome =="

    genome_seqs = json.load(open(db_path+"ref.json"))

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

    for read_file in my_files:
        no_my_files+=1
        random_id = ".tmp-"+str(random.randint(1000000,9999999))
        if asktag=="Y":

            #----------------------------------------------------------------
            outfile1='Trimed_BS.fa'+random_id
            outfile2='Trimed_C2T.fa'+random_id
            outfile3='Trimed_G2A.fa'+random_id

            outf1=open(tmp_path+outfile1,'w')
            outf2=open(tmp_path+outfile2,'w')
            outf3=open(tmp_path+outfile3,'w')

            #----------------------------------------------------------------
            read_inf=open(tmp_path+read_file,"r")
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

            #print "detected data format: %s"%(input_format)
            #----------------------------------------------------------------
            id=""
            seq=""
            seq_ready="N"
            for line in fileinput.input(tmp_path+read_file):
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
                        seq=''.join(["N" for x in range(cut2-cut1+1)])


                    #---------  trimmed_raw_BS_read  ------------------
                    outf1.write('>%s'%id+"\n")
                    outf1.write('%s'%seq+"\n")

                    #---------  FW_C2T  ------------------
                    FWseq=copy.deepcopy(seq)
                    FWseq=FWseq.replace("C","T")
                    outf2.write('>%s'%id+"\n")
                    outf2.write('%s'%FWseq+"\n")
                    #---------  RC_G2A  ------------------
                    RC_seq=copy.deepcopy(seq)
                    RC_seq=RC_seq.replace("G","A")
                    outf3.write('>%s'%id+"\n")
                    outf3.write('%s'%RC_seq+"\n")


            fileinput.close()

            outf1.close()
            outf2.close()
            outf3.close()

            Popen('rm %s &'%(tmp_path+read_file),shell=True)

           #--------------------------------------------------------------------------------
            # Bowtie mapping
            #-------------------------------------------------------------------------------
            WC2T="W_C2T_m"+str(indexname)+".mapping"+random_id
            CC2T="C_C2T_m"+str(indexname)+".mapping"+random_id
            WG2A="W_G2A_m"+str(indexname)+".mapping"+random_id
            CG2A="C_G2A_m"+str(indexname)+".mapping"+random_id

            bowtie_map1=Popen('%sbowtie -e %d --nomaqround --norc -k 2 --quiet --best --suppress 2,5,6 -p 2 %sW_C2T  -f %s%s %s%s '%(bowtie_path,40*int_no_mismatches,db_path,tmp_path,outfile2,tmp_path,WC2T),shell=True)
            bowtie_map2=Popen('%sbowtie -e %d --nomaqround --norc -k 2 --quiet --best --suppress 2,5,6 -p 2 %sC_C2T  -f %s%s %s%s '%(bowtie_path,40*int_no_mismatches,db_path,tmp_path,outfile2,tmp_path,CC2T),shell=True)
            bowtie_map3=Popen('%sbowtie -e %d --nomaqround --norc -k 2 --quiet --best --suppress 2,5,6 -p 2 %sW_G2A  -f %s%s %s%s '%(bowtie_path,40*int_no_mismatches,db_path,tmp_path,outfile3,tmp_path,WG2A),shell=True)
            bowtie_map4=Popen('%sbowtie -e %d --nomaqround --norc -k 2 --quiet --best --suppress 2,5,6 -p 2 %sC_G2A  -f %s%s %s%s '%(bowtie_path,40*int_no_mismatches,db_path,tmp_path,outfile3,tmp_path,CG2A),shell=True)
            bowtie_map1.wait()
            bowtie_map2.wait()
            bowtie_map3.wait()
            bowtie_map4.wait()

            Popen('rm %s%s %s%s &'%(tmp_path,outfile2, tmp_path,outfile3),shell=True)

            #--------------------------------------------------------------------------------
            # Post processing
            #--------------------------------------------------------------------------------

            ali_path='./'

            ali_file1=WC2T # mapped to Watson strand, has methylation info of Watson strand
            ali_file2=CG2A # mapped to Crick strand, has methylation info of Watson strand
            ali_file3=WG2A # mapped to Crick strand, has methylation info of Crick strand
            ali_file4=CC2T # mapped to Watson strand, has methylation info of Crick strand

            #----------------------------------------------------------------

            FW_C2T_U,FW_C2T_R=extract_mapping(ali_file1)
            RC_G2A_U,RC_G2A_R=extract_mapping(ali_file2)

            FW_G2A_U,FW_G2A_R=extract_mapping(ali_file3)
            RC_C2T_U,RC_C2T_R=extract_mapping(ali_file4)

            #----------------------------------------------------------------
            # get uniq-hit reads
            #----------------------------------------------------------------
            Union_set=set(FW_C2T_U.keys()) | set(RC_G2A_U.keys()) | set(FW_G2A_U.keys()) | set(RC_C2T_U.keys())

            Unique_FW_C2T=set() # +
            Unique_RC_G2A=set() # +
            Unique_FW_G2A=set() # -
            Unique_RC_C2T=set() # -


            for x in Union_set:
                list=[]
                for d in [FW_C2T_U,RC_G2A_U,FW_G2A_U,RC_C2T_U]:
                    mis_lst=d.get(x,[99])
                    mis=int(mis_lst[0])
                    list.append(mis)
                for d in [FW_C2T_R,RC_G2A_R,FW_G2A_R,RC_C2T_R]:
                    mis=d.get(x,99)
                    list.append(mis)
                mini=min(list)
                if list.count(mini)==1:
                    mini_index=list.index(mini)
                    if mini_index==0:
                        Unique_FW_C2T.add(x)
                    elif mini_index==1:
                        Unique_RC_G2A.add(x)
                    elif mini_index==2:
                        Unique_FW_G2A.add(x)
                    elif mini_index==3:
                        Unique_RC_C2T.add(x)

            Union_set=set()

            FW_C2T_R={}
            RC_C2T_R={}
            FW_G2A_R={}
            RC_G2A_R={}

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
            n1=len(Unique_FW_C2T)
            n2=len(Unique_RC_G2A)
            n3=len(Unique_FW_G2A)
            n4=len(Unique_RC_C2T)

            numbers_premapped_lst[0]+=n1
            numbers_premapped_lst[1]+=n2
            numbers_premapped_lst[2]+=n3
            numbers_premapped_lst[3]+=n4

            Unique_FW_C2T=set()
            Unique_RC_C2T=set()
            Unique_FW_G2A=set()
            Unique_RC_G2A=set()

            '''
            n1234=n1+n2+n3+n4
            logoutf.write("# %10d W_C2T bs-unique reads"%(n1)+"\n")
            logoutf.write("# %10d C_G2A bs-unique reads"%(n2)+"\n")
            logoutf.write("# %10d W_G2A bs-unique reads"%(n3)+"\n")
            logoutf.write("# %10d C_C2T bs-unique reads"%(n4)+"\n")
            logoutf.write("# -- %d + strand bs-unique reads"%(n1+n2)+"\n")
            logoutf.write("# -- %d - strand bs-unique reads"%(n3+n4)+"\n")
            logoutf.write("# -- %d bs-unique reads"%(n1234)+"\n")
            '''

            #------ Original BS reads --------------------------------------
            original_bs_reads={}
            bs_reads_file=outfile1

            for line in fileinput.input(tmp_path+bs_reads_file):
                if line[0]==">":
                    header=line[1:-1]
                else:
                    l=line.split()
                    original_bs_reads[header]=l[0]
            fileinput.close()

            Popen('rm %s%s &'%(tmp_path,outfile1),shell=True)

            #----------------------------------------------------------------

            nn=0
            for ali_unique_lst, ali_dic in [(FW_C2T_uniq_lst,FW_C2T_U),(RC_G2A_uniq_lst,RC_G2A_U),(FW_G2A_uniq_lst,FW_G2A_U),(RC_C2T_uniq_lst,RC_C2T_U)]:
                nn+=1
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
            Popen(' rm %s%s %s%s %s%s %s%s &'%(tmp_path,WC2T, tmp_path,WG2A, tmp_path,CC2T, tmp_path,CG2A),shell=True)

        if asktag=="N":
            #----------------------------------------------------------------
            outfile1='Trimed_BS.fa'+random_id
            outfile2='Trimed_C2T.fa'+random_id

            outf1=open(tmp_path+outfile1,'w')
            outf2=open(tmp_path+outfile2,'w')

            n=0
            #----------------------------------------------------------------
            read_inf=open(tmp_path+read_file,"r")
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
            for line in fileinput.input(tmp_path+read_file):
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
                    outf1.write('>%s'%(id)+"\n")
                    outf1.write('%s'%(seq)+"\n")

                    #---------  FW_C2T  ------------------
                    FWseq=copy.deepcopy(seq)
                    FWseq=FWseq.replace("C","T")
                    outf2.write('>%s'%(id)+"\n")
                    outf2.write('%s'%(FWseq)+"\n")

            fileinput.close()

            outf1.close()
            outf2.close()
            Popen('rm %s &'%(tmp_path+read_file),shell=True)

            #--------------------------------------------------------------------------------
            # Bowtie mapping
            #--------------------------------------------------------------------------------
            WC2T="W_C2T_m"+str(indexname)+".mapping"+random_id
            CC2T="C_C2T_m"+str(indexname)+".mapping"+random_id

            b1='%sbowtie -e %d --nomaqround --norc --best --quiet -k 2 --suppress 2,5,6 -p 3 %sW_C2T -f %s%s %s%s '%(bowtie_path,40*int_no_mismatches,db_path,tmp_path,outfile2,tmp_path,WC2T)
            b2='%sbowtie -e %d --nomaqround --norc --best --quiet -k 2 --suppress 2,5,6 -p 3 %sC_C2T -f %s%s %s%s '%(bowtie_path,40*int_no_mismatches,db_path,tmp_path,outfile2,tmp_path,CC2T)
            #print b1
            #print b2
            bowtie_map1=Popen(b1,shell=True)
            bowtie_map2=Popen(b2,shell=True)
            bowtie_map1.wait()
            bowtie_map2.wait()
            Popen(' rm %s%s &'%(tmp_path,outfile2),shell=True)

            #--------------------------------------------------------------------------------
            # Post processing
            #--------------------------------------------------------------------------------

            ali_path='./'

            ali_file1=WC2T # mapped to Watson strand, has methylation info of Watson strand
            ali_file4=CC2T # mapped to Watson strand, has methylation info of Crick strand

            #----------------------------------------------------------------

            FW_C2T_U,FW_C2T_R=extract_mapping(ali_file1)
            RC_C2T_U,RC_C2T_R=extract_mapping(ali_file4)

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

            Union_set=set()

            FW_C2T_R={}
            RC_C2T_R={}

            FW_C2T_uniq_lst=[[FW_C2T_U[u][1],u] for u in Unique_FW_C2T]
            RC_C2T_uniq_lst=[[RC_C2T_U[u][1],u] for u in Unique_RC_C2T]
            FW_C2T_uniq_lst.sort()
            RC_C2T_uniq_lst.sort()
            FW_C2T_uniq_lst=[x[1] for x in FW_C2T_uniq_lst]
            RC_C2T_uniq_lst=[x[1] for x in RC_C2T_uniq_lst]

            #----------------------------------------------------------------

            n1=len(Unique_FW_C2T)
            n2=len(Unique_RC_C2T)

            n12=n1+n2

            numbers_premapped_lst[0]+=n1
            numbers_premapped_lst[1]+=n2


            Unique_FW_C2T=set() # +
            Unique_RC_C2T=set() # -

            '''
            logoutf.write("# -- %d + strand bs-unique reads"%(n1)+"\n")
            logoutf.write("# -- %d - strand bs-unique reads"%(n2)+"\n")
            logoutf.write("# -- %d bs-unique reads"%(n12)+"\n")
            '''


            #------ Original BS reads --------------------------------------
            original_bs_reads={}
            bs_reads_file=outfile1


            for line in fileinput.input(tmp_path+bs_reads_file):
                if line[0]==">":
                    header=line[1:-1]
                else:
                    l=line.split()
                    original_bs_reads[header]=l[0]
            fileinput.close()

            Popen('rm %s%s &'%(tmp_path,outfile1),shell=True)

            #----------------------------------------------------------------

            nn=0
            for ali_unique_lst, ali_dic in [(FW_C2T_uniq_lst,FW_C2T_U),(RC_C2T_uniq_lst,RC_C2T_U)]:
                nn+=1
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
            Popen(' rm %s%s %s%s &'%(tmp_path,WC2T, tmp_path,CC2T),shell=True)


    #----------------------------------------------------------------

    outf.close()
    Popen('rm -rf %s '% tmp_path ,shell=True)

    logoutf.write("Number of raw reads: %d \n"% all_raw_reads)
    if all_raw_reads >0:
        logoutf.write("Number of reads having adapter removed: %d \n" % all_trimed )
        logoutf.write("Number of unique-hits reads for post-filtering: %d\n" % all_mapped)
        if asktag=="Y":
            logoutf.write(" ---- %7d FW reads mapped to Watson strand (before post-filtering)"%(numbers_premapped_lst[0])+"\n")
            logoutf.write(" ---- %7d RC reads mapped to Watson strand (before post-filtering)"%(numbers_premapped_lst[1])+"\n")
            logoutf.write(" ---- %7d FW reads mapped to Crick strand (before post-filtering)"%(numbers_premapped_lst[2])+"\n")
            logoutf.write(" ---- %7d RC reads mapped to Crick strand (before post-filtering)"%(numbers_premapped_lst[3])+"\n")
        elif asktag=="N":
            logoutf.write(" ---- %7d FW reads mapped to Watson strand (before post-filtering)"%(numbers_premapped_lst[0])+"\n")
            logoutf.write(" ---- %7d FW reads mapped to Crick strand (before post-filtering)"%(numbers_premapped_lst[1])+"\n")

        logoutf.write("Post-filtering %d uniqlely aligned reads with mismatches <= %s"%(all_mapped_passed,str(indexname))+"\n")
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
        logoutf.write(" mCG %1.3f%%"%(100*float(mC_lst[0])/n_CG)+'\n')
        logoutf.write(" mCHG %1.3f%%"%(100*float(mC_lst[1])/n_CHG)+'\n')
        logoutf.write(" mCHH %1.3f%%"%(100*float(mC_lst[2])/n_CHH)+'\n')

    logoutf.write("----------------------------------------------"+"\n")
    logoutf.write("------------------- END --------------------"+"\n")
    elapsed("=== END %s ===" % main_read_file)


    logoutf.close()
