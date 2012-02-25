import fileinput, string,os, gzip,copy, time, shelve,subprocess, random, math
import json
from subprocess import Popen
from utils import *

#----------------------------------------------------------------
def extract_mapping(ali_file):
    U={}
    R={}
    header0=""
    family=[]
    for line in fileinput.input(ali_file):
        l=line.split()
        header=l[0][:-2]
        chr=str(l[1])
        location=int(l[2])
        #no_hits=int(l[4])
        #-------- mismatchs -----------
        if len(l)==4:
            no_mismatch=0
        elif len(l)==5:
            no_mismatch=l[4].count(":")
        else:
            print l
        #------------------------
        if header != header0:
            #--------------------
            if header0 != "":
                # --- output ----
                if len(family)==1:
                    U[header0]=family[0]
                else:
                    if family[0][0]<family[1][0]:
                        U[header0]=family[0]
                    elif family[1][0]<family[0][0]:
                        U[header0]=family[1]
                    else:
                        R[header0]=family[0][0]
                family=[]
                # ---------------
            header0=header
            family=[[no_mismatch,chr,location]]
            member=1
        elif header == header0:
            if member==1:
                family[-1][0]+=no_mismatch
                family[-1].append(location)
                member=2
            elif member==2:
                family.append([no_mismatch,chr,location])
                member=1
        #------------------------------

    fileinput.close()
    return U,R


#----------------------------------------------------------------

from optparse import OptionParser


if __name__ == '__main__':

    parser = OptionParser()
    parser.set_defaults(infilename_1="seq_1.txt")
    parser.add_option("-1", "--input_1", type="string", dest="infilename_1",help="Input your read file end 1 (FORMAT: sequences, illumina fastq, qseq)", metavar="FILE")

    parser.set_defaults(infilename_2="2.txt")
    parser.add_option("-2", "--input_2", type="string", dest="infilename_2",help="Input your read file end 2 (FORMAT: sequences, illumina fastq, qseq)", metavar="FILE")

    parser.set_defaults(taginfo="Y")
    parser.add_option("-t", "--tag", type="string", dest="taginfo",help="Yes for undirectional library, no for directional library [Y]", metavar="TAG")

    parser.set_defaults(cutnumber1=1)
    parser.add_option("-s","--start_base",type = "int",dest = "cutnumber1", help="The first base of your read to be mapped [1]")

    parser.set_defaults(cutnumber2=36)
    parser.add_option("-e","--e",type = "int",dest = "cutnumber2", help="The last cycle number of your read to be mapped [36]")

    parser.set_defaults(min_insert_size=-1)
    parser.add_option("-i","--minins",type = "int",dest = "min_insert_size", help="The minimum insert size for valid paired-end alignments")

    parser.set_defaults(max_insert_size=400)
    parser.add_option("-x","--maxins",type = "int",dest = "max_insert_size", help="The maximum insert size for valid paired-end alignments [400]")

    parser.set_defaults(adapterfilename="")
    parser.add_option("-a", "--adapter", type="string", dest="adapterfilename",help="Input text file of the adapter sequence(s). for directional library 1 seq is allowed. For undirectional library 2 seqs (A then B) are allowed. Input one line per sequence", metavar="FILE")

    parser.set_defaults(bowtiepath=default_bowtie_path)
    parser.add_option("-p", "--path", type="string", dest="bowtiepath",help="Path to Bowtie [%s]" % default_bowtie_path, metavar="PATH")

    parser.set_defaults(dbpath = reference_genome_path)
    parser.add_option("-d", "--db", type="string", dest="dbpath",help="Path to the reference genome library (generated in preprocessing genome) [%s]" % reference_genome_path, metavar="DBPATH")

    parser.add_option("-g", "--genome", type="string", dest="genome",help="Name of the reference genome (the same as the reference genome file in the preprocessing step) [ex. chr21_hg18]")

    parser.set_defaults(indexname=3)
    parser.add_option("-m", "--mis",type = "int", dest="indexname",help="Number of mismatches (0, 1, 2, 3, unlimited< read length) [3]")

    parser.set_defaults(no_split=2000000)
    parser.add_option("-l", "--split_line",type = "int", dest="no_split",help="Number of lines per split (the read file will be split into small files for mapping. The result will be merged.[2000000]")


    parser.add_option("-o", "--output", type="string", dest="outfilename",help="The name of output file [INFILE.bspe]", metavar="OUTFILE")

    #----------------------------------------------------------------

    # parse the options
    (options, args) = parser.parse_args()

    # if no options were given by the user, print help and exit
    import sys
    if len(sys.argv) == 1:
        print parser.print_help()
        exit(0)


    main_read_file_1=options.infilename_1
    main_read_file_2=options.infilename_2



    asktag=options.taginfo

    adapter_file=options.adapterfilename

    min_insert=options.min_insert_size
    max_insert=options.max_insert_size
    max_string="-X %d" % max_insert
    if min_insert >=0:
        min_string="-I %d" % min_insert
    else:
        min_string=""


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

    db_path = os.path.join(options.dbpath, genome + '_' + asktag)

    if not os.path.isfile(os.path.join(db_path,'ref.json')):
        error(genome + ' cannot be found in ' + options.dbpath +'. Please, run the Preprocessing_genome to create it.')



    #----------------------------------------------------------------
    adapter=""
    adapterA=""
    adapterB=""
    if adapter_file !="":
        adapter_inf=open(adapter_file,"r")
        if asktag=="N": #<--- directional library
            adapter=adapter_inf.readline()
            adapter_inf.close()
            adapter=adapter.rstrip("\n")
        elif asktag=="Y":#<--- undirectional library
            adapterA=adapter_inf.readline()
            adapterB=adapter_inf.readline()
            adapter_inf.close()
            adapterA=adapterA.rstrip("\n")
            adapterB=adapterB.rstrip("\n")

    #----------------------------------------------------------------

    outfilename = options.outfilename or main_read_file_1+'.bspe'

    outf = open(outfilename ,'w')
    logoutf = open(outfilename + '.log_BS_Seeker_PE','w')

    outf_u1=open(outfilename+'.Un_E1',"w")
    outf_u2=open(outfilename+'.Un_E2',"w")

    #----------------------------------------------------------------

    logoutf.write("End 1 filename: %s"% main_read_file_1 +'\n')
    logoutf.write("End 2 filename: %s"% main_read_file_2 +'\n')
    logoutf.write("The first base (for mapping): %d"% cut1 +"\n")
    logoutf.write("The last base (for mapping): %d"% cut2 +"\n")

    logoutf.write("-------------------------------- "+'\n')
    logoutf.write("Undirectional library: %s" % asktag +"\n")
    if min_insert >= 0:
        logoutf.write("Min insert size: %d" % min_insert + '\n')
    logoutf.write("Max insert size: %d" % max_insert +'\n')
    logoutf.write("Bowtie path: %s"% bowtie_path + '\n')
    logoutf.write("Reference genome library path: %s"% db_path +'\n')
    logoutf.write("Number of mismatches allowed: %s"% indexname +'\n')

    if adapter_file !="":
        if asktag=="Y":
            logoutf.write("Adapters to be removed from 3' of the reads:"+'\n')
            logoutf.write("-- A: %s" % adapterA +'\n')
            logoutf.write("-- B: %s" % adapterB +'\n')
        elif asktag=="N":
            logoutf.write("Adapter to be removed from 3' of the reads:"+'\n')
            logoutf.write("-- %s" % adapter +'\n')

    logoutf.write("-------------------------------- "+'\n')

    #----------------------------------------------------------------

    tmp_path = main_read_file_1 + '-TMP'
    if not os.path.exists(tmp_path):
        os.mkdir(tmp_path)
    tmp_path += '/'


    #----------------------------------------------------------------
    # splitting the 2 big read files

    input_fname1 = os.path.split(main_read_file_1)[1]
    input_fname2 = os.path.split(main_read_file_2)[1]

    splitting1=Popen('split -l %d %s %s%s-E1-'%(no_small_lines,main_read_file_1,tmp_path,input_fname1),shell=True)
    splitting2=Popen('split -l %d %s %s%s-E2-'%(no_small_lines,main_read_file_2,tmp_path,input_fname2),shell=True)
    splitting1.wait()
    splitting2.wait()

    dirList=os.listdir(tmp_path)
    my_files = zip(sorted(filter(lambda fname: fname.startswith("%s-E1-" % input_fname1), dirList)),
                   sorted(filter(lambda fname: fname.startswith("%s-E2-" % input_fname2), dirList)))


    #--- Reference genome -------------------------------------------------------------
    print "\n"
    print "== Reading reference genome =="
    genome_seqs = json.load(open(os.path.join(db_path,"ref.json")))

    logoutf.write("G %d ref sequence(s)"%(len(genome_seqs))+"\n")

    #---- Stats ------------------------------------------------------------
    all_raw_reads=0
    all_trimed=0
    all_mapped=0
    all_mapped_passed=0

    all_unmapped=0

    all_FW_FR_pairs=0
    all_FW_RF_pairs=0
    all_RC_FR_pairs=0
    all_RC_RF_pairs=0

    numbers_premapped_lst=[0,0,0,0]
    numbers_mapped_lst=[0,0,0,0]

    mC_lst=[0,0,0]
    uC_lst=[0,0,0]

    no_my_files=0
    start = time.time()


    #----------------------------------------------------------------
    print "== Start mapping =="

    for read_file_1, read_file_2 in my_files:

        no_my_files+=1

        random_id=".tmp-"+str(random.randint(1000000,9999999))
        if asktag=="Y":

            #----------------------------------------------------------------
            outfile_1BS = 'Trimed_BS_1.fa'+random_id
            outfile_1FCT= 'Trimed_FCT_1.fa'+random_id
            outfile_1RCT= 'Trimed_RCT_1.fa'+random_id
            outfile_2BS = 'Trimed_BS_2.fa'+random_id
            outfile_2FCT= 'Trimed_FCT_2.fa'+random_id
            outfile_2RCT= 'Trimed_RCT_2.fa'+random_id

            read_inf=open(tmp_path+read_file_1,"r")
            oneline=read_inf.readline()
            l=oneline.split()
            input_format=""

            #if len(l)==5: # old solexa format
            #	input_format="old Solexa Seq file"

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

            print "Detected data format: %s" % input_format

            #----------------------------------------------------------------
            read_file_list=[read_file_1,read_file_2]
            outfile_BS_list=[outfile_1BS,outfile_2BS]
            outfile_FCT_list=[outfile_1FCT,outfile_2FCT]
            outfile_RCT_list=[outfile_1RCT,outfile_2RCT]
            n_list=[0,0]


            for f in range(2):
                read_file=read_file_list[f]
                outf_BS=open(tmp_path+outfile_BS_list[f],'w')
                outf_FCT=open(tmp_path+outfile_FCT_list[f],'w')
                outf_RCT=open(tmp_path+outfile_RCT_list[f],'w')
                n=n_list[f]

                id=""
                seq=""
                seq_ready="N"
                for line in fileinput.input(tmp_path+read_file):
                    l=line.split()
                    if input_format=="old Solexa Seq file":
                        n+=1
                        id=str(n)
                        id=id.zfill(12)
                        seq=l[4]
                        seq_ready="Y"
                    elif input_format=="list of sequences":
                        n+=1
                        seq_ready="N"
                        id=str(n)
                        id=id.zfill(12)
                        seq=l[0]
                        seq_ready="Y"
                    elif input_format=="FastQ":
                        m_fastq=math.fmod(n_fastq,4)
                        n_fastq+=1
                        seq_ready="N"
                        if m_fastq==0:
                            n+=1
                            id=str(n)
                            id=id.zfill(12)
                            seq=""
                        elif m_fastq==1:
                            seq=l[0]
                            seq_ready="Y"
                        else:
                            seq=""
                    elif input_format=="Illumina GAII qseq file":
                        n+=1
                        seq_ready="N"
                        id=str(n)
                        id=id.zfill(12)
                        seq=l[8]
                        seq_ready="Y"
                    elif input_format=="fasta":
                        m_fasta=math.fmod(n_fasta,2)
                        n_fasta+=1
                        seq_ready="N"
                        if m_fasta==0:
                            n+=1

                            #id=str(n)
                            #id=id.zfill(12)
                            id=l[0][1:]
                            id=id.zfill(17)
                            seq=""
                        elif m_fasta==1:
                            seq=l[0]
                            seq_ready="Y"
                        else:
                            seq=""
                    #----------------------------------------------------------------
                    if seq_ready=="Y":
                        seq=seq[cut1-1:cut2] #<----------------------selecting 0..52 from 1..72  -e 52
                        seq=seq.upper()
                        seq=seq.replace(".","N")

                        #--striping BS adapter from 3' read --------------------------------------------------------------
                        if (adapterA !="") and (adapterB !=""):
                            signature=adapterA[:6]
                            if signature in seq:
                                signature_pos=seq.index(signature)
                                if seq[signature_pos:] in adapterA:
                                    seq=seq[:signature_pos]#+"".join(["N" for x in range(len(seq)-len(signature_pos))])
                                    all_trimed+=1
                            else:
                                signature=adapterB[:6]
                                if signature in seq:
                                    #print id,seq,signature;
                                    signature_pos=seq.index(signature)
                                    if seq[signature_pos:] in adapterB:
                                        seq=seq[:signature_pos]#+"".join(["N" for x in range(len(seq)-len(signature_pos))])
                                        all_trimed+=1

                        if len(seq)<=4:
                            seq=''.join(["N" for x in xrange(cut2-cut1+1)])
                        #---------  trimmed_raw_BS_read  ------------------
                        outf_BS.write('%s\t%s'%(id,seq)+"\n")

                        #---------  FW_C2T  ------------------
                        FWseq=copy.deepcopy(seq)
                        FWseq1=FWseq.replace("C","T")
                        outf_FCT.write('>%s' % id +"\n")
                        outf_FCT.write('%s' % FWseq1 +"\n")
                        #---------  RC_G2A  ------------------
                        RCseq=FWseq.replace("G","A")
                        outf_RCT.write('>%s'% id +"\n")
                        outf_RCT.write('%s' % RCseq +"\n")

                n_list[f]=n
                outf_BS.close()
                outf_FCT.close()
                outf_RCT.close()

                fileinput.close()


            #print "All input end 1: %d , end 2: %d "%(n_list[0],n_list[1]);
            all_raw_reads+=n
            elapsed1 = (time.time() - start)

            #--------------------------------------------------------------------------------
            # Bowtie mapping
            #--------------------------------------------------------------------------------
            WC2T_fr="W_C2T_fr_m"+str(indexname)+".mapping"+random_id
            WC2T_rf="W_C2T_rf_m"+str(indexname)+".mapping"+random_id
            CC2T_fr="C_C2T_fr_m"+str(indexname)+".mapping"+random_id
            CC2T_rf="C_C2T_rf_m"+str(indexname)+".mapping"+random_id

            bowtie_map1=Popen('%sbowtie -e %d --nomaqround --fr -k 2 %s %s --quiet --best --suppress 2,5,6 -p 2 %s -f -1 %s%s -2 %s%s %s%s '%(bowtie_path,
                                                                                                                                                   40*int_no_mismatches,
                                                                                                                                                   min_string,
                                                                                                                                                   max_string,
                                                                                                                                                   os.path.join(db_path,'W_C2T'),
                                                                                                                                                   tmp_path,
                                                                                                                                                   outfile_1FCT,tmp_path,outfile_2RCT,tmp_path,WC2T_fr),shell=True)

            bowtie_map2=Popen('%sbowtie -e %d --nomaqround --fr -k 2 %s %s --quiet --best --suppress 2,5,6 -p 2 %s -f -1 %s%s -2 %s%s %s%s '%(bowtie_path,40*int_no_mismatches,min_string, max_string,
                                                                                                                                                   os.path.join(db_path,'C_C2T'),
                                                                                                                                                   tmp_path,outfile_1FCT,tmp_path,outfile_2RCT,tmp_path,CC2T_fr),shell=True)

            bowtie_map3=Popen('%sbowtie -e %d --nomaqround --fr -k 2 %s %s --quiet --best --suppress 2,5,6 -p 2 %s -f -1 %s%s -2 %s%s %s%s '%(bowtie_path,40*int_no_mismatches,min_string, max_string,
                                                                                                                                                   os.path.join(db_path,'W_C2T')
                                                                                                                                                   ,tmp_path,outfile_2FCT,tmp_path,outfile_1RCT,tmp_path,WC2T_rf),shell=True)

            bowtie_map4=Popen('%sbowtie -e %d --nomaqround --fr -k 2 %s %s --quiet --best --suppress 2,5,6 -p 2 %s -f -1 %s%s -2 %s%s %s%s '%(bowtie_path,40*int_no_mismatches,min_string, max_string,
                                                                                                                                                   os.path.join(db_path,'C_C2T'),
                                                                                                                                                   tmp_path,outfile_2FCT,tmp_path,outfile_1RCT,tmp_path,CC2T_rf),shell=True)

            bowtie_map1.wait()
            bowtie_map2.wait()
            bowtie_map3.wait()
            bowtie_map4.wait()

            Popen('rm %s %s %s %s &'%(tmp_path+outfile_1FCT,tmp_path+outfile_2FCT,tmp_path+outfile_1RCT,tmp_path+outfile_2RCT),shell=True)
            #--------------------------------------------------------------------------------
            # Post processing
            #--------------------------------------------------------------------------------



            ali_file1=WC2T_fr # mapped to Watson strand, has methylation info of Watson strand
            ali_file2=WC2T_rf # mapped to Watson strand, has methylation info of Watson strand
            ali_file3=CC2T_fr # mapped to Crick strand, has methylation info of Crick strand
            ali_file4=CC2T_rf # mapped to Crick strand, has methylation info of Crick strand

            #----------------------------------------------------------------

            FW_C2T_fr_U,FW_C2T_fr_R=extract_mapping(tmp_path+ali_file1)
            FW_C2T_rf_U,FW_C2T_rf_R=extract_mapping(tmp_path+ali_file2)
            RC_C2T_fr_U,RC_C2T_fr_R=extract_mapping(tmp_path+ali_file3)
            RC_C2T_rf_U,RC_C2T_rf_R=extract_mapping(tmp_path+ali_file4)

            Popen('rm %s %s %s %s &'%(tmp_path+WC2T_fr, tmp_path+WC2T_rf,tmp_path+CC2T_fr,tmp_path+CC2T_rf),shell=True)

            #----------------------------------------------------------------
            # get uniq-hit reads
            #----------------------------------------------------------------
            Union_set=set(FW_C2T_fr_U.keys()) | set(FW_C2T_rf_U.keys()) | set(RC_C2T_fr_U.keys()) | set(RC_C2T_rf_U.keys())

            Unique_FW_fr_C2T=set() # +
            Unique_FW_rf_C2T=set() # +
            Unique_RC_fr_C2T=set() # -
            Unique_RC_rf_C2T=set() # -


            for x in Union_set:
                list=[]
                for d in [FW_C2T_fr_U,FW_C2T_rf_U,RC_C2T_fr_U,RC_C2T_rf_U]:
                    mis_lst=d.get(x,[99])
                    mis=int(mis_lst[0])
                    list.append(mis)
                for d in [FW_C2T_fr_R,FW_C2T_rf_R,RC_C2T_fr_R,RC_C2T_rf_R]:
                    mis=d.get(x,99)
                    list.append(mis)
                mini=min(list)
                if list.count(mini)==1:
                    mini_index=list.index(mini)
                    if mini_index==0:
                        Unique_FW_fr_C2T.add(x)
                    elif mini_index==1:
                        Unique_FW_rf_C2T.add(x)
                    elif mini_index==2:
                        Unique_RC_fr_C2T.add(x)
                    elif mini_index==3:
                        Unique_RC_rf_C2T.add(x)

            Union_set=set()

            FW_C2T_fr_R={}
            FW_C2T_rf_R={}
            RC_C2T_fr_R={}
            RC_C2T_rf_R={}

            FW_C2T_fr_uniq_lst=[[FW_C2T_fr_U[u][1],u] for u in Unique_FW_fr_C2T]
            FW_C2T_rf_uniq_lst=[[FW_C2T_rf_U[u][1],u] for u in Unique_FW_rf_C2T]
            RC_C2T_fr_uniq_lst=[[RC_C2T_fr_U[u][1],u] for u in Unique_RC_fr_C2T]
            RC_C2T_rf_uniq_lst=[[RC_C2T_rf_U[u][1],u] for u in Unique_RC_rf_C2T]

            FW_C2T_fr_uniq_lst.sort()
            FW_C2T_rf_uniq_lst.sort()
            RC_C2T_fr_uniq_lst.sort()
            RC_C2T_rf_uniq_lst.sort()
            FW_C2T_fr_uniq_lst=[x[1] for x in FW_C2T_fr_uniq_lst]
            FW_C2T_rf_uniq_lst=[x[1] for x in FW_C2T_rf_uniq_lst]
            RC_C2T_fr_uniq_lst=[x[1] for x in RC_C2T_fr_uniq_lst]
            RC_C2T_rf_uniq_lst=[x[1] for x in RC_C2T_rf_uniq_lst]
            #----------------------------------------------------------------

            n1=len(Unique_FW_fr_C2T)
            n2=len(Unique_FW_rf_C2T)
            n3=len(Unique_RC_fr_C2T)
            n4=len(Unique_RC_rf_C2T)
            n12=n1+n2+n3+n4

            numbers_premapped_lst[0]+=n1
            numbers_premapped_lst[1]+=n2
            numbers_premapped_lst[2]+=n3
            numbers_premapped_lst[3]+=n4

            Unique_FW_fr_C2T=set()
            Unique_FW_rf_C2T=set()
            Unique_RC_fr_C2T=set()
            Unique_RC_rf_C2T=set()

            #logoutf.write("U -- %d FW-RC strand bs-unique pairs (mapped to Watson)"%(n1)+"\n")
            #logoutf.write("U -- %d RC-FW strand bs-unique pairs (mapped to Crick)"%(n2)+"\n")
            #logoutf.write("U -- %d bs-unique pairs"%(n12)+"\n")
            #logoutf.write("-------------------------------- "+'\n')

            #print "# %10d FW-RC bs-unique reads (mapped to Watson)"%(n1);
            #print "# %10d RC-FW bs-unique reads (mapped to Watson)"%(n2);
            #print "# %10d FW-RC bs-unique reads (mapped to Crick)"%(n3);
            #print "# %10d RC-FW bs-unique reads (mapped to Crick)"%(n4);
            #------ Original BS reads --------------------------------------
            original_bs_reads_1={}
            original_bs_reads_2={}

            original_bs_reads_lst=[original_bs_reads_1,original_bs_reads_2]
            outfile_BS_list=[outfile_1BS,outfile_2BS]

            for i in range(2):
                original_bs_reads=original_bs_reads_lst[i]
                outfile_BS=outfile_BS_list[i]
                for line in fileinput.input(tmp_path+outfile_BS):
                    l=line.split()
                    if len(l)==2:
                        original_bs_reads[str(l[0])]=str(l[1])
                fileinput.close()

            #print "# raw BS reads end 1: %d ,end 2: %d"%(len(original_bs_reads_1),len(original_bs_reads_2));
            #logoutf.write("# raw BS reads end 1: %d ,end 2: %d"%(len(original_bs_reads_1),len(original_bs_reads_2))+"\n")

            Popen('rm %s %s &'%(tmp_path+outfile_1BS,tmp_path+outfile_2BS),shell=True)

            #----------------------------------------------------------------

            nn=0
            for ali_unique_lst, ali_dic in [(FW_C2T_fr_uniq_lst,FW_C2T_fr_U),(FW_C2T_rf_uniq_lst,FW_C2T_rf_U),(RC_C2T_fr_uniq_lst,RC_C2T_fr_U),(RC_C2T_rf_uniq_lst,RC_C2T_rf_U)]:
                nn+=1

                mapped_chr0=""
                for header in ali_unique_lst:

                    l=ali_dic[header]

                    mapped_chr=str(l[1])
                    mapped_location_1=int(l[2])
                    mapped_location_2=int(l[3])
                    #-------------------------------------
                    if mapped_chr != mapped_chr0:
                        my_gseq=genome_seqs[mapped_chr]
                        chr_length=len(my_gseq)
                        mapped_chr0=mapped_chr
                    #-------------------------------------
                    all_mapped+=1

                    if nn==1: 							# FW-RC mapped to + strand:
                        original_BS_1=original_bs_reads_1[header]
                        original_BS_2=reverse_compl_seq(original_bs_reads_2[header])
                        FR="+FR"
                        mapped_location_1 += 1
                        origin_genome_long_1=my_gseq[mapped_location_1-2-1:mapped_location_1+len(original_BS_1)+2-1]
                        origin_genome_long_1=origin_genome_long_1.upper()
                        origin_genome_1=origin_genome_long_1[2:-2]
                        mapped_strand_1="+"

                        mapped_location_2 += 1
                        origin_genome_long_2=my_gseq[mapped_location_2-2-1:mapped_location_2+len(original_BS_2)+2-1]
                        origin_genome_long_2=origin_genome_long_2.upper()
                        origin_genome_2=origin_genome_long_2[2:-2]
                        mapped_strand_2="+"

                    elif nn==2: 							# RC-FW mapped to + strand:
                        original_BS_1=original_bs_reads_2[header]
                        original_BS_2=reverse_compl_seq(original_bs_reads_1[header])
                        FR="+RF"
                        mapped_location_1 += 1
                        origin_genome_long_1=my_gseq[mapped_location_1-2-1:mapped_location_1+len(original_BS_1)+2-1]
                        origin_genome_long_1=origin_genome_long_1.upper()
                        origin_genome_1=origin_genome_long_1[2:-2]
                        mapped_strand_1="+"

                        mapped_location_2 += 1
                        origin_genome_long_2=my_gseq[mapped_location_2-2-1:mapped_location_2+len(original_BS_2)+2-1]
                        origin_genome_long_2=origin_genome_long_2.upper()
                        origin_genome_2=origin_genome_long_2[2:-2]
                        mapped_strand_2="+"


                    elif nn==3: 						# FW-RC mapped to - strand:
                        original_BS_1=original_bs_reads_1[header]
                        original_BS_2=reverse_compl_seq(original_bs_reads_2[header])

                        FR="-FR"
                        mapped_location_1=chr_length-mapped_location_1-len(original_BS_1)+1
                        origin_genome_long_1=my_gseq[mapped_location_1-2-1:mapped_location_1+len(original_BS_1)+2-1]
                        origin_genome_long_1=reverse_compl_seq(origin_genome_long_1)
                        origin_genome_1=origin_genome_long_1[2:-2]
                        mapped_strand_1="-"

                        mapped_location_2=chr_length-mapped_location_2-len(original_BS_2)+1
                        origin_genome_long_2=reverse_compl_seq(my_gseq[mapped_location_2-2-1:mapped_location_2+len(original_BS_2)+2-1])
                        origin_genome_long_2=origin_genome_long_2.upper()
                        origin_genome_2=origin_genome_long_2[2:-2]
                        mapped_strand_2="-"

                    elif nn==4: 						# RC-FW mapped to - strand:
                        original_BS_1=original_bs_reads_2[header]
                        original_BS_2=reverse_compl_seq(original_bs_reads_1[header])

                        FR="-RF"
                        mapped_location_1=chr_length-mapped_location_1-len(original_BS_1)+1
                        origin_genome_long_1=my_gseq[mapped_location_1-2-1:mapped_location_1+len(original_BS_1)+2-1]
                        origin_genome_long_1=reverse_compl_seq(origin_genome_long_1)
                        origin_genome_1=origin_genome_long_1[2:-2]
                        mapped_strand_1="-"

                        mapped_location_2=chr_length-mapped_location_2-len(original_BS_2)+1
                        origin_genome_long_2=reverse_compl_seq(my_gseq[mapped_location_2-2-1:mapped_location_2+len(original_BS_2)+2-1])
                        origin_genome_long_2=origin_genome_long_2.upper()
                        origin_genome_2=origin_genome_long_2[2:-2]
                        mapped_strand_2="-"

                    N_mismatch_1=N_MIS(original_BS_1,origin_genome_1)
                    N_mismatch_2=N_MIS(original_BS_2,origin_genome_2)
                    if max(N_mismatch_1,N_mismatch_2) <= int(indexname) :
                        all_mapped_passed+=1
                        numbers_mapped_lst[nn-1]+=1
                        #---- unmapped -------------------------
                        del original_bs_reads_1[header]
                        del original_bs_reads_2[header]
                        #---------------------------------------
                        mapped_location_1=str(mapped_location_1).zfill(10)
                        mapped_location_2=str(mapped_location_2).zfill(10)

                        coordinate_1=mapped_chr+mapped_strand_1+mapped_location_1
                        coordinate_2=mapped_chr+mapped_strand_2+mapped_location_2
                        output_genome_1=origin_genome_long_1[0:2]+"_"+origin_genome_1+"_"+origin_genome_long_1[-2:]
                        output_genome_2=origin_genome_long_2[0:2]+"_"+origin_genome_2+"_"+origin_genome_long_2[-2:]

                        methy_1=methy_seq(original_BS_1,output_genome_1)
                        methy_2=methy_seq(original_BS_2,output_genome_2)

                        mC_lst,uC_lst=mcounts(methy_1,mC_lst,uC_lst)
                        mC_lst,uC_lst=mcounts(methy_2,mC_lst,uC_lst)

                        #---STEVE FILTER----------------
                        condense_seq_1=methy_1.replace('-','')
                        STEVE_1=0
                        if "ZZZ" in condense_seq_1:
                            STEVE_1=1

                        condense_seq_2=methy_2.replace('-','')
                        STEVE_2=0
                        if "ZZZ" in condense_seq_2:
                            STEVE_2=1

                        outf.write('%s/1	%2d 	%3s	%s	%s	%s	%s	%d'%(header,N_mismatch_1,FR,coordinate_1,output_genome_1,original_BS_1,methy_1,STEVE_1)+"\n")
                        outf.write('%s/2	%2d 	%3s	%s	%s	%s	%s	%d'%(header,N_mismatch_2,FR,coordinate_2,output_genome_2,original_BS_2,methy_2,STEVE_2)+"\n")

            print "--> %s %s (%d/%d) "%(read_file_1,read_file_2,no_my_files,len(my_files));
            #----------------------------------------------------------------
            #	output unmapped pairs
            #----------------------------------------------------------------

            unmapped_lst=original_bs_reads_1.keys()
            unmapped_lst.sort()

            for u in unmapped_lst:
                outf_u1.write("%s"%(original_bs_reads_1[u])+"\n")
                outf_u2.write("%s"%(original_bs_reads_2[u])+"\n")

            all_unmapped+=len(unmapped_lst)





        if asktag=="N":

            #----------------------------------------------------------------
            outfile_1BS = 'Trimed_BS_1.fa'+random_id
            outfile_1FCT= 'Trimed_FCT_1.fa'+random_id
            outfile_2BS = 'Trimed_BS_2.fa'+random_id
            outfile_2FCT= 'Trimed_FCT_2.fa'+random_id

            read_inf=open(tmp_path+read_file_1,"r")
            oneline=read_inf.readline()
            l=oneline.split()
            input_format=""

            #if len(l)==5: # old solexa format
            #	input_format="old Solexa Seq file"

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

            print "Detected data format: %s" % input_format

            #----------------------------------------------------------------
            read_file_list=[read_file_1,read_file_2]
            outfile_BS_list=[outfile_1BS,outfile_2BS]
            outfile_FCT_list=[outfile_1FCT,outfile_2FCT]
            n_list=[0,0]


            for f in range(2):
                read_file=read_file_list[f]
                outf_BS=open(tmp_path+outfile_BS_list[f],'w')
                outf_FCT=open(tmp_path+outfile_FCT_list[f],'w')
                n=n_list[f]

                id=""
                seq=""
                seq_ready="N"
                for line in fileinput.input(tmp_path+read_file):
                    l=line.split()
                    if input_format=="old Solexa Seq file":
                        n+=1
                        id=str(n)
                        id=id.zfill(12)
                        seq=l[4]
                        seq_ready="Y"
                    elif input_format=="list of sequences":
                        n+=1
                        seq_ready="N"
                        id=str(n)
                        id=id.zfill(12)
                        seq=l[0]
                        seq_ready="Y"
                    elif input_format=="FastQ":
                        m_fastq=math.fmod(n_fastq,4)
                        n_fastq+=1
                        seq_ready="N"
                        if m_fastq==0:
                            n+=1
                            id=str(n)
                            id=id.zfill(12)
                            seq=""
                        elif m_fastq==1:
                            seq=l[0]
                            seq_ready="Y"
                        else:
                            seq=""
                    elif input_format=="Illumina GAII qseq file":
                        n+=1
                        seq_ready="N"
                        id=str(n)
                        id=id.zfill(12)
                        seq=l[8]
                        seq_ready="Y"
                    elif input_format=="fasta":
                        m_fasta=math.fmod(n_fasta,2)
                        n_fasta+=1
                        seq_ready="N"
                        if m_fasta==0:
                            n+=1
                            #id=str(n)
                            #id=id.zfill(12)
                            id=l[0][1:]
                            id=id.zfill(17)
                            seq=""
                        elif m_fasta==1:
                            seq=l[0]
                            seq_ready="Y"
                        else:
                            seq=""
                    #----------------------------------------------------------------
                    if seq_ready=="Y":
                        seq=seq[cut1-1:cut2] #<----------------------selecting 0..52 from 1..72  -e 52
                        seq=seq.upper()
                        seq=seq.replace(".","N")

                        #--striping BS adapter from 3' read --------------------------------------------------------------
                        if (adapterA !="") and (adapterB !=""):
                            signature=adapterA[:6]
                            if signature in seq:
                                signature_pos=seq.index(signature)
                                if seq[signature_pos:] in adapterA:
                                    seq=seq[:signature_pos]#+"".join(["N" for x in range(len(seq)-len(signature_pos))])
                                    all_trimed+=1
                            else:
                                signature=adapterB[:6]
                                if signature in seq:
                                    #print id,seq,signature;
                                    signature_pos=seq.index(signature)
                                    if seq[signature_pos:] in adapterB:
                                        seq=seq[:signature_pos]#+"".join(["N" for x in range(len(seq)-len(signature_pos))])
                                        all_trimed+=1

                        if len(seq)<=4:
                            seq=''.join(["N" for x in range(cut2-cut1+1)])
                        #---------  trimmed_raw_BS_read  ------------------
                        outf_BS.write('%s\t%s'%(id,seq)+"\n")

                        #---------  FW_C2T  ------------------
                        FWseq=copy.deepcopy(seq)
                        if f==0:
                            FWseq1=FWseq.replace("C","T")
                            outf_FCT.write('>%s'% id +"\n")
                            outf_FCT.write('%s' % FWseq1 +"\n")
                        elif f==1:
                            RCseq=reverse_compl_seq(FWseq)
                            RCseq=RCseq.replace("C","T")
                            outf_FCT.write('>%s'% id +"\n")
                            outf_FCT.write('%s'% RCseq +"\n")

                n_list[f]=n
                outf_BS.close()
                outf_FCT.close()

                fileinput.close()


            #print "All input end 1: %d , end 2: %d "%(n_list[0],n_list[1]);
            all_raw_reads+=n
            elapsed1 = (time.time() - start)

            #--------------------------------------------------------------------------------
            # Bowtie mapping
            #--------------------------------------------------------------------------------
            WC2T_fr="W_C2T_fr_m"+str(indexname)+".mapping"+random_id
            CC2T_fr="C_C2T_fr_m"+str(indexname)+".mapping"+random_id

            bowtie_map1=Popen('%sbowtie -e %d --nomaqround --norc --ff -k 2 %s %s --quiet --best --suppress 2,5,6 -p 4 %s -f -1 %s%s -2 %s%s %s%s '%(bowtie_path,40*int_no_mismatches,min_string, max_string,
                                                                                                                                                          os.path.join(db_path,'W_C2T'),
                                                                                                                                                          tmp_path,outfile_1FCT,tmp_path,outfile_2FCT,tmp_path,WC2T_fr),shell=True)

            bowtie_map2=Popen('%sbowtie -e %d --nomaqround --norc --ff -k 2 %s %s --quiet --best --suppress 2,5,6 -p 4 %s -f -1 %s%s -2 %s%s %s%s '%(bowtie_path,40*int_no_mismatches,min_string, max_string,
                                                                                                                                                          os.path.join(db_path,'C_C2T'),
                                                                                                                                                          tmp_path,outfile_1FCT,tmp_path,outfile_2FCT,tmp_path,CC2T_fr),shell=True)


            bowtie_map1.wait()
            bowtie_map2.wait()

            Popen('rm %s %s &'%(tmp_path+outfile_1FCT,tmp_path+outfile_2FCT),shell=True)
            #--------------------------------------------------------------------------------
            # Post processing
            #--------------------------------------------------------------------------------



            ali_file1=WC2T_fr # mapped to Watson strand, has methylation info of Watson strand
            ali_file3=CC2T_fr # mapped to Crick strand, has methylation info of Crick strand

            #----------------------------------------------------------------

            FW_C2T_fr_U,FW_C2T_fr_R=extract_mapping(tmp_path+ali_file1)
            RC_C2T_fr_U,RC_C2T_fr_R=extract_mapping(tmp_path+ali_file3)

            #Popen('rm %s %s &'%(tmp_path+WC2T_fr, tmp_path+CC2T_fr),shell=True)

            #----------------------------------------------------------------
            # get uniq-hit reads
            #----------------------------------------------------------------
            Union_set=set(FW_C2T_fr_U.keys()) | set(RC_C2T_fr_U.keys())

            Unique_FW_fr_C2T=set() # +
            Unique_RC_fr_C2T=set() # -


            for x in Union_set:
                list=[]
                for d in [FW_C2T_fr_U,RC_C2T_fr_U]:
                    mis_lst=d.get(x,[99])
                    mis=int(mis_lst[0])
                    list.append(mis)
                for d in [FW_C2T_fr_R,RC_C2T_fr_R]:
                    mis=d.get(x,99)
                    list.append(mis)
                mini=min(list)
                if list.count(mini)==1:
                    mini_index=list.index(mini)
                    if mini_index==0:
                        Unique_FW_fr_C2T.add(x)
                    elif mini_index==1:
                        Unique_RC_fr_C2T.add(x)

            Union_set=set()

            FW_C2T_fr_R={}
            RC_C2T_fr_R={}

            FW_C2T_fr_uniq_lst=[[FW_C2T_fr_U[u][1],u] for u in Unique_FW_fr_C2T]
            RC_C2T_fr_uniq_lst=[[RC_C2T_fr_U[u][1],u] for u in Unique_RC_fr_C2T]

            FW_C2T_fr_uniq_lst.sort()
            RC_C2T_fr_uniq_lst.sort()
            FW_C2T_fr_uniq_lst=[x[1] for x in FW_C2T_fr_uniq_lst]
            RC_C2T_fr_uniq_lst=[x[1] for x in RC_C2T_fr_uniq_lst]
            #----------------------------------------------------------------

            n1=len(Unique_FW_fr_C2T)
            n3=len(Unique_RC_fr_C2T)
            n12=n1+n3

            numbers_premapped_lst[0]+=n1
            numbers_premapped_lst[1]+=n3

            Unique_FW_fr_C2T=set()
            Unique_RC_fr_C2T=set()

            #logoutf.write("U -- %d FW-RC strand bs-unique pairs (mapped to Watson)"%(n1)+"\n")
            #logoutf.write("U -- %d RC-FW strand bs-unique pairs (mapped to Crick)"%(n2)+"\n")
            #logoutf.write("U -- %d bs-unique pairs"%(n12)+"\n")
            #logoutf.write("-------------------------------- "+'\n')

            #print "# %10d FW-RC bs-unique reads (mapped to Watson)"%(n1);
            #print "# %10d RC-FW bs-unique reads (mapped to Watson)"%(n2);
            #print "# %10d FW-RC bs-unique reads (mapped to Crick)"%(n3);
            #print "# %10d RC-FW bs-unique reads (mapped to Crick)"%(n4);
            #------ Original BS reads --------------------------------------
            original_bs_reads_1={}
            original_bs_reads_2={}

            original_bs_reads_lst=[original_bs_reads_1,original_bs_reads_2]
            outfile_BS_list=[outfile_1BS,outfile_2BS]

            for i in range(2):
                original_bs_reads=original_bs_reads_lst[i]
                outfile_BS=outfile_BS_list[i]
                for line in fileinput.input(tmp_path+outfile_BS):
                    l=line.split()
                    if len(l)==2:
                        original_bs_reads[str(l[0])]=str(l[1])
                fileinput.close()

            #print "# raw BS reads end 1: %d ,end 2: %d"%(len(original_bs_reads_1),len(original_bs_reads_2));
            #logoutf.write("# raw BS reads end 1: %d ,end 2: %d"%(len(original_bs_reads_1),len(original_bs_reads_2))+"\n")

            Popen('rm %s %s &'%(tmp_path+outfile_1BS,tmp_path+outfile_2BS),shell=True)

            #----------------------------------------------------------------

            nn=0
            for ali_unique_lst, ali_dic in [(FW_C2T_fr_uniq_lst,FW_C2T_fr_U),(RC_C2T_fr_uniq_lst,RC_C2T_fr_U)]:
                nn+=1
                mapped_chr0=""
                for header in ali_unique_lst:
                    l=ali_dic[header]
                    mapped_chr=str(l[1])
                    mapped_location_1=int(l[2])
                    mapped_location_2=int(l[3])
                    #-------------------------------------
                    if mapped_chr != mapped_chr0:
                        my_gseq=genome_seqs[mapped_chr]
                        chr_length=len(my_gseq)
                        mapped_chr0=mapped_chr
                    #-------------------------------------
                    all_mapped+=1

                    if nn==1: 							# FW-RC mapped to + strand:
                        original_BS_1=original_bs_reads_1[header]
                        original_BS_2=reverse_compl_seq(original_bs_reads_2[header])
                        FR="+FR"
                        mapped_location_1 += 1
                        origin_genome_long_1=my_gseq[mapped_location_1-2-1:mapped_location_1+len(original_BS_1)+2-1]
                        origin_genome_long_1=origin_genome_long_1.upper()
                        origin_genome_1=origin_genome_long_1[2:-2]
                        mapped_strand_1="+"

                        mapped_location_2 += 1
                        origin_genome_long_2=my_gseq[mapped_location_2-2-1:mapped_location_2+len(original_BS_2)+2-1]
                        origin_genome_long_2=origin_genome_long_2.upper()
                        origin_genome_2=origin_genome_long_2[2:-2]
                        mapped_strand_2="+"

                    elif nn==2: 						# FW-RC mapped to - strand:
                        original_BS_1=original_bs_reads_1[header]
                        original_BS_2=reverse_compl_seq(original_bs_reads_2[header])

                        FR="-FR"
                        mapped_location_1=chr_length-mapped_location_1-len(original_BS_1)+1
                        origin_genome_long_1=my_gseq[mapped_location_1-2-1:mapped_location_1+len(original_BS_1)+2-1]
                        origin_genome_long_1=reverse_compl_seq(origin_genome_long_1)
                        origin_genome_1=origin_genome_long_1[2:-2]
                        mapped_strand_1="-"

                        mapped_location_2=chr_length-mapped_location_2-len(original_BS_2)+1
                        origin_genome_long_2=reverse_compl_seq(my_gseq[mapped_location_2-2-1:mapped_location_2+len(original_BS_2)+2-1])
                        origin_genome_long_2=origin_genome_long_2.upper()
                        origin_genome_2=origin_genome_long_2[2:-2]
                        mapped_strand_2="-"


                    N_mismatch_1=N_MIS(original_BS_1,origin_genome_1)
                    N_mismatch_2=N_MIS(original_BS_2,origin_genome_2)
                    if max(N_mismatch_1,N_mismatch_2) <= int(indexname):
                        numbers_mapped_lst[nn-1]+=1
                        all_mapped_passed+=1
                        #---- unmapped -------------------------
                        del original_bs_reads_1[header]
                        del original_bs_reads_2[header]
                        #---------------------------------------
                        mapped_location_1=str(mapped_location_1).zfill(10)
                        mapped_location_2=str(mapped_location_2).zfill(10)

                        coordinate_1=mapped_chr+mapped_strand_1+mapped_location_1
                        coordinate_2=mapped_chr+mapped_strand_2+mapped_location_2
                        output_genome_1=origin_genome_long_1[0:2]+"_"+origin_genome_1+"_"+origin_genome_long_1[-2:]
                        output_genome_2=origin_genome_long_2[0:2]+"_"+origin_genome_2+"_"+origin_genome_long_2[-2:]

                        methy_1=methy_seq(original_BS_1,output_genome_1)
                        methy_2=methy_seq(original_BS_2,output_genome_2)

                        mC_lst,uC_lst=mcounts(methy_1,mC_lst,uC_lst)
                        mC_lst,uC_lst=mcounts(methy_2,mC_lst,uC_lst)

                        #---STEVE FILTER----------------
                        condense_seq_1=methy_1.replace('-','')
                        STEVE_1=0
                        if "ZZZ" in condense_seq_1:
                            STEVE_1=1

                        condense_seq_2=methy_2.replace('-','')
                        STEVE_2=0
                        if "ZZZ" in condense_seq_2:
                            STEVE_2=1

                        outf.write('%s/1	%2d 	%3s	%s	%s	%s	%s	%d'%(header,N_mismatch_1,FR,coordinate_1,output_genome_1,original_BS_1,methy_1,STEVE_1)+"\n")
                        outf.write('%s/2	%2d 	%3s	%s	%s	%s	%s	%d'%(header,N_mismatch_2,FR,coordinate_2,output_genome_2,original_BS_2,methy_2,STEVE_2)+"\n")

            print "--> %s %s (%d/%d) "%(read_file_1,read_file_2,no_my_files,len(my_files));
            #----------------------------------------------------------------
            #	output unmapped pairs
            #----------------------------------------------------------------

            unmapped_lst=original_bs_reads_1.keys()
            unmapped_lst.sort()

            for u in unmapped_lst:
                outf_u1.write("%s"%(original_bs_reads_1[u])+"\n")
                outf_u2.write("%s"%(original_bs_reads_2[u])+"\n")

            all_unmapped+=len(unmapped_lst)


    #==================================================================================================
    outf.close()

    outf_u1.close()
    outf_u2.close()
    Popen('rm -r %s '% tmp_path,shell=True)

    logoutf.write("-------------------------------- "+'\n')
    logoutf.write("O Number of raw BS-read pairs: %d ( %d bp)"%(all_raw_reads,cut2-cut1+1)+"\n")
    logoutf.write("O Number of ends trimmed for adapter: %d"% all_trimed+"\n")

    if all_raw_reads >0:

        logoutf.write("O Number of unique-hits read pairs for post-filtering: %d" % all_mapped + "\n")
        if asktag=="Y":
            logoutf.write("O -- %7d FW-RC pairs mapped to Watson strand (before post-filtering)"%(numbers_premapped_lst[0])+"\n")
            logoutf.write("O -- %7d RC-FW pairs mapped to Watson strand (before post-filtering)"%(numbers_premapped_lst[1])+"\n")
            logoutf.write("O -- %7d FW-RC pairs mapped to Crick strand (before post-filtering)"%(numbers_premapped_lst[2])+"\n")
            logoutf.write("O -- %7d RC-FW pairs mapped to Crick strand (before post-filtering)"%(numbers_premapped_lst[3])+"\n")
        elif asktag=="N":
            logoutf.write("O -- %7d FW-RC pairs mapped to Watson strand (before post-filtering)"%(numbers_premapped_lst[0])+"\n")
            logoutf.write("O -- %7d FW-RC pairs mapped to Crick strand (before post-filtering)"%(numbers_premapped_lst[1])+"\n")


        logoutf.write("O --- %d uniqlely aligned pairs, where each end has mismatches <= %s"%(all_mapped_passed,str(indexname))+"\n")
        if asktag=="Y":
            logoutf.write("O ----- %7d FW-RC pairs mapped to Watson strand"%(numbers_mapped_lst[0])+"\n")
            logoutf.write("O ----- %7d RC-FW pairs mapped to Watson strand"%(numbers_mapped_lst[1])+"\n")
            logoutf.write("O ----- %7d FW-RC pairs mapped to Crick strand"%(numbers_mapped_lst[2])+"\n")
            logoutf.write("O ----- %7d RC-FW pairs mapped to Crick strand"%(numbers_mapped_lst[3])+"\n")
        elif asktag=="N":
            logoutf.write("O ----- %7d FW-RC pairs mapped to Watson strand"%(numbers_mapped_lst[0])+"\n")
            logoutf.write("O ----- %7d FW-RC pairs mapped to Crick strand"%(numbers_mapped_lst[1])+"\n")

        logoutf.write("O Mapability= %1.4f%%"%(100*float(all_mapped_passed)/all_raw_reads)+"\n")
        logoutf.write("O Unmapped read pairs: %d"% all_unmapped+"\n")


        n_CG=mC_lst[0]+uC_lst[0]
        n_CHG=mC_lst[1]+uC_lst[1]
        n_CHH=mC_lst[2]+uC_lst[2]

        logoutf.write("-------------------------------- "+'\n')
        logoutf.write("Methylated C in mapped reads "+'\n')
        logoutf.write(" mCG %1.3f%%"%((100*float(mC_lst[0])/n_CG) if n_CG != 0 else 0)+'\n')
        logoutf.write(" mCHG %1.3f%%"%((100*float(mC_lst[1])/n_CHG) if n_CHG != 0 else 0)+'\n')
        logoutf.write(" mCHH %1.3f%%"%((100*float(mC_lst[2])/n_CHH) if n_CHH != 0 else 0)+'\n')

    logoutf.write("----------------------------------------------"+"\n")
    logoutf.write("------------------- END --------------------"+"\n")
    logoutf.close()
    elapsed("=== END %s %s ===" % (main_read_file_1, main_read_file_2))
