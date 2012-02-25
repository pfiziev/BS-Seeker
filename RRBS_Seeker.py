import fileinput, string,os, gzip,copy, time, shelve,subprocess, random, math, os.path
from subprocess import Popen

#----------------------------------------------------------------
def reverse_compl_seq(strseq):
	strseq=strseq.upper()
	rc_strseq=strseq.translate(string.maketrans("ATCG", "TAGC"))[::-1]
	return rc_strseq;
	
def compl_seq(strseq):
	strseq=strseq.upper()
	c_strseq=strseq.translate(string.maketrans("ATCG", "TAGC"))
	return c_strseq;
	
def N_MIS(r,g):
	combo=[]
	if len(r)==len(g):
		combo=[r[i]+g[i] for i in range(len(r)) if r[i]!=g[i] and r[i] !="N" and g[i]!="N"]
		combo=[x for x in combo if x !="TC"]
	return len(combo);
	

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
			print l;
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
	return U,R;
	


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
	return out_serial,out_start,out_end;


	
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
					
		m_seq=m_seq+xx
	return m_seq;

def mcounts(mseq,mlst,ulst):
	out_mlst=[mlst[0]+mseq.count("X"),mlst[1]+mseq.count("Y"),mlst[2]+mseq.count("Z")]
	out_ulst=[ulst[0]+mseq.count("x"),ulst[1]+mseq.count("y"),ulst[2]+mseq.count("z")]
	return out_mlst,out_ulst;
	
#----------------------------------------------------------------

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i", "--input", type="string", dest="infilename",help="Input your read file name (FORMAT: sequences, illumina fastq, qseq,fasta)", metavar="FILE")

parser.set_defaults(taginfo="CGG/TGG")
parser.add_option("-t", "--tag", type="string",dest="taginfo",help="Msp-I tag: CGG TGG CGA or CGG/TGG (both)", metavar="TAG")

parser.set_defaults(cutnumber2=100)
parser.add_option("-e","--e",type = "int",dest = "cutnumber2", help="The last cycle number of your read to be mapped [36]") 

parser.set_defaults(adapterfilename="")
parser.add_option("-a", "--adapter", type="string",dest="adapterfilename",help="Input your adapter file", metavar="FILE")

parser.set_defaults(bowtiepath="~/bowtie-0.12.7/")
parser.add_option("-p", "--path", dest="bowtiepath",help="Path to Bowtie [~/bowtie-0.12.7/]", metavar="PATH")

parser.set_defaults(dbpath="/u/home/mcdb/paoyang/genome/mouse/RRBS_reference_genome_100_300/")
parser.add_option("-d", "--db", dest="dbpath",help="Path to Reference genome library (generated in preprocessing genome) [./RRBS_reference_genome_100_300/]", metavar="DBPATH")


parser.set_defaults(indexname=3)
parser.add_option("-m", "--mis", dest="indexname",help="Number of mismatches (0, 1, 2, 3) [3]")

parser.set_defaults(no_split=4000000)
parser.add_option("-l", "--split_line",type = "int", dest="no_split",help="Number of lines per split (the read file will be split into small files for mapping. The result will be merged. [4000000]")

parser.set_defaults(outfilename="RRBS_SEEKER_OUTPUT.txt")
parser.add_option("-o", "--output", dest="outfilename",help="The name of output file [RRBS_SEEKER_OUTPUT.txt]", metavar="OUTFILE")


(options, args) = parser.parse_args()
                  
main_read_file=options.infilename

mytag=options.taginfo
mytag=mytag.upper()
mytag_lst=mytag.split("/")

adapter_file=options.adapterfilename

out_filename=options.outfilename
out_filename=str(out_filename)

cut2=options.cutnumber2

no_small_lines=options.no_split

indexname=options.indexname
int_no_mismatches=min(int(indexname),cut2)
indexname=str(int_no_mismatches)


bowtie_path=options.bowtiepath
bowtie_path=str(bowtie_path)
if bowtie_path[-1] !="/":
	bowtie_path=bowtie_path+"/"

db_path=options.dbpath
db_path=str(db_path)
if db_path[-1] !="/":
	db_path=db_path+"/"

regions_file="RRBS_mapable_regions.shelve"

reads_path="./"

#----------------------------------------------------------------
os.system('mkdir ./%s-TMP/'%(main_read_file))
tmp_path='./%s-TMP/'%(main_read_file)

#----------------------------------------------------------------
adapter_seq=""
if adapter_file !="":
	adapter_inf=open(adapter_file,"r")
	adapter_seq=adapter_inf.readline()
	adapter_inf.close()
	adapter_seq=adapter_seq.rstrip("\n")

#----------------------------------------------------------------
# splitting the big read file


splitting=Popen('split -l %d %s %s%s-s-'%(no_small_lines,main_read_file,tmp_path,main_read_file),shell=True)
splitting.wait()

dirList=os.listdir(tmp_path)
my_files=[]
for splitted_file in dirList:
	if splitted_file.startswith("%s-s-"%(main_read_file)):
		my_files.append(splitted_file)

my_files.sort()


#----------------------------------------------------------------
# output files

outfile=out_filename
outf=open(outfile,'w')

logoutfile='log_RRBS_Seeker_SE_'+out_filename
logoutf=open(logoutfile,'w')



#----------------------------------------------------------------
print "Read filename: %s"%(main_read_file)
print "Starting Msp-1 tag: %s"%(mytag)
print "The last cycle (for mapping): %d"%(cut2)
print "Bowtie path: %s"%(bowtie_path)
print "Reference genome library path: %s"%(db_path)
print "Number of mismatches allowed: %s"%(indexname)

logoutf.write("I Read filename: %s"%(main_read_file)+"\n")
logoutf.write("I Starting Msp-1 tag: %s"%(mytag)+"\n")
logoutf.write("I The last cycle (for mapping): %d"%(cut2)+"\n")
logoutf.write("I Bowtie path: %s"%(bowtie_path)+"\n")
logoutf.write("I Reference genome library path: %s"%(db_path)+"\n")
logoutf.write("I Number of mismatches allowed: %s"%(indexname)+"\n")
logoutf.write("I adapter seq: %s"%(adapter_seq)+"\n")
logoutf.write("----------------------------------------------"+"\n")

'''
if adapter_file !="":
	print "Adapter to be removed from 3'read: %s"%(adapter.rstrip("\n"));
'''	
#--- Reference genome -------------------------------------------------------------

print "== Reading reference genome ==";

genome_path=db_path  # need full path
genome_file="ref.shelve"

genome_seqs={}
d = shelve.open(genome_path+genome_file,'r')
for chr in d:
	genome_seqs[chr]=d[chr]
	logoutf.write("G ref seq: %s %12d bp"%(chr,len(genome_seqs[chr]))+"\n")
d.close()
logoutf.write("G %d ref sequence(s)"%(len(genome_seqs))+"\n")
logoutf.write("----------------------------------------------"+"\n")


#--- Mapable regions -------------------------------------------------------------
FW_regions={}
RC_regions={}
d2 = shelve.open(genome_path+regions_file,'r')
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
	FW_temp_regions={}
	RC_temp_regions={}
	
d2.close()
#print "# %d mapable regions"%(n_mapable_regions);
logoutf.write("G %d mapable fragments"%(n_mapable_regions)+"\n")
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
start = time.time()

#----------------------------------------------------------------
print "== Start mapping ==";
for read_file in my_files:
	no_my_files+=1
	random_id="-"+main_read_file+".tmp-"+str(random.randint(1000000,9999999))
	outfile1='Trimed_BS.fa'+random_id
	outfile2='Trimed_C2T.fa'+random_id
	
	outf1=open(tmp_path+outfile1,'w')
	outf2=open(tmp_path+outfile2,'w')
	
	#--- Checking input format ------------------------------------------	
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
	
	#----------------------------------------------------------------	
	id=""
	seq=""
	seq_ready=0
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
	bowtie_map1=Popen(b1,shell=True)
	bowtie_map2=Popen(b2,shell=True)
	bowtie_map1.wait()
	bowtie_map2.wait()
	Popen('rm %s%s &'%(tmp_path,outfile2),shell=True)
	
	#print "bowtie",time.time() - start;
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
	#print "unique",time.time() - start;
	
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
	bs_reads_file=outfile1
	
	
	for line in fileinput.input(tmp_path+bs_reads_file):
		if line[0]==">":
			header=line[1:-1]
		else:
			l=line.split()
			original_bs_reads[header]=l[0]
			
	fileinput.close()	
	
	Popen('rm %s%s &'%(tmp_path,outfile1),shell=True)
	#print "read bs fasta",time.time() - start;
	
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
				mapped_location=mapped_location+1 # 1 based (after plus 1)
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
					my_region_serial, my_region_start, my_region_end =my_mapable_region(FW_chr_regions,mapped_location,"+FW") 
				elif FR=="-FW":
					my_region_serial, my_region_start, my_region_end =my_mapable_region(RC_chr_regions,mapped_location+len(original_BS),"-FW")
				#---------------------------------------------
				N_mismatch=N_MIS(original_BS,origin_genome)
				if N_mismatch<= int(indexname) and my_region_serial != 0:
					all_mapped_passed+=1
					#---------------------------------------------
					
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
	
					outf.write('%s	%2d	%3s	%s	%s	%s	%s	%d	%d	%d	%d'%(header,N_mismatch,FR,coordinate,output_genome,original_BS,methy,my_region_serial, my_region_start, my_region_end,STEVE)+"\n")
	
	
	print "--> %s (%d/%d) "%(read_file,no_my_files,len(my_files));
	Popen('rm %s%s %s%s&'%(tmp_path,WC2T, tmp_path,CC2T),shell=True)
	#print "output",time.time() - start;
					
outf.close()
Popen('rm -r %s '%(tmp_path),shell=True)

#----------------------------------------------------------------
		
	
logoutf.write("O Number of raw reads: %d "%(all_raw_reads)+"\n")
if all_raw_reads >0:
	logoutf.write("O Number of CGG/TGG tagged reads: %d (%1.3f)"%(all_tagged,float(all_tagged)/all_raw_reads)+"\n")
	for kk in range(len(n_mytag_lst)):
		logoutf.write("O Number of raw reads with %s tag: %d (%1.3f)"%(mytag_lst[kk],n_mytag_lst[mytag_lst[kk]],float(n_mytag_lst[mytag_lst[kk]])/all_raw_reads)+"\n")
	logoutf.write("O Number of CGG/TGG reads having adapter removed: %d "%(all_tagged_trimed)+"\n")
	logoutf.write("O Number of unique-hits reads for post-filtering: %d"%(all_mapped)+"\n")		

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
elapsed = (time.time() - start)
logoutf.write("T Time: %s"%(elapsed)+"\n")
logoutf.write("------------------- END --------------------"+"\n")
print "=== END ===";

print elapsed;
logoutf.close()
