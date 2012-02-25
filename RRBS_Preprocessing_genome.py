#----------------------------------------------------------------
#----------------------------------------------------------------




# Part 1



#----------------------------------------------------------------
#----------------------------------------------------------------
import fileinput, string,os, gzip,copy, time, shelve,subprocess, random, math
from subprocess import Popen
	
#--------------------------------------------------------------------------------
from optparse import OptionParser

parser = OptionParser()

parser.set_defaults(low_b=75)
parser.add_option("-l", "--low", dest="low_b",help="lower bound")
parser.set_defaults(up_b=280)
parser.add_option("-u", "--up", dest="up_b",help="upper bound")


parser.set_defaults(filename="mouse_mm9.fa")
parser.add_option("-f", "--file", dest="filename",help="Input your reference genome file (fasta format)", metavar="FILE")

#parser.set_defaults(bowtiepath="~/bowtie-0.12.7/")
parser.set_defaults(bowtiepath="/u/home/mcdb/paoyang/bowtie-0.12.7/")
#parser.set_defaults(bowtiepath="/Users/Shared/LAB-Applications/Bowtie-0.12.7-OSX10.4quadUniversalSJC/")

parser.add_option("-p", "--path", dest="bowtiepath",help="Path to Bowtie [~/bowtie-0.10.0/]", metavar="PATH")


(options, args) = parser.parse_args()
                  

#--------------------------------------------------------------------------------
low_bound=options.low_b
low_bound=int(low_bound)

up_bound=options.up_b
up_bound=int(up_bound)

fasta_file=options.filename
fasta_file=str(fasta_file)

frag_range=str(low_bound)+"_"+str(up_bound)

outfile="RRBS_mapable_genome.fa"
mapfile="RRBS_mapable_genome.shelve"

mapable_region_shelve="RRBS_mapable_regions.shelve"
mapable_region_file="RRBS_mapable_regions.txt"

start = time.time()

os.system('mkdir ./RRBS_reference_genome_%s'%(frag_range))
ref_path='./RRBS_reference_genome_%s/'%(frag_range)
genome_path="./"
logfile="log_RRBS_MAPABLE_"+fasta_file[:-3]+"_"+str(low_bound)+"_"+str(up_bound)+".txt"
ref_log=open(ref_path+logfile,"w")

bowtie_path=options.bowtiepath
bowtie_path=str(bowtie_path)
if bowtie_path[-1] !="/":
	bowtie_path=bowtie_path+"/"
#--------------------------------------------------------------------------------

ref_log.write("---------------- Extracting mapable regions and mapable genome (1-2) ----------------"+"\n")
ref_log.write("\n")

#--- Read-in Reference genome ----------------------------------------------------------
start = time.time()

g=""
header=""
genome={}
all_base=0
n=0
for line in fileinput.input(genome_path+fasta_file):
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
g=g.upper()
ref_log.write("reference seq: %s (renamed as %s ) %d bp"%(header,short_header,len(g))+"\n")
genome[short_header]=g
all_base+=len(g)
ref_log.write("--- In tatal %d reference seqs ==> %d bp"%(len(genome),all_base)+"\n")


#--- positive strand -------------------------------------------------------------
outf = open(ref_path+outfile,"w")
d2   = shelve.open(ref_path+mapfile,'n')
d3   = shelve.open(ref_path+mapable_region_shelve,'n')
f3   = open(ref_path+mapable_region_file,"w")

all_L=0
all_mappable_length=0
all_unmappable_length=0
all_chr=genome.keys()
all_chr.sort()
seq=""

no_mapable_region=0
for chr in all_chr:
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
		DD=(CCGG_sites[j+1]-CCGG_sites[j])-4 # NOT including both CCGG
		if DD>=low_bound and DD<=up_bound:
			CCGG_CCGG.append([CCGG_sites[j],CCGG_sites[j+1]+3]) # leftmost <--> rightmost
			mapable_seq=seq[CCGG_sites[j]:CCGG_sites[j+1]+4]
			no_mapable_region+=1
			
			
			chr_regions.append([CCGG_sites[j],CCGG_sites[j+1]+3,no_mapable_region,mapable_seq])
			# start_position, end_position, serial, sequence
			f3.write("%s\t%d\t%d\t%d\t%s"%(chr,no_mapable_region,CCGG_sites[j],CCGG_sites[j+1]+3,mapable_seq)+"\n")
			
	CCGG_sites=[]
	d3[chr]=chr_regions
	chr_regions=[]
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
			if m>=p1 and m<p2+1:
				map_seq+=seq[m]
				mappable_length+=1
				mark="open"
			else:
				if mark=="close":
					map_seq+="-"
					unmappable_length+=1
				elif mark=="open": # the last eligible base
					dump=CCGG_CCGG.pop(0)
					if len(CCGG_CCGG)>0:
						pair=CCGG_CCGG[0]
						p1=pair[0]
						p2=pair[1]
						if m>=p1 and m<p2+1:
							map_seq+=seq[m]
							mappable_length+=1
							mark="open"
						else:
							map_seq+="-"
							unmappable_length+=1
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
		
	CCGG_CCGG=[]
	seq=""
	del genome[chr];
	#-----------------------------------
	d2[chr]=map_seq
	outf.write(">%s"%(chr)+"\n")
	for ii in range(0,L,50):
		y=min(ii+50,L)
		outf.write("%s"%(map_seq[ii:y])+"\n")
	
	#-----------------------------------
	map_seq=""
	ref_log.write("# %s : all %d : %d (unmapable -) %d (mapable) (%1.5f)"%(chr,L,unmappable_length,mappable_length,float(mappable_length)/L)+"\n")
	all_L+=L
	all_mappable_length+=mappable_length
	all_unmappable_length+=unmappable_length
	
ref_log.write("# total %d chromosome seqs ==> %d : %d (unmapable -) %d (mapable) (%1.5f)" %(len(all_chr),all_L,all_unmappable_length,all_mappable_length,float(all_mappable_length)/all_L)+"\n")

ref_log.write("#       %d eligible fragments"%(no_mapable_region)+"\n")


d2.close()
d3.close()
outf.close()
Popen('nohup gzip %s &'%(ref_path+outfile),shell=True)

#----------------------------------------------------------------
#----------------------------------------------------------------




# Part 2



#----------------------------------------------------------------
ref_log.write("\n")
ref_log.write("----------------         Pre-processing mapable genome         (2-2) ----------------"+"\n")
#----------------------------------------------------------------
def reverse_compl_seq(strseq):
	strseq=strseq.upper()
	rc_strseq=strseq.translate(string.maketrans("ATCG", "TAGC"))[::-1]
	return rc_strseq;
	
#----------------------------------------------------------------
parser = OptionParser()

fasta_file="RRBS_mapable_genome.fa"


print "Reference genome file: %s"%(fasta_file);
print "Bowtie path: %s"%(bowtie_path);
#---------------------------------------------------------------

#os.system('mkdir ./reference_genome')
#ref_path='./RRBS_reference_genome_%s/'
#---------------------------------------------------------------
# 1. First get the complementary genome (also do the reverse)
# 2. Then do CT and GA conversions
#---------------------------------------------------------------
FW_genome={}
header=""
g=''
n=0

#ref_log=open(ref_path+"log_PREPROCESSING_GNM_"+fasta_file+".txt","w")
refd = shelve.open(ref_path+"refname.shelve",'n')

for line in fileinput.input(ref_path+fasta_file):
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
			
g=g.upper()
short_header=str(n).zfill(4)
#print "reference seq: %s (renamed as %s) %d bp"%(header,short_header,len(g));
ref_log.write("Pre-processing reference seq: %s ( %d bp)"%(short_header,len(g))+"\n")
refd[short_header]=[header,len(g)]
FW_genome[short_header]=g
g=""
refd.close()

FW_lst=FW_genome.keys()
FW_lst.sort()

#---------------- Python shleve -----------------------------------------------
d = shelve.open(ref_path+"ref.shelve",'n')
for chr_id in FW_genome:
	d[chr_id]=FW_genome[chr_id]
d.close()
#---------------- Reverse complement (Crick strand) ----------------------------
header=""
RC_genome={}
for header in FW_lst:
	g=FW_genome[header]
	g=reverse_compl_seq(g)
	RC_genome[header]=g
RC_lst=RC_genome.keys()
RC_lst.sort()



#---------------- 4 converted fasta -------------------------------------------

outf=open(ref_path+'W_C2T.fa','w')
for header in FW_lst:
	outf.write('>%s'%(header)+"\n")
	g=FW_genome[header]
	g=g.replace("c","t")
	g=g.replace("C","T")
	outf.write('%s'%(g)+'\n')
outf.close()
print 'end 4-1';

outf=open(ref_path+'C_C2T.fa','w')
for header in RC_lst:
	outf.write('>%s'%(header)+"\n")
	g=RC_genome[header]
	g=g.replace("c","t")
	g=g.replace("C","T")
	outf.write('%s'%(g)+'\n')
outf.close()
print 'end 4-2';

outf=open(ref_path+'W_G2A.fa','w')
for header in FW_lst:
	outf.write('>%s'%(header)+"\n")
	g=FW_genome[header]
	g=g.replace("g","a")
	g=g.replace("G","A")
	outf.write('%s'%(g)+'\n')
outf.close()
print 'end 4-3';
FW_lst={}

outf=open(ref_path+'C_G2A.fa','w')
for header in RC_lst:
	outf.write('>%s'%(header)+"\n")
	g=RC_genome[header]
	g=g.replace("g","a")
	g=g.replace("G","A")
	outf.write('%s'%(g)+'\n')
outf.close()
print 'end 4-4';
RC_lst={}
#---------------- bowtie libraries -------------------------------------------
map1=Popen('nohup %sbowtie-build -f ./RRBS_reference_genome_%s/W_C2T.fa ./RRBS_reference_genome_%s/W_C2T > ./RRBS_reference_genome_%s/W_C2T.log'%(bowtie_path,frag_range,frag_range,frag_range),shell=True)
map2=Popen('nohup %sbowtie-build -f ./RRBS_reference_genome_%s/W_G2A.fa ./RRBS_reference_genome_%s/W_G2A > ./RRBS_reference_genome_%s/W_G2A.log'%(bowtie_path,frag_range,frag_range,frag_range),shell=True)
map3=Popen('nohup %sbowtie-build -f ./RRBS_reference_genome_%s/C_C2T.fa ./RRBS_reference_genome_%s/C_C2T > ./RRBS_reference_genome_%s/C_C2T.log'%(bowtie_path,frag_range,frag_range,frag_range),shell=True)
map4=Popen('nohup %sbowtie-build -f ./RRBS_reference_genome_%s/C_G2A.fa ./RRBS_reference_genome_%s/C_G2A > ./RRBS_reference_genome_%s/C_G2A.log'%(bowtie_path,frag_range,frag_range,frag_range),shell=True)
map1.wait()
map2.wait()
map3.wait()
map4.wait()
os.system('rm ./RRBS_reference_genome_%s/W_C2T.fa ./RRBS_reference_genome_%s/W_G2A.fa ./RRBS_reference_genome_%s/C_C2T.fa ./RRBS_reference_genome_%s/C_G2A.fa'%(frag_range,frag_range,frag_range,frag_range))

elapsed = (time.time() - start)
ref_log.write("Time: %s"%(elapsed)+"\n")
ref_log.close()
print ("Time: %s"%(elapsed));
print 'END';

