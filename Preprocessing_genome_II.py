import fileinput, string,os, operator, shelve, time, subprocess
from subprocess import Popen
from optparse import OptionParser
from utils import *

'''
#----------------------------------------------------------------
def reverse_compl_seq(strseq):
	strseq=strseq.upper()
	rc_strseq=strseq.translate(string.maketrans("ATCG", "TAGC"))[::-1]
	return rc_strseq;
	
#----------------------------------------------------------------
'''

if __name__ == '__main__':

	parser = OptionParser()
	parser.add_option("-f", "--file", dest="filename",help="Input your reference genome file (fasta)", metavar="FILE")

	parser.set_defaults(taginfo="Y")
	parser.add_option("-t", "--tag", dest="taginfo",help="Yes for undirectional lib, no for directional [Y]", metavar="TAG")

	parser.set_defaults(bowtiepath = default_bowtie_path)
	parser.add_option("-p", "--path", dest="bowtiepath",help="Path to Bowtie [~/bowtie-0.10.0/]", metavar="PATH")

	(options, args) = parser.parse_args()

	# if no options were given by the user, print help and exit
	import sys
	if len(sys.argv) == 1:
		print parser.print_help()
		exit(0)



	fasta_file=options.filename
	fasta_file=str(fasta_file)

	asktag=options.taginfo
	asktag=str(asktag)
	asktag=asktag.upper()

	bowtie_path=options.bowtiepath
	bowtie_path=str(bowtie_path)
	if bowtie_path[-1] !="/":
		bowtie_path=bowtie_path+"/"

	print "Reference genome file: %s"%(fasta_file);
	print "BS reads from undirectional/directional library: %s"%(asktag);
	print "Bowtie path: %s"%(bowtie_path);
	#---------------------------------------------------------------
	start = time.time()

	os.system('mkdir ./reference_genome')
	ref_path='./reference_genome/'
	#---------------------------------------------------------------
	# 1. First get the complementary genome (also do the reverse)
	# 2. Then do CT and GA conversions
	#---------------------------------------------------------------
	FW_genome={}
	header=""
	g=''
	n=0

	ref_log=open(ref_path+fasta_file+".log","w")
	refd = shelve.open(ref_path+"refname.shelve",'n')

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
				print "reference seq: %s (renamed as %s ) %d bp"%(header,short_header,len(g));
				ref_log.write("reference seq: %s (renamed as %s ) %d bp"%(header,short_header,len(g))+"\n")
				refd[short_header]=[header,len(g)]
				FW_genome[short_header]=g

				g=""
				header=l[0][1:]
				n+=1
				short_header=str(n).zfill(4)
				
	g=g.upper()
	short_header=str(n).zfill(4)
	print "reference seq: %s (renamed as %s) %d bp"%(header,short_header,len(g));
	ref_log.write("reference seq: %s (renamed as %s ) %d bp"%(header,short_header,len(g))+"\n")
	refd[short_header]=[header,len(g)]
	FW_genome[short_header]=g
	g=""
	refd.close()
	ref_log.close()

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


	if asktag=="Y":
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
		map1=Popen('nohup %sbowtie-build -f ./reference_genome/W_C2T.fa ./reference_genome/W_C2T > ./reference_genome/W_C2T.log'%(bowtie_path),shell=True)
		map2=Popen('nohup %sbowtie-build -f ./reference_genome/W_G2A.fa ./reference_genome/W_G2A > ./reference_genome/W_G2A.log'%(bowtie_path),shell=True)
		map3=Popen('nohup %sbowtie-build -f ./reference_genome/C_C2T.fa ./reference_genome/C_C2T > ./reference_genome/C_C2T.log'%(bowtie_path),shell=True)
		map4=Popen('nohup %sbowtie-build -f ./reference_genome/C_G2A.fa ./reference_genome/C_G2A > ./reference_genome/C_G2A.log'%(bowtie_path),shell=True)
		map1.wait()
		map2.wait()
		map3.wait()
		map4.wait()
		os.system('rm ./reference_genome/W_C2T.fa ./reference_genome/W_G2A.fa ./reference_genome/C_C2T.fa ./reference_genome/C_G2A.fa')
		
	elif asktag=="N":
		#---------------- 2 converted fasta -------------------------------------------
		
		outf=open(ref_path+'W_C2T.fa','w')
		for header in FW_lst:
			outf.write('>%s'%(header)+"\n")
			g=FW_genome[header]
			g=g.replace("c","t")
			g=g.replace("C","T")
			outf.write('%s'%(g)+'\n')
		outf.close()
		print 'end 2-1';
		FW_lst={}
		
		outf=open(ref_path+'C_C2T.fa','w')
		for header in RC_lst:
			outf.write('>%s'%(header)+"\n")
			g=RC_genome[header]
			g=g.replace("c","t")
			g=g.replace("C","T")
			outf.write('%s'%(g)+'\n')
		outf.close()
		print 'end 2-2';
		
		RC_lst={}
		#---------------- bowtie libraries -------------------------------------------
		map1=Popen('nohup %sbowtie-build -f ./reference_genome/W_C2T.fa ./reference_genome/W_C2T > ./reference_genome/W_C2T.log'%(bowtie_path),shell=True)
		map3=Popen('nohup %sbowtie-build -f ./reference_genome/C_C2T.fa ./reference_genome/C_C2T > ./reference_genome/C_C2T.log'%(bowtie_path),shell=True)
		map1.wait()
		map3.wait()
		os.system('rm ./reference_genome/W_C2T.fa ./reference_genome/C_C2T.fa')	
	elapsed = (time.time() - start)
	print ("Time: %s"%(elapsed));
	print 'END';

