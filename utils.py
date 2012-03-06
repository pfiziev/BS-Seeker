import os

#----------------------------------------------------------------
import datetime
import re
import shutil
import types
from itertools import izip

# test comment2


_rc_dict = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G'}
def reverse_compl_seq(strseq):
    return ''.join(_rc_dict.get(c, c) for c in reversed(strseq.upper()))


def N_MIS(r,g):
    mismatches = 0
    if len(r)==len(g):
        for i in xrange(len(r)):
            if r[i] != g[i] and r[i] != "N" and g[i] != "N" and not(r[i] == 'T' and g[i] == 'C'):
                mismatches += 1
    return mismatches


#----------------------------------------------------------------

def next_nuc(seq, pos, n):
    """ Returns the nucleotide that is n places from pos in seq. Skips gap symbols.
    """
    i = pos + 1
    while i < len(seq):
        if seq[i] != '-':
            n -= 1
            if n == 0: break
        i += 1
    return seq[i]



def methy_seq(read, genome):
    H = ['A', 'C', 'T']
    m_seq = ''
    xx = "-"
    for i in xrange(len(read)):

        if genome[i] == '-':
            continue

        elif read[i] != 'C' and read[i] != 'T':
            xx = "-"

        elif read[i] == "T" and genome[i] == "C": #(unmethylated):
            nn1 = next_nuc(genome, i, 1)
            if nn1 == "G":
                xx = "x"
            elif nn1 in H :
                nn2 = next_nuc(genome, i, 2)
                if nn2 == "G":
                    xx = "y"
                elif nn2 in H :
                    xx = "z"

        elif read[i] == "C" and genome[i] == "C": #(methylated):
            nn1 = next_nuc(genome, i, 1)

            if nn1 == "G":
                xx = "X"
            elif nn1 in H :
                nn2 = next_nuc(genome, i, 2)

                if nn2 == "G":
                    xx = "Y"
                elif nn2 in H:
                    xx = "Z"
        else:
            xx = "-"
        m_seq += xx

    return m_seq

def _methy_seq(r, g_long):
    H = ['A', 'C', 'T']
    g_long = g_long.replace("_", "")
    m_seq = ''
    xx = "-"
    for i in xrange(len(r)):
        if r[i] not in ['C', 'T']:
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



def mcounts(mseq,mlst,ulst):
    out_mlst=[mlst[0]+mseq.count("X"),mlst[1]+mseq.count("Y"),mlst[2]+mseq.count("Z")]
    out_ulst=[ulst[0]+mseq.count("x"),ulst[1]+mseq.count("y"),ulst[2]+mseq.count("z")]
    return out_mlst,out_ulst

#-------------------------------------------------------------------------------------

# set a reasonable defaults
def find_location(program):
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)
    for path in os.environ["PATH"].split(os.pathsep):
        if is_exe(os.path.join(path, program)):
            return path

    return None

BOWTIE = 'bowtie'
BOWTIE2 = 'bowtie2'
SOAP = 'soap'

supported_aligners = [
                      BOWTIE,
                      BOWTIE2,
                     # SOAP
                    ]

aligner_options_prefixes = { BOWTIE : '--bt-',
                             BOWTIE2 : '--bt2-',
                             SOAP   : '--soap-' }
aligner_path = dict((aligner, find_location(aligner) or default_path) for aligner, default_path in [(BOWTIE,'~/bowtie-0.12.7/'),
                                                                                                    (BOWTIE2, '~/bowtie-0.12.7/'),
#                                                                                                    (SOAP, '~/soap2.21release/')
                                                                                                    ])
#
#default_bowtie_path = find_location('bowtie') or "~/bowtie-0.12.7/"
#default_bowtie2_path = find_location('bowtie2') or "~/bowtie-0.12.7/"
#default_soap_path = find_location('soap') or "~/soap2.21release/"


def process_aligner_output(filename, pair_end = False):

    QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL = range(11)

    m = re.search(r'-('+'|'.join(supported_aligners) +')-TMP', filename)
    if m is None:
        error('The temporary folder path should contain the name of one of the supported aligners: ' + filename)

    format = m.group(1)

    input = open(filename)
    def parse_bowtie(line):
        l = line.split()

        header = l[0]
        chr = l[1]
        location = int(l[2])
        #-------- mismatches -----------
        if len(l) == 4:
            no_mismatch = 0
        elif len(l) == 5:
            no_mismatch = l[4].count(":")
        else:
            print l
        return header, chr, location, no_mismatch, None

    def parse_SAM(line):
        buf = line.split()

        flag = int(buf[FLAG])
        # skip reads that are not mapped
        if flag & 0x4:
            return None, None, None, None, None, None

        mismatches = int([buf[i][5:] for i in xrange(11, len(buf)) if buf[i][:5] == 'NM:i:'][0]) # get the edit distance

        # add the soft clipped nucleotides to the number of mismatches
        cigar_string = buf[CIGAR]
        cigar = parse_cigar(cigar_string)
        if cigar[0][0] == 'S':
            mismatches += cigar[0][1]
        if cigar[-1][0] == 'S':
            mismatches += cigar[-1][1]

        return (buf[QNAME], # read ID
                buf[RNAME], # reference ID
                int(buf[POS]) - 1, # position, 0 based (SAM is 1 based)
                mismatches,    # number of mismatches
                cigar_string, # the cigar string
                flag & 0x40 # true if it is the first mate in a pair, false if it is the second mate
                )


    if format == BOWTIE:
        if pair_end:
            for line in input:
                header1, chr, location1, no_mismatch1, _ = parse_bowtie(line)
                header2,   _, location2, no_mismatch2, _ = parse_bowtie(input.next())

                # flip the location info if the second pair comes first in the alignment file
                if header1[-1] == '2':
                    location1, location2 = location2, location1

                yield header1[:-2], chr, no_mismatch1 + no_mismatch2, location1, None, location2, None
        else:
            # single end
            for line in input:
                yield parse_bowtie(line)

    elif format == BOWTIE2:
        if pair_end:
            for line in input:
                header1, chr1, location1, no_mismatch1, cigar_string1,        _ = parse_SAM(line)
                header2,    _, location2, no_mismatch2, cigar_string2, mate_no2 = parse_SAM(input.next())


                if header1 and header2:
                    # flip the location info if the second mate comes first in the alignment file
                    if mate_no2:
                        location1, location2 = location2, location1
                        cigar_string1, cigar_string2 = cigar_string2, cigar_string1


                    yield header1, chr1, no_mismatch1 + no_mismatch2, location1, cigar_string1, location2, cigar_string2
        else:
            for line in input:
                header, chr, location, no_mismatch, cigar_string, _ = parse_SAM(line)
                if header is not None:
                    yield header, chr, location, no_mismatch, cigar_string

    input.close()


def parse_cigar(cigar_string):
    i = 0
    prev_i = 0
    cigar = []
    tags = {'S', 'M', 'D', 'I'}

    while i < len(cigar_string):
        if cigar_string[i] in tags:
            cigar.append((cigar_string[i], int(cigar_string[prev_i:i])))
            prev_i = i + 1
        i += 1
    return cigar


def reverse_cigar_string(cigar_string):
    return ''.join('%d%s'%(count, tag) for tag, count in reversed(parse_cigar(cigar_string)))

def get_read_start_end_and_genome_length(cigar_string):
    cigar = parse_cigar(cigar_string)
    r_start = cigar[0][1] if cigar[0][0] == 'S' else 0
    r_end = r_start
    g_len = 0
    for edit_op, count in cigar:
        if edit_op == 'M':
            r_end += count
            g_len += count
        elif edit_op == 'I':
            r_end += count
        elif edit_op == 'D':
            g_len += count
    return r_start, r_end, g_len # return the start and end in the read and the length of the genomic sequence



def cigar_to_alignment(cigar_string, read_seq, genome_seq):
    """ Reconstruct the pairwise alignment based on the CIGAR string and the two sequences
    """
    # first, parse the CIGAR string
    cigar = parse_cigar(cigar_string)

    # reconstruct the alignment
    r_pos = cigar[0][1] if cigar[0][0] == 'S' else 0
    g_pos = 0
    r_aln = ''
    g_aln = ''
    for edit_op, count in cigar:
        if edit_op == 'M':
            r_aln += read_seq[r_pos : r_pos + count]
            g_aln += genome_seq[g_pos : g_pos + count]
            r_pos += count
            g_pos += count
        elif edit_op == 'D':
            r_aln += '-'*count
            g_aln += genome_seq[g_pos : g_pos + count]
            g_pos += count
        elif edit_op == 'I':
            r_aln += read_seq[r_pos : r_pos + count]
            g_aln += '-'*count
            r_pos += count

    return r_aln, g_aln




reference_genome_path = os.path.join(os.path.split(globals()['__file__'])[0],'reference_genomes')



def error(msg):
    print 'ERROR: %s' % msg
    exit(1)


global_stime = datetime.datetime.now()
def elapsed(msg = None):
    print "[%s]" % msg if msg is not None else "+", "Last:" , datetime.datetime.now() - elapsed.stime, '\tTotal:', datetime.datetime.now() - global_stime
    elapsed.stime = datetime.datetime.now()

elapsed.stime = datetime.datetime.now()

def clear_dir(path):
    """ If path does not exist, it creates a new directory.
        If path points to a directory, it deletes all of its content.
        If path points to a file, it raises an exception."""

    if os.path.exists(path):
        if not os.path.isdir(path):
            error("%s is a file. Please, delete it manually!" % path)
        else:
            for the_file in os.listdir(path):
                file_path = os.path.join(path, the_file)
                try:
                    if os.path.isfile(file_path):
                        os.unlink(file_path)
                    elif os.path.isdir(file_path):
                        shutil.rmtree(file_path)
                except Exception, e:
                    print e
    else:
        os.mkdir(path)


def delete_files(*filenames):
    return
    """ Deletes a number of files. filenames can contain generator expressions and/or lists, too"""

    for fname in filenames:
        if type(fname) in [list, types.GeneratorType]:
            delete_files(*list(fname))
        elif os.path.isdir(fname):
            shutil.rmtree(fname)
        else:
            os.remove(fname)

def split_file(filename, output_prefix, nlines):
    """ Splits a file (equivalend to UNIX split -l ) """
    fno = 0
    lno = 0
    input = open(filename, 'r')
    output = None
    for l in input:
        if lno == 0:
            fno += 1
            if output is not None: output.close()
            output = open('%s%d' % (output_prefix, fno), 'w')
            lno = nlines
        output.write(l)
        lno -= 1
    output.close()
    input.close()


def detect_format(filename):
    read_inf = open(filename, "r")
    oneline = read_inf.readline()
    l = oneline.split()
    input_format = ""

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
