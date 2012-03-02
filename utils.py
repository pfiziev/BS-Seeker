import os

#----------------------------------------------------------------
import datetime
import shutil
import types
from itertools import izip

# test comment2


_rc_dict = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G'}
def reverse_compl_seq(strseq):
    return ''.join(_rc_dict.get(c, c) for c in reversed(strseq.upper()))


def N_MIS(r,g):
    l = 0
    if len(r)==len(g):
        for i in xrange(len(r)):
            if r[i] != g[i] and r[i] != "N" and g[i] != "N" and not(r[i] == 'T' and g[i] == 'C'):
                l += 1
    return i


#----------------------------------------------------------------

def methy_seq(r,g_long):
    H=['A','C','T']
    g_long=g_long.replace("_","")
    m_seq=str()
    xx="-"
    for i in xrange(len(r)):
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
#                      BOWTIE2,
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


def process_aligner_output(filename, format = BOWTIE):
    input = open(filename)
    for l in input:
        if format == BOWTIE:
            yield l
        elif format == BOWTIE2:
            QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL = l.split()

            # skip reads that are not mapped
            if FLAG & 0x4:
                continue

            yield '\t'.join([
                             QNAME, # read ID
                             RNAME, # reference ID
                             POS, # position

                            ])
    pass


def cigar_to_alignment(cigar_string, read_seq, genome_seq, genome_start):
    """ Reconstruct the pairwise alignment based on the CIGAR string and the two sequences
    """
    # first, parse the CIGAR string
    i = 0
    prev_i = 0
    cigar = []
    tags = {'S', 'M', 'D', 'I'}

    while i < len(cigar_string):
        if cigar_string[i] in tags:
            cigar.append((cigar_string[i], int(cigar_string[prev_i:i])))
            prev_i = i + 1
        i += 1

    # reconstruct the alignment
    r_pos = cigar[0][1] if cigar[0][0] == 'S' else 0
    g_pos = genome_start
    r_aln, g_aln = '', ''
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
    """ Deletes a number of files. filenames can contain generator expressions and/or lists, too"""

    for fname in filenames:
        if type(fname) in [list, types.GeneratorType]:
            delete_files(*list(fname))
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
