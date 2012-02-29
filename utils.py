import os

#----------------------------------------------------------------
import datetime
import shutil
import types


# test comment2

def reverse_compl_seq(strseq):
    rc_dict = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G'}
    return ''.join(rc_dict.get(c, c) for c in reversed(strseq.upper()))


def N_MIS(r,g):
    combo=[]
    if len(r)==len(g):
        combo=[r[i]+g[i] for i in range(len(r)) if r[i]!=g[i] and r[i] !="N" and g[i]!="N"]
        combo=[x for x in combo if x !="TC"]
    return len(combo)

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

default_bowtie_path = find_location('bowtie') or "~/bowtie-0.12.7/"

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
