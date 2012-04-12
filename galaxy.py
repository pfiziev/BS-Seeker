import tempfile

__author__ = 'pf'
from subprocess import Popen
#from utils import *
from collections import defaultdict
import sys, shutil, os, re

BUILD = 'build'
BSSEEKER2 = 'bsseeker2'
PATH = 'exec'
PATH_TAG = PATH+'-path'
ARG_TYPES = [BUILD, BSSEEKER2, PATH]

USAGE = """
%(script)s is a wrapper script for bs_seeker2-build.py and bs_seeker2.py that is intended to be used with the Galaxy web platform.
The script takes command line parameters and runs bs_seeker2.py and bs_seeker2-build.py, if neccessary.

The parameters that are related to bs_seeker2-build.py must be prefixed with --%(build_tag)s.
The parameters that are related to bs_seeker2.py must be prefixed with --%(bs_seeker_tag)s.
Additionally, the path to BS Seeker has to be specified via the --%(path_tag)s option.

For example:

    python %(script)s --%(path_tag)s /mnt/Data/UCLA/Matteo/BS-Seeker --build-f data/arabidopsis/genome/Arabidopsis.fa --bsseeker2-i data/arabidopsis/BS6_N1try2L7_seq.txt.fa --bsseeker2-o data/arabidopsis/BS6_N1try2L7_seq.txt.fa.test_output

This will run build the genome in Arabidopsis.fa and put the indexes in a temporary directory. bs_seeker2.py will be run on the
newly created genome index. I.e. the following two commands will be run in a shell:

    python /mnt/Data/UCLA/Matteo/BS-Seeker/bs_seeker2-build.py --db /tmp/tmpg8Eq1o -f /mnt/Data/UCLA/Matteo/bck_BS-Seeker/data/arabidopsis/genome/Arabidopsis.fa

    python /mnt/Data/UCLA/Matteo/BS-Seeker/bs_seeker2.py --db /tmp/tmpg8Eq1o -o /mnt/Data/UCLA/Matteo/bck_BS-Seeker/data/arabidopsis/BS6_N1try2L7_seq.txt.fa.test_output -i /mnt/Data/UCLA/Matteo/bck_BS-Seeker/data/arabidopsis/BS6_N1try2L7_seq.txt.fa -g Arabidopsis.fa


The temporary directory will be deleted after the wrapper exits.


If no options related to bs_seeker2-build are passed, no genome index will be built and the corresponding pre-built genome index will be used
instead. No temporary files and directories will be created.

""" % { 'script' : os.path.split(__file__)[1], 'build_tag' :BUILD, 'bs_seeker_tag' : BSSEEKER2, 'path_tag' : PATH_TAG}


def error(msg):
    print 'ERROR: %s' % msg
    exit(1)


if __name__ == '__main__':
    if len(sys.argv) == 1:
        error('No parameters\n\n'+USAGE)


    # Parse command line arguments
    args = defaultdict(dict)
    arg_key = None
    arg_val = None
    arg_type = None

    for arg in sys.argv[1:]:
        if arg.startswith('--'):
            try:
                arg_type, arg_key = re.match(r'--(\w+)(.*)', arg).groups()
                if arg_type not in ARG_TYPES:
                    raise Exception("Bad argument: %s. arg_type (%s) must be one of: %s." % (arg, arg_type, ', '.join(ARG_TYPES)))
                if not arg_key or arg_key[0] !=  '-':
                    raise Exception("Bad argument: %s. arg_key (%s) must start with - or --." % (arg, arg_key))

            except Exception, e:
                error(str(e) + '\n\n' + USAGE)

            args[arg_type][arg_key] = ''
        else:
            args[arg_type][arg_key] = arg

    path_to_bs_seeker = args.get('exec', {'-path' : None})['-path']
    if path_to_bs_seeker is None:
        error('You have to specify the path to BS Seeker 2 via --%s\n\n' % PATH_TAG + USAGE)

    tempdir = None
    def run_prog(prog, params):
        cwd, _ = os.path.split(__file__)
        cmd = 'python %(prog)s %(params)s' % {
                       'prog'      : os.path.join(cwd, prog),
                       'params'    : ' '.join('%s %s' % (arg_key, arg_val) for arg_key, arg_val in params.items())
                       }
        print 'exec:', cmd

        return_code = Popen(args = cmd, shell = True).wait()
        if return_code:
            if tempdir:
                shutil.rmtree(tempdir)
            error("%s exitted with error code %d" % (prog, return_code))

    if BUILD in args:
        tempdir = tempfile.mkdtemp()
        args[BUILD]['--db'] = tempdir
        args[BSSEEKER2]['--db'] = tempdir
        run_prog(os.path.join(path_to_bs_seeker, 'bs_seeker2-build.py'), args[BUILD])

    run_prog(os.path.join(path_to_bs_seeker, 'bs_seeker2.py'), args['bsseeker2'])

    if tempdir:
        shutil.rmtree(tempdir)