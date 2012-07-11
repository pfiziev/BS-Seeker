import json

__author__ = 'pf'
import random
import string

_rc_trans = string.maketrans('ACGT', 'TGCA')
def reverse_compl_seq(strseq):
    return strseq.translate(_rc_trans)[::-1]


def generate_reads(out, cell_type, seq, nreads):
    nucs = 'ACTG'
    
    for i in xrange(nreads):
        
        read_seq = None
        while read_seq is None or 'N' in read_seq:
            start = random.randint(0, len(seq) - READ_LENGTH)
            end = start + READ_LENGTH
            read_seq = seq[start:end]
        
        strand = '+'
        if random.random() < 0.5:
            read_seq = reverse_compl_seq(read_seq)
            strand = '-'
            
        read_seq = ''.join(n if random.random() > ERROR_RATE else (random.choice([k for k in nucs if k != n])) for n in read_seq)
            
        
        out.write('>read_%s_%d_%s_%d_%d\n%s\n' % (cell_type, i, strand, start, end, read_seq))
        
        


if __name__ == '__main__':

    METHYLATION_RATE = 0.3
    COVERAGE = 3
    CELL_TYPES = 4
    READ_LENGTH = 75
    ERROR_RATE = 0.05

    genome = ''.join([l.strip() for l in list(open('chr22.fa')) if l[0] != '>']).upper()
    
    nreads = COVERAGE*len(genome)/READ_LENGTH

    out = open('chr22_reads.fa','w')
    print "reading genome done"
    methylation = {}


    def methylate(seq, is_fwd_strand):
        methylated_nucs = []
        for i in xrange(len(seq)):
            if seq[i] != 'C':
                methylated_nucs.append(seq[i])
            else:
                pos = i if is_fwd_strand else len(seq) - i

                if pos not in methylation:
                    methylation[pos] = {'C' : 0, 'T' : 0}

                if random.random() < METHYLATION_RATE:
                    methylated_nucs.append('C')
                    methylation[pos]['C'] += 1
                else:
                    methylated_nucs.append('T')
                    methylation[pos]['T'] += 1

        return ''.join(methylated_nucs)

    for cell_type in xrange(CELL_TYPES):
        print cell_type

        generate_reads( out,
                        'fwd_%d' % cell_type,
                        methylate(genome, True),
                        nreads/2 )
                        
        generate_reads( out,
                        'rev_%d' % cell_type,
                        methylate(reverse_compl_seq(genome), False),
                        nreads/2 )

    for i in methylation:
        methylation[i]['level'] = float(methylation[i]['C'])/(methylation[i]['C'] + methylation[i]['T']) if (methylation[i]['C'] + methylation[i]['T']) > 0 else 0
    json.dump(methylation.items(), open('chr22_methylation.json','w'))
        
    
