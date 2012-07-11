import sys


import pysam

def bs_seeker(outf):
    for a in pysam.Samfile(outf, 'rb'):
        _, original_strand, ct, rno, strand, start, end  = a.qname.split('_')
        yield a.qname, original_strand, strand, start, end, a.pos

def bismark(outf):

    for l in open(outf):
        if l[0] == '@':
            continue
        buf = l.split()
        _, original_strand, ct, rno, strand, start, end  = buf[0].split('_')
        yield  buf[0], original_strand, strand, start, end, buf[3]





def rmapbs(outf):
    pass


if __name__ == '__main__':


    if len(sys.argv) != 3:
        print "usage: %s bs_seeker|bismark|rmapbs output-file" % __file__
        exit(1)

    total = 0
    tp = 0
    if sys.argv[1] == 'bs_seeker':
        parse_output = bs_seeker
    elif sys.argv[1] == 'bismark':
        parse_output = bismark
    elif sys.argv[1] == 'rmapbs':
        parse_output = rmapbs
    else:
        print "usage: %s bs_seeker|bismark|rmapbs output-file" % __file__
        exit(1)

    chr22_len = 49691432

    for rid, original_strand, strand, start, end, mapped_pos in parse_output(sys.argv[2]):
        start, end, mapped_pos = map(int, [start, end, mapped_pos])
        if original_strand == 'rev':
            start = chr22_len - start
            end = chr22_len - end
            start, end = end, start
        if abs(start - mapped_pos) < 20:
            tp += 1
#        else:
#            print rid
        total += 1
    print 'total:', total, 'tp:', tp, '(%.2lf)' % (100*float(tp)/total)

