def hxb2_translate(hxb2_seq):
    ''' return an array of strings such that xlate[n] is the HXB2 name for the n'th site in the alignment
    here, n is zero-based position, but the HXB2 position is one based.  To get a sense of what this
    looks like: HXB2 = "mrv-k--eky..." -> xlate = [001,002,003,003a,004,004a,004b,005,006,007]"
    '''
    abc="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
    xlate = []
    n=0
    nabc=0
    for c in hxb2_seq:
        if c != '-':
            n += 1
            nabc = 0
            xlate.append( "%03d" % (n,) )
        else:
            xlate.append( "%03d%1s" % (n,abc[nabc]) )
            nabc += 1
    return xlate        

class hxb2:

    def __init__(self,seq):
        self.seq = seq
        self.xlate = hxb2_translate(seq)

    def seq(self):
        return self.seq

    def fname(self,n,c):
        name=""
        if self.seq[n] != "-":
            name += self.seq[n]
        name += self.xlate[n]
        name += c
        return name
