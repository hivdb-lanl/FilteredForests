## sequence utilities

import config
import re
from collections import Counter
from SequenceSample import SequenceSample

def find_alldash_columns(seqs):
    alldashcol=[]
    maxstrlen = max( [len(seq.str) for seq in seqs] )
    for n in range(maxstrlen):
        alldash=True
        for m in range(len(seqs)):
            if len(seqs[m].str)>n and seqs[m].str[n] != '-':
                alldash=False
                break
        if alldash:
            alldashcol.append(n)
    return alldashcol

def read_tbl(filename,rmdash=False):
    '''
    read_tbl: read .tbl file, return list of sequences
    rmdash=False: remove "-" from string; useful for unaligned
    '''
    re_comment    = re.compile('^\#.*')
    re_leadwhite  = re.compile('^\s*')
    re_trailwhite = re.compile('\s*$')
    re_badchar    = re.compile('[\#\$\*X]')
    re_dash       = re.compile('-')

    seq_samples=[]
    if not filename:
        return seq_samples

    with open(filename,'rU') as tbl:
        for line in tbl.readlines():
            line = re_comment.sub('',line)
            line = re_leadwhite.sub('',line)
            line = re_trailwhite.sub('',line)
            if not line:
                continue
            tokens = line.split()
            if len(tokens) != 2:
                print "Invalid line in tbl file:"
                print "[",line,"]"
                continue
            (name,str) = tokens[:2]
            if str:
                str = re_badchar.sub('x',str)
                if rmdash:
                    str = re_dash.sub('',str)
                sample = SequenceSample(name,str)
                seq_samples.append(sample)

    return seq_samples

def write_tbl(seqlist,filename,names=None):
    with open(filename,"w") as h:
        for s in seqlist:
            if names and s.name not in names:
                continue
            h.write("%s %s\n" % (s.name,s.str))

def write_mase(seqlist,filename,names=None):
    WIDTH=70
    with open(filename,"w") as h:
        for s in seqlist:
            if names and s.name not in names:
                continue
            h.write(";\n%s\n" % s.name)
            str = s.str
            while len(str)>0:
                if len(str)>WIDTH:
                    h.write("%s\n" % str[:WIDTH])
                    str = str[WIDTH:]
                else:
                    h.write("%s\n" % str)
                    str = ""

def read_mase(filename):
    re_comment = re.compile('^;.*')
    re_leadwhite  = re.compile('^\s*')
    re_trailwhite = re.compile('\s*$')
    re_white = re.compile('\s*')

    seq_samples=[]
    cur_name=''
    cur_str=''
    with open(filename,'rU') as mase:
        for line in mase.readlines():
            line = re_leadwhite.sub('',line)
            line = re_trailwhite.sub('',line)
            ## This should not be necessary!
            if re_white.match(line):
                print "Ignoring whitespace in MASE file: ",filename
                line = re_white.sub('',line)
            if cur_name:
                if re_comment.match(line):
                    sample = SequenceSample(cur_name,cur_str)
                    seq_samples.append(sample)
                    cur_name=''
                    cur_str=''
                else:
                    cur_str += line
            else:
                line = re_comment.sub('',line)
                if not line:
                    continue
                cur_name=line
    if cur_name:
        ## append the last sequence in the file
        sample = SequenceSample(cur_name,cur_str)
        seq_samples.append(sample)
    return seq_samples

def seq_filter(seqs,reverse=False,tmr=0,R=0,rmx=True):
    '''
    reverse=False: reverse all the strings
    tmr=0: filter out strings with too many repeats 
    R=0: if R=(lo,hi) then truncate strings to this range
    rmx: filter out strings that have no epitopes 
    '''
    if reverse:
        for s in seqs:
            s.str = s.str[::-1]
    if tmr:
        seqs = seq_filter_repeats(seqs,tmr)

    if R and len(R)==2:
        if config.VERBOSE:
            print "Truncate to range:",R
        seqs = truncate_strings(seqs,R)

    if rmx:
        badseqs=[]
        for s in seqs:
            e = episeq_from_string(s.str,rmx=True)
            if len(e)==0:
                if config.VERBOSE:
                    print "No epitopes, so ignoring seq:",s.name
                badseqs.append(s)
        if badseqs:
            print "Ignoring",len(badseqs),"bad sequences"
            newseqs=[]
            for s in seqs:
                if s not in badseqs:
                    newseqs.append(s)
            seqs = newseqs

    return seqs

def seq_report_stats(seqs):
    '''
    print out various statistics on the sequences
    '''
    print "Read",len(seqs),"sequences"
    seqlen = [len(s.str) for s in seqs]
    print "Sequences have between",min(seqlen),"and",max(seqlen),"characters"
    maxlen = max(seqlen)

    ## Histogram of sequence lengths
    seqlend=Counter()
    for s in seqs:
        lenstr = len(s.str)
        seqlend[lenstr] += 1
    for lenstr in sorted(seqlend.keys()):
        print "len = %d, number of sequences: %d"%(lenstr,seqlend[lenstr])


def seq_filter_repeats(seqs,tmr):
    """ Return a list of SequenceSample's which avoids those
    sequences with too many repeats (more than tmr)
    """

    if tmr == 0:
        return seqs

    EPIMER=config.EPIMER
    sslist = []
    for s in seqs:
        episeq = episeq_from_string(s.str,rmdash=True)
        repeats_in_str = 0
        for epi in set(episeq):
            if episeq.count(epi) > 1:
                repeats_in_str += 1
        if repeats_in_str < tmr:
            sslist.append(s)
    if len(seqs) > len(sslist):
        print "Read",len(seqs),"sequences,",
        print "Removed",len(seqs)-len(sslist),"sequences,",
        print "Remaining",len(sslist),"sequences."
    return sslist            

def episeq_from_string(str,rmx=True,rmdash=True):
    """
    input a string and get a list of epitopes of length config.EPIMER
    rmdash: remove '-' from string before creating list of epitopes
    rmx: remove epitopes that contain 'x'
    """
    EPIMER=config.EPIMER
    if rmdash:
        re_dash=re.compile('-')
        str = re_dash.sub('',str)
    if rmx:
        episeq = [str[n:n+EPIMER] \
                   for n in range(len(str)-EPIMER+1)\
                   if 'x' not in str[n:n+EPIMER]]
    else:
        episeq = [str[n:n+EPIMER] \
                   for n in range(len(str)-EPIMER+1)]

    return episeq


def episeq_to_string(episeq):
    ## Convert a consistent list [ABC, -BC, BCD, CDE ...] to string "A-BC..."
    s = ""

    for epi in episeq:
        s += epi[0]     ## first character of each epitope
    s += episeq[-1][1:] ## all but first character [1:] of last epitope [-1]
    return s

def check_duplicated_epi(s):
    dups=set()
    pcount = Counter()
    p = episeq_from_string(s,rmdash=True)
    for epi in p:
        pcount[epi] += 1
        if pcount[epi]>1:
            print "Duplicated epitope: ",epi
            dups.add(epi)
    return dups


def epitope_tally(seqs,rmx=True,rmdash=True):
    ''' Input is a list of SequenceSamples
    Returns a Counter object, epicount, 
    for which epicount[epitope] is the number of 
    sequences in which the given epitope appears
    '''
    epicount = Counter()
    for s in seqs:
        episeq = episeq_from_string(s.str,rmx=rmx,rmdash=rmdash)
        for epi in set(episeq):
            epicount[epi] += 1
    return epicount

def remove_bad_epitopes(epicount):
    for epi in list(epicount.keys()):
        if 'x' in epi:
            epicount.pop(epi)

def truncate_strings(seqs,lohi):
    lo,hi = lohi
    if lo>hi:
        raise RuntimeError, "ERROR: lo>hi"
    for seq in seqs:
        if lo > len(seq.str):
            raise RuntimeError, "ERROR! lo too high: %d > %d" % (lo,len(seq.str))
        seq.str = seq.str[lo:hi]
    return seqs



def make_posn_adjust_map(str):
    """ returns an array posn_adjust[] such that
    long_offset = posn_adjust[short_offset]
    where long_offset is position in the original string (str)
    and short_offset is position in dash-removed string
    """
    posn_adjust = []
    for n_long in range(len(str)):
        if str[n_long] != "-":
            posn_adjust.append(n_long)
    return posn_adjust

def make_epi_posn_dict(long_str):
    re_dash=re.compile('-')
    short_str = re_dash.sub('',long_str)
    posn_adjust = make_posn_adjust_map(long_str)

    ## dictionary of positions of epitopes in the long string
    ## (what if same epitope appears in two places??)
    posn_in_long_string_of_epi={}

    for off_short in range(len(short_str)-EPIMER+1):
        off_long = posn_adjust[off_short]
        end_long = posn_adjust[off_short+EPIMER-1]
        epi = long_str[off_long:end_long+1]
        if (re_x.find(epi)):
            continue
        posn_in_long_string_of_epi[epi] = off_long

    return posn_in_long_string_of_epi
    


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("filename",help="input tbl filename");
    args = parser.parse_args()
    
    seqs = read_tbl(args.filename)

    print "Read",len(seqs),"sequences"
    print "Name of first sequence:",seqs[0].name
    print "Name of last sequence: ",seqs[-1].name
    print "Sequence length:",len(seqs[0].str)
    print "Sequence length:",len(seqs[1].str)
    print "Sequence length:",len(seqs[2].str)
    print "Sequence length:",len(seqs[-1].str)
    print "Last sequence:",seqs[-1].str

    
    posn_adj = make_posn_adjust_map(seqs[0].str)
    print posn_adj

    import pylab
    pylab.plot(posn_adj)
    pylab.show()
