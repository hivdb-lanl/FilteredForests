#!/usr/local/bin/python2.7

""" Make training file from 'XXXX-Train.csv' and 'XX-sig.txt'
"""

from __future__ import division,print_function

import re
import sys
import math
import argparse

import util

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    paa = parser.add_argument
    paa("--input",help="Input csv data file")
    paa("--key",help="Input signature key table")
    paa("--abname",help="Sspecify name of antibody")
    paa("--contactonly","-K",action="store_true",help="Only include contact=1 features")
    paa("--sumonly","-S",action="store_true",help="Only include count of +1s and count of -1s")
    paa("-v","--verbose",action="count",help="verbose")
    args = parser.parse_args()

    def vprint(*pargs,**kwargs):
        if args.verbose:
            print(*pargs,file=sys.stderr,**kwargs)

    rows,header = util.read_table(args.input)
    sigkey = util.read_sigtable(args.key)

    sig_column = -1
    for i,h in enumerate(header):
        if re.search(".*signatures*",h):
            sig_column = i
        
    if sig_column == -1:
        raise RuntimeError("signatures column not identified")

    vprint("Signature column: %s"%(sig_column,))

    if args.abname:
        if args.abname not in header:
            raise RuntimeError("Antibody [%s] not in header: %s"%(args.abname,repr(header)))
        else:
            ab_column = [i for i in range(len(header)) if header[i]==args.abname]
            ab_column = ab_column[0]


        vprint("Antibody: Column[%d]=%s"%(ab_column,args.abname))
    else:
        args.abname=""

    #KeyNames = " ".join(["Key%02d"%(i,) for i in range(len(sigkey.posn))])))
    KeyNames = ""
    if args.sumonly:
        KeyNames = "Count+1 Count-1"
    else:
        for i,n in enumerate(sigkey.posn):
            if args.contactonly==True and sigkey.contact[i]==0:
                continue ## don't keep
            KeyNames += " Key%02d"%(i,)

    if args.contactonly:
        print("%s %s %s" %
              (header[0],args.abname,
               KeyNames))
    else:
        print("%s %s %s %s" %
              (header[0],args.abname,
               " ".join(header[1:sig_column]),KeyNames))


    for r in rows:
        name = r[0]
        in_dat = []
        ab_out = r[ab_column] if args.abname else ""
        if re.search("NA",ab_out): # == "NA":
            continue

        ## Convert to log10
        if ab_out == "100":
            ab_out = "2"
        else:
            ab_out = str( math.log10(float(ab_out)) )

        if not args.contactonly:
            for i in range(1,sig_column):
                in_dat.append(r[i])

        seq = r[sig_column]
        vals = util.seq_to_sigkey(sigkey,seq,contact=args.contactonly,sum=args.sumonly)
        in_dat.extend(vals)
        print("%s %s %s" % (name,ab_out," ".join(map(str,in_dat))))


        


            


        
