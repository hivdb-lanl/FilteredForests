#!/usr/bin/env python

from __future__ import division, print_function

import numpy as np
import re

def read_table(filename):

    with open(filename,"rU") as f:
        header = f.readline().strip().split(',')
        nheader = len(header)
        rows=[]
        for line in f:
            txt = line.strip().split(',')
            if len(txt) != nheader:
                raise RuntimeError("mismatched line: [%s] should have %d items" %(line,nheader))
            rows.append(txt)

    return rows,header


class SigKey:
    pass

def read_sigtable(filename):
    ''' similar to read_keytable in keytable.py '''
    
    sigkey = SigKey()

    sigkey.sen = []
    sigkey.res = []
    sigkey.contact = []
    sigkey.hxb2 = []
    sigkey.posn = []
    with open(filename,"rU") as f:
        header = f.readline().strip().split()
        nprev=0
        for line in f:
            (n,h,_,_,s,r,c) = line.split()
            n = int(n)
            if n <= nprev:
                raise RuntimeError("Indices out of order")
            nprev = n

            sigkey.posn.append(int(n))
            sigkey.sen.append(s)
            sigkey.res.append(r)
            sigkey.contact.append(int(c))
            sigkey.hxb2.append(int(h))

    return sigkey

def seq_to_sigkey(key,seq,contact=False,sum=False):
    '''convert a sequence into a list of signature values: +1, 0, -1 
       if sum==True, then return only two values, count of +1s and count of -1s '''

    if len(seq) < max(key.posn):
        raise RuntimeError("seq %s is too short" % (seq,))

    if sum:
        valkeep = [0,0]
    else:
        valkeep = []
    for i,n in enumerate(key.posn):
        if contact==True and key.contact[i]==0:
            continue ## don't keep
        val=0
        if seq[n-1] in key.sen[i]:
            val = 1
        if seq[n-1] in key.res[i]:
            val = -1
        if sum:
            if val ==  1: valkeep[0] += 1
            if val == -1: valkeep[1] += 1
        else:
            valkeep.append(val)

    return valkeep

class snxyClass:
    pass

def snxy_fromfile(filename,sans_y=False,sigonly=False):
    '''
    in addition to x,y from a file,
    get sample names and feature names
    return structure snxy that has components:
      snames: list of N sample names
      xfnames: list of F feature names
      yfname: name for the y feature
      x: 2D array of N rows and F columns of feature values
      y: 1D array of N rows of output values
    '''

    with open(filename,"rU") as f:
        lines = f.readlines()

    header,lines = lines[0],lines[1:]
    fnames = header.split()
    fnames = fnames[1:] ## remove VirusID,first is output Ab name, rest are feature names

    ## Produce numpy arrays x,y
    nfeatures = 0
    badlines=[]
    for ndx,line in enumerate(lines):
        if sans_y:
            try:
                (name,x) = line.split(None,1)
            except ValueError:
                raise RuntimeError("Bad line: %s"%(line,))
           #(name,x) = line.split(None,1)
        else:
            (name,y,x) = line.split(None,2)
            if y=="NA":
                badlines.append(ndx)
        if nfeatures==0:
            nfeatures = len(x.split())
        elif nfeatures != len(x.split()):
            raise RuntimeWarning("Inconsistent number of features [%s]" % (name,))

    tmplines = []
    for ndx,line in enumerate(lines):
        if ndx not in badlines:
            tmplines.append(line)
    lines = tmplines


    nsamples = len(lines)
    x = np.zeros((nsamples,nfeatures))
    if not sans_y:
        y = np.zeros((nsamples,))
    snames = []
    for n,line in enumerate(lines):
        if not sans_y:
            (name,yval,xval) = line.split(None,2)
            y[n] = yval
        else:
            try:
                (name,xval) = line.split(None,1)
            except ValueError:
                raise RuntimeError("Bad line: %s"%(line,))

        snames.append(name)
        try:
            xval = re.sub("NA","NaN",xval)
            x[n,:] = [float(xi) for xi in xval.split()]
        except ValueError:
            raise RuntimeError("bad line: %s"%(xval,))

    if sans_y:
        yfname,xfnames = None,fnames
    else:
        yfname,xfnames = fnames[0],fnames[1:]
    
    if len(xfnames)>0 and len(xfnames) != len(x[0,:]):
        print("fnames:",fnames)
        raise RuntimeError("Mismatched count of feature names")

    snxy = snxyClass()
    
    snxy.snames = snames
    snxy.xfnames = xfnames
    snxy.yfname = yfname
    snxy.y = y if not sans_y else None
    snxy.x = x

    return snxy



        
    






        
                                   
                    
            
        
