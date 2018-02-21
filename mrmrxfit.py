#!/usr/local/bin/python2.7

from __future__ import division,print_function

import re
import sys
import time
import math
import numpy as np
import scipy.stats as stats

import argparse

#import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import Counter

from sklearn import utils

from sequtil import read_tbl

import util
import mlutil

from SequenceSample import SequenceSample

from mrmr import mrmr
from hxb2 import hxb2


def get_hxb2_sequence(seqf,hxb2name="HXB2"):
    hxb2_seq=None
    for s in seqf:
        if s.name == hxb2name:
            hxb2_seq = s.str
            break
    if not hxb2_seq:
        raise RuntimeError("Failed to find HXB2 sequence in file [%s]" % (args.train,))
    return hxb2_seq

def select_sequences_fromfile(filename,seq_select_names=[]):
    seq_avail = read_tbl(filename)
    seq_avail_names = [s.name for s in seq_avail]
    seq_select = []

    ## if not specified, then return all sequences
    if len(seq_select_names)==0:
        seq_select_names = seq_avail_names

    for sname in seq_select_names:
        if sname not in seq_avail_names:
            raise RuntimeError("Sequence name [%s] not available in file [%s]" % (sname,filename))
        ndx = seq_avail_names.index(sname)
        seq_select.append(seq_avail[ndx])
    return seq_select

def make_features(seqf,hx,minc=0):
    feature_set = make_feature_set(seqf,hx)
    features = make_feature_dict(seqf,feature_set,hx)
    if minc:
        remove_lopsided_features(features,minc)
    return features


def make_feature_set(seqf,hx):
    feature_set = set()
    for s in seqf:
        if s.name == "HXB2":
            continue
        for n,c in enumerate(s.str):
            fname = hx.fname(n,c)
            feature_set.add(fname)
    return feature_set

def make_feature_dict(seqf,feature_set,hx):
    feature = dict()
    for f in feature_set:
        feature[f] = [False]*len(seqf)

    for ns,s in enumerate(seqf):
        for n,c in enumerate(s.str):
            fname = hx.fname(n,c)
            if fname in feature_set:
                feature[fname][ns] = True

    return feature

def remove_lopsided_features(features,minc):
    bad_features=[]
    for f in features:
        for tf in [True,False]:
            if sum( [feature[f][n]==tf for n in range(len(seqf)) ] ) < minc:
                bad_features.append(f)
    for f in bad_features:
        feature.pop(f)
    return bad_features


def write_features(filename,f_list,features,seq_names):
    with open(filename,"w") as file:
        file.write("%20s"%"VirusID")
        for f in f_list:
            file.write("\t%6s" % (f,) )
        file.write("\n")
        for n,sname in enumerate(seq_names):
            file.write("%20s" % (sname,))
            for f in f_list:
                p = 1 if feature[f][n] else -1
                file.write("\t%5d" % (p,))
            file.write("\n")

## Read yfile, get output values and associated sample names
def read_yfile(yfile,ycolumn):

    snxy = util.snxy_fromfile(yfile,sans_y=True)
    col_ndx = snxy.xfnames.index(ycolumn)
    y = snxy.x[:,col_ndx]
    ok_ndx = np.isfinite(y)
    y = y[ok_ndx]
    snames = np.array(snxy.snames)[ok_ndx]
    ## binary version of y
    ybool = np.array(y >= 100, dtype=np.bool) ## Note hardcoded 100
    y = np.log10(y)                           ## Log_10

    return snames,y,ybool

def make_features_numeric(features,flist):
    #nseq = len(features[flist[0]])
    nseq = next(len(features[f]) for f in features)
    xf = np.empty((nseq,len(flist)),dtype=np.float)
    for nf,f in enumerate(flist):
        if f in features:
            xf[:,nf] = features[f]
        else:
            xf[:,nf] = [False] * nseq
    return xf


if __name__ == "__main__":

    ## read command line options
    parser = argparse.ArgumentParser()
    paa = parser.add_argument
    paa("--train","-T",required=True,help="Input training data file")
    paa("--yfile",required=True,help="File with IC50 for various ab's")
    paa("--ab",required=True,help="Name of antibody of interest")
    paa("--rc",default="r",
        help="rc=c for classification, rc=r for regression")
    paa("--holdout",help="Holdout testing file")
    paa("--yholdout",help="Holdout IC50 values")
    paa("--minc",type=int,default=0,help="mininum count of +1s and -1s; else delete feature")
    paa("-n","--nfeatures",type=int,default=5,help="number of features")
    paa("--onlynfeatures",action="store_true",help="Only n features, not 1,2,...,n")
    paa("--model",default="RFX",help="Choose: LIN,SVR,LASSO,RF,RFX")
    paa("--cv",type=int,default=0,help="Cross validation")
    paa("-F",help="filename for writing out features")
    paa("--savefig",help="Save plots to file named <savefig>-xxx.eps")
    paa("--noshow",action="store_true",help="Do not show plots")
    paa("-v","--verbose",action="count",help="verbose")
    args = parser.parse_args()

    def vprint(*pargs,**kwargs):
        if args.verbose:
            print(*pargs,file=sys.stderr,**kwargs)

    def savefigfile(s):
        return args.savefig + "-" + s if args.savefig else None

    def show(ytru,yest,title=None,save=None):
        if args.noshow and not args.savefig:
            return
        if args.rc == "c":
            return
        plt.figure()
        plt.plot(ytru,yest,".")
        plt.plot([-3,2],[-3,2],'g--')
        plt.xlim([-3.2,2.1])
        plt.xlabel("Measured log10 IC50")
        plt.ylabel("Estimated log10 IC50")
        if title:
            plt.title(title)
        if save:
            plt.savefig(save + ".eps")

    ## Read sequence file, get HXB2 sequence
    hxb2_seq = select_sequences_fromfile(args.train,["HXB2"])[0]
    hx = hxb2(hxb2_seq.str)

    ## Read yfile, get y values and seq_names
    seq_names,yval,ybin = read_yfile(args.yfile,args.ab)

    seq_names,yval,ybin = utils.shuffle(seq_names,yval,ybin)

    ## Read sequence file, get sequences assocatied w/ seq_names
    seqf = select_sequences_fromfile(args.train,seq_names)

    ## Make features
    to = time.clock()
    feature = make_features(seqf,hx,minc=args.minc)
    vprint("make %d features: %f sec" % (len(feature),time.clock()-to,))

    ### Ok, now lets go compute some mutual informations
    ### So we can select best mRMR features
    to = time.clock()
    M = mrmr(ybin,feature)
    vprint("initialize mrmr: %f sec" % (time.clock()-to,))

    to = time.clock()
    flist = M.get_features(args.nfeatures)
    fcorr = [M.corr(f) for f in flist]
    vprint("get_features: %f sec" % (time.clock() - to,))
    #print("Features:",flist[:10])
    #print("Correlation:",fcorr[:10])


    if args.F:
        write_features(args.F,flist,feature,[s.name for s in seqf])
    else:
        for n,f in enumerate(flist):
            print("Ab: %15s  feature[%03d]: %-6s  mu: %f" % 
                  (args.ab,n+1,f,M.mutinfo_y(f)))

    ## Convert features into a numeric array
    xf = make_features_numeric(feature,flist)


    if args.holdout:
        ## Read holdout data, make h_y and h_x
        h_seq_names,h_yval,h_ybin = read_yfile(args.yholdout,args.ab)
        h_seqf = select_sequences_fromfile(args.holdout,h_seq_names)
        h_feature = make_features(h_seqf,hx)
        h_xf = make_features_numeric(h_feature,flist)

    if args.cv:
        cv = mlutil.cv_model(args.cv)

    evalfcn = mlutil.evalrgrsr_summarize
    if args.rc=="c":
        evalfcn = mlutil.evalclsfr_summarize
        yval   = np.array([int(y) for y in   ybin])
        if args.holdout:
            h_yval = np.array([int(y) for y in h_ybin])

    for nf in range(args.nfeatures):
        if args.onlynfeatures and nf+1 < args.nfeatures:
            continue

        rgrsr = mlutil.train_model(yval,xf[:,:nf+1],rc=args.rc,model=args.model)
        yest = rgrsr.predict(xf[:,:nf+1])

        f_is = evalfcn(yval,yest,pre="%02d is"%(nf+1,))

        if nf+1 == args.nfeatures and args.rc=="r":
            show(yval,yest,title="In-sample",save=savefigfile("ins"))

        if args.holdout:
            h_yest = rgrsr.predict(h_xf[:,:nf+1])
            f_ho = evalfcn(h_yval,h_yest,pre="%02d ho"%(nf+1,))
            if nf+1 == args.nfeatures:
                show(h_yval,h_yest,title="Holdout",save=savefigfile("ho"))

        if args.cv:

            ## This version of cross-validation is cheating ... because we are using the 
            ## full dataset to choose the features.  We go ahead and do this here as a kind of
            ## baseline/benchmark; later in the code, we'll do it correctly.

            yest_bucket=[]
            ytru_bucket=[]

            for train,test in cv.split(yval):
                cv_rgrsr = mlutil.train_model(yval[train],xf[train,:nf+1], rc=args.rc, model=args.model)
                yest = cv_rgrsr.predict(xf[:,:nf+1])
                for n in test:
                    yest_bucket.append(yest[n])
                    ytru_bucket.append(yval[n])

            if nf+1 == args.nfeatures:
                show(ytru_bucket,yest_bucket,title="BadX-Val",
                     save=savefigfile("badxval"))

            f_xv = evalfcn(ytru_bucket,yest_bucket,pre="%02d xv"%(nf+1,))

        if nf+1 == args.nfeatures:
            ## identify the best features
            rgrsr.oob_score_
            ndx = np.argsort(rgrsr.feature_importances_)
            ntop = len(flist) ## min([10,len(flist)])
            bestfeaturelist = [flist[ndx[-k]] for k in range(1,1+ntop)]
            corrfeaturelist = np.sign([fcorr[ndx[-k]] for k in range(1,1+ntop)])
            ## we'll wait till the end to print out the best features

    ## OK, Now we are going to do CV correctly!
    if args.cv:
        train_list = []
        test_list = []
        cv_xf_list = []
        for train,test in cv.split(yval):
            train_list.append(train)
            test_list.append(test)

            cv_feature = dict()
            for f in feature:
                cv_feature[f] = [feature[f][n] for n in range(len(seqf)) if n in train]
            cv_M = mrmr(ybin[train],cv_feature)
            cv_flist = cv_M.get_features(nf+1)
            cv_xf = make_features_numeric(feature,cv_flist)
            cv_xf_list.append( cv_xf )

        for nf in range(args.nfeatures):
            if args.onlynfeatures and nf+1 < args.nfeatures:
                continue

            yest_bucket=[]
            ytru_bucket=[]

            for ncv,(train,test) in enumerate(zip(train_list,test_list)):
                cv_xf = cv_xf_list[ncv]
                cv_rgrsr = mlutil.train_model(yval[train],cv_xf[train,:nf+1], rc=args.rc, model=args.model)
                yest = cv_rgrsr.predict(cv_xf[:,:nf+1])
                for n in test:
                    yest_bucket.append(yest[n])
                    ytru_bucket.append(yval[n])
        
            if nf+1 == args.nfeatures:
                show(ytru_bucket,yest_bucket,title="Cross-Val %02d"%(nf+1,),
                     save=savefigfile("crossval"))

            f_cv = evalfcn(ytru_bucket,yest_bucket,pre="%02d cv"%(nf+1,))


    ## Print out Best Features at the very end
    ## corrfeaturelist and bestfeaturelist were computed above, on the last iteration
    Bestfeatures = []
    for c,b in zip(corrfeaturelist,bestfeaturelist):
        if c < 0:
            Bestfeatures.append( b.upper() )
        else:
            Bestfeatures.append( b.lower() )

    ntop = min([10,len(flist)]) ## HARDCODED '10' -- no more than 10 top features
    print("%15s %2d/%2d %s"%(args.ab,ntop,len(flist)," ".join(["%5d"%(corrfeaturelist[k],) for k in range(ntop)])))
    print("%15s %2d/%2d %s"%(args.ab,ntop,len(flist)," ".join(Bestfeatures[:ntop])))


    if not args.noshow:
        plt.show()



    


