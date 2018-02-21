#!/usr/local/bin/python2.7
""" Make training file from 'XXXX-Train.csv' and 'XX-sig.txt'
"""

from __future__ import division,print_function

import re
import sys
import math
import numpy as np
import scipy.stats as stats

import argparse

#import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import util
import mlutil

from sklearn import metrics,model_selection

import warnings
warnings.filterwarnings("ignore",r"Variables are collinear.")
warnings.filterwarnings("ignore",r"invalid value encountered in long_scalars")

TOPVAL=2 ## top value (corresponds to "negative" classification)

def train_model(y,x,rc="r",model="LIN"):
    rgrsr = mlutil.get_model(rc=rc,model=model)
    if rc == "r":
        rgrsr.fit(x,y)
    elif rc == "c":
        yy = [int(yc!=TOPVAL) for yc in y]
        rgrsr.fit(x,yy)
    else:
        raise RuntimeError("rc=%s should be r or c"%(rc,))

    return rgrsr

## evalXXXXX_summarize: wrap the versions in mlutil

def evalclsfr_summarize(ytru,yest,pre=""):
    ytru = [int(y != TOPVAL) for y in ytru]
    f = mlutil.evalclsfr_summarize(ytru,yest,quiet=True)
    return f

def evalrgrsr_summarize(ytru,yest,pre=""):
    return mlutil.evalrgrsr_summarize(ytru,yest,quiet=True)



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    paa = parser.add_argument
    paa("--train",required=True,help="Input training data file (reqd)")
    paa("--holdout",help="Holdout testing file")
    paa("--cv",type=int,default=0,help="Cross validation folds: 0=no, 1=LOO")
    paa("--predict",help="Apply model to make predictions for this data file")
    paa("--rc",default="r",
        help="rc=c for classification, rc=r for regression")
    paa("--model",default="RFX",help="Choose: LIN,SVR,LASSO,RF,RFX")
    paa("--noshow",action="store_true",help="do not make plots")
    paa("--savefig",help="Save plots to file named <savefig>-xxx.eps")
    paa("--prefix",default="//",
        help="string used to prefix output lines (format: protein/antibody/featuretype)")
    paa("--title",help="title for plots")
    paa("--top",action="store_true",help="List top 5 features")
    paa("--xmin",type=float,default=-3.0,help="lower bound of plot")
    paa("-v","--verbose",action="count",help="verbose")
    args = parser.parse_args()

    eval_summarize = evalrgrsr_summarize
    if args.rc == "c":
        args.noshow = True
        eval_summarize = evalclsfr_summarize

    def vprint(*pargs,**kwargs):
        if args.verbose:
            print(*pargs,file=sys.stderr,**kwargs)

    def prefix(s):
        return args.prefix + " " + s if args.prefix else s

    def savefigfile(s):
        return args.savefig + "-" + s if args.savefig else None

    def show(ytru,yest,title=None,save=None):
        if args.rc=="c" or (args.noshow and not args.savefig):
            return
        plt.figure()
        plt.subplot(111,aspect='equal')
        plt.plot(ytru,yest,'.')
        #plt.xlim(xmax=2.1)
        #plt.ylim(ymax=2.1)
        xmin=-3
        if args.xmin:
            xmin=args.xmin
        plt.xlim([xmin-0.1,2.1])
        plt.ylim([xmin-0.1,2.1])
        plt.plot([xmin,2],[xmin,2],'g--',linewidth=1.5)
        plt.plot(np.log10([50,50]),[xmin,2],'r-',linewidth=0.5)
        #plt.axis('equal')
        plt.xlabel("observed")
        plt.ylabel("estimated")
        if args.title:
            plt.title(args.title)
        elif title:
            plt.title(title)
        if save:
            plt.savefig(save + ".eps")

    snxy = util.snxy_fromfile(args.train)

    rgrsr = train_model(snxy.y,snxy.x,rc=args.rc,model=args.model)
    yest = rgrsr.predict(snxy.x)

    if args.top:
        rgrsr.oob_score_
        ndx = np.argsort(rgrsr.feature_importances_)
        print("Top 5 Features",ndx[-5:])
        for k in range(1,1+5):
            print( "%2d %s %8.4f" % (
            ndx[-k],snxy.xfnames[ndx[-k]],rgrsr.feature_importances_[ndx[-k]]))

    f_is = eval_summarize(snxy.y,yest)

    #show(snxy.y,yest,title="In-sample: %s"%(prefix(""),),save=savefigfile('insampl'))

    if args.cv:
        yest_bucket=[]
        ytru_bucket=[]

        if args.cv == 1:
            cv = model_selection.LeaveOneOut()
        else:
            cv = model_selection.KFold(n_splits=args.cv) 

        for train,test in cv.split(snxy.x):
            cv_rgrsr = train_model(snxy.y[train],snxy.x[train],rc=args.rc,model=args.model)
            yest = cv_rgrsr.predict(snxy.x)
            for n in test:
                yest_bucket.append(yest[n])
                ytru_bucket.append(snxy.y[n])

        f_cv = eval_summarize(ytru_bucket,yest_bucket)
        show(ytru_bucket,yest_bucket,title="Cross-validation: %s"%(prefix(""),),save=savefigfile('crosval'))
    else:
        f_cv = eval_summarize([],[])

    if args.holdout:
        h_snxy = util.snxy_fromfile(args.holdout)

        for trn,tst in zip(snxy.xfnames,h_snxy.xfnames):
            if trn != tst:
                raise RuntimeError("Mismatch: %s != %s"%(trn,tst))

        yest = rgrsr.predict(h_snxy.x)

        f_ho = eval_summarize(h_snxy.y,yest)
    else:
        f_ho = eval_summarize([],[])

    if args.predict:
        p_snxy = util.snxy_fromfile(args.predict,sans_y=True)
        yest = rgrsr.predict(p_snxy.x)
        for n,y in zip(p_snxy.snames,yest):
            print("%25s %8.3f"%(n,y))

    pfx = args.prefix.split("/")
    if args.rc == "r":
        print(" Prot Antibody        Features  | Mean Abs Error |   Classic R^2  |        p-value")
        print("                                   Ins Xval Hout    Ins Xval Hout       Ins    Xval    Hout")
        print("%5s %-15s %-9s | %4.2f %4.2f %4.2f | %4.2f %4.2f %4.2f | %7.1e %7.1e %7.1e" %
              (pfx[0],pfx[1],pfx[2],f_is.ae,f_cv.ae,f_ho.ae,f_is.crr,f_cv.crr,f_ho.crr,f_is.pval,f_cv.pval,f_ho.pval))
    if args.rc == "c":
        print(" Prot Antibody        Features  |    Accuracy    |  Default-accy  |        MCC        |        p-value")
        print("                                   Ins Xval Hout    Ins Xval Hout     Ins  Xval  Hout       Ins    Xval    Hout")
        print("%5s %-15s %-9s | %4.2f %4.2f %4.2f | %4.2f %4.2f %4.2f | %5.2f %5.2f %5.2f | %7.1e %7.1e %7.1e" %
              (pfx[0],pfx[1],pfx[2],f_is.accy,f_cv.accy,f_ho.accy,f_is.tpos,f_cv.tpos,f_ho.tpos,f_is.mcc,f_cv.mcc,f_ho.mcc,f_is.pval,f_cv.pval,f_ho.pval))


    if not args.noshow:
        plt.show()



            


        
