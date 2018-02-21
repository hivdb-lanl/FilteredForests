
import numpy as np
import scipy.stats as stats

from sklearn import linear_model,ensemble,svm,model_selection
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn import metrics

def get_model(rc="r",model="LIN"):
    
    if rc == "r":
        ## REGRESSORS
        if model == "LIN": 
            rgrsr = linear_model.LinearRegression()
        elif model == "LASSO":
            rgrsr = linear_model.LassoLars(alpha=0.0)
        elif model == "RF":
            rgrsr = ensemble.RandomForestRegressor(n_estimators=120,
                                                   bootstrap=True,
                                                   max_depth=None,
                                                   oob_score=True)
        elif model == "RFX":
            rgrsr = ensemble.ExtraTreesRegressor(n_estimators=120,
                                                 bootstrap=True,
                                                 max_depth=None,
                                                 oob_score=True,
                                                 max_features=None)
        elif model == "RFXX": 
            rgrsr = ensemble.ExtraTreesRegressor(n_estimators=12,
                                                 bootstrap=True,
                                                 max_depth=None,
                                                 max_features=None)
        elif model == "SV":  
            rgrsr = svm.SVR(epsilon=0.0,tol=0.1)
        else:
            raise RuntimeError("Invalid model [%s]"%(model,))

    if rc == "c":
        ## CLASSIFIERS
        if model == "LIN":
            rgrsr =  LinearDiscriminantAnalysis()
        elif model == "RF":
            rgrsr =  ensemble.RandomForestClassifier(n_estimators=120,
                                                     bootstrap=True,
                                                     max_depth=None,
                                                     oob_score=True)
        elif model == "RFX":
            rgrsr =  ensemble.ExtraTreesClassifier(n_estimators=120,
                                                   bootstrap=True,
                                                   max_depth=None,
                                                   oob_score=True,
                                                   max_features=None)
        elif model == "RFXX":
            rgrsr =  ensemble.ExtraTreesClassifier(n_estimators=12,
                                                   bootstrap=True,
                                                   max_depth=None,
                                                   max_features=None)
        elif model == "SV":
            rgrsr =  svm.SVC()
        else:
            raise RuntimeError("Invalid model [%s]"%(model,))

    return rgrsr


def train_model(y,x,rc="r",model="LIN"):
    rgrsr = get_model(rc=rc,model=model)
    if rc == "r":
        rgrsr.fit(x,y)
    elif rc == "c":
        #yy = [int(yc!=TOPVAL) for yc in y]
        rgrsr.fit(x,y)
    else:
        raise RuntimeError("rc=%s should be r or c"%(rc,))

    return rgrsr


def cv_model(nsplits):
    if nsplits <= 0:
        cv = None
    elif nsplits==1:
        cv = model_selection.LeaveOneOut()
    else:
        cv = model_selection.KFold(n_splits=nsplits) 
    return cv


TOPVAL=2 ## top value (corresponds to "negative" classification)

class FitScore:
    pass

def evalclsfr_mcc(ytru,yest):
    mcc = metrics.matthews_corrcoef(ytru,yest)
    return mcc

def evalclsfr_contigency(ytru,yest):
    M = metrics.confusion_matrix(ytru_bucket,yest_bucket)
    _,pval = stats.fisher_exact(M)
    [[TN, FP],[FN, TP]] = M
    print("[%d %d][%d %d]: %f"%(TN,FP,FN,TP,pval))

def evalclsfr_summarize(ytru,yest,pre="",quiet=False):
    f = FitScore()
    if len(ytru) == 0:
        f.mcc = f.accy = f.tpos = 0
        f.pval = 1
        return f

    M = metrics.confusion_matrix(ytru,yest)
    #print("M=",M)
    _,pval = stats.fisher_exact(M)
    [[TN, FP],[FN, TP]] = M
    if TN==0 and FN==0 or (FP==0 and TP==0):
        mcc=0
    else:
        mcc = metrics.matthews_corrcoef(ytru,yest)
    if not quiet:
        print("%s: mcc= %8.4f [TN FP][FN TP]=[ %3d %3d ][ %3d %3d ]: p= %8.2e"%(pre,mcc,TN,FP,FN,TP,pval))
    
    f.mcc = mcc
    f.pval = pval
    f.accy = (TP+TN)/(TP+TN+FP+FN)
    f.tpos = max([TP+FN,TN+FP])/(TP+TN+FP+FN)
    return f


def evalrgrsr_classic_r2(ytru,yest):
    _,_,r,_,_ = stats.linregress(ytru,yest)
    return r*r

def evalrgrsr_r2(ytru,yest):
    rr = metrics.r2_score(ytru,yest)
    if rr < -1.0e6:
        rr = -9.9999
    return rr

def evalrgrsr_kendalltau(ytru,yest):
    return stats.kendalltau(ytru,yest)

def evalrgrsr_abserr(ytru,yest):
    ''' average abslute error '''
    mae = np.sum(np.abs(np.array(yest)-np.array(ytru)))/len(ytru)
    if mae > 1.0e6:
        mae = 99.9999
    return mae

def evalrgrsr_summarize(ytru,yest,pre="",quiet=False):
    f = FitScore()
    if len(ytru)==0:
        f.crr = f.ae = 0
        f.pval = 1
        return f

    ytru = np.array(ytru)
    yest = np.array(yest)

    ndx_tn = [i for i,y in enumerate(ytru) if y<TOPVAL]

    rr = evalrgrsr_r2(ytru,yest)
    rrtn = evalrgrsr_r2(ytru[ndx_tn],yest[ndx_tn])

    crr = evalrgrsr_classic_r2(ytru,yest)
    crrtn = evalrgrsr_classic_r2(ytru[ndx_tn],yest[ndx_tn])

    ae = evalrgrsr_abserr(ytru,yest)
    aetn = evalrgrsr_abserr(ytru[ndx_tn],yest[ndx_tn])

    rho,pval = evalrgrsr_kendalltau(ytru,yest)
    rho_tn,pval_tn = evalrgrsr_kendalltau(ytru[ndx_tn],yest[ndx_tn])

    #print(pre,rr,crr,ae,rrtn,crrtn,aetn,rho,pval)
    #print("%s: r2= %7.4f cr2 = %7.4f ae= %7.4f p=%8.2e [w/o: r2 %9.3f cr2: %6.3f ae %7.4f kt %7.4f p %8.2e]" 
    #      % (pre,rr,crr,ae,pval,rrtn,crrtn,aetn,rho_tn,pval_tn))
    if not quiet:
        print("%s: mae= %7.4f cr2= %7.4f p= %8.2e"%(pre,ae,crr,pval))

    f.crr = crr
    f.ae = ae
    f.pval = pval
    return f
