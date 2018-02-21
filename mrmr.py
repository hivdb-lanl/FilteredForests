from __future__ import print_function, division
from collections import Counter
import numpy as np

def entropy(z):
    zz = z[z>0]
    return np.sum( -zz*np.log(zz) )

def mutual_info(ta,tb):
    z = Counter( zip(ta,tb) )
    z = np.array([[z[False,False], z[False,True]],
                  [z[True,False], z[True,True]]],dtype=np.float)
    assert( np.sum(z) > 0 )
    z = z/np.sum(z)
    za = np.sum(z,axis=0)
    zb = np.sum(z,axis=1)
    mu_info = entropy(za)+entropy(zb)-entropy(z)
    corr = (z[0,0]+z[1,1]-z[1,0]-z[0,1])/np.sum(z)
    return mu_info,corr

class mrmr:
    
    def __init__(self,y,xf):
        ''' y: list of label values
            xf: dictionary of lists; key is feature name, list is feature values
            values for y and xf lists are binary, each element is either True or False
        '''
        self.y = y
        self.xf = xf
        self.mu = dict()
        self.mu_y = dict()
        self.corr_y = dict()

        for f in self.xf:
            self.mu_y[f],self.corr_y[f] = mutual_info(self.y,self.xf[f])

        self.flist = sorted(self.mu_y, key=self.mu_y.get, reverse=True)
    
    def mutinfo(self,fa,fb):
        if fa+fb not in self.mu:
            self.mu[fa+fb],_ = mutual_info(self.xf[fa],self.xf[fb])
            self.mu[fb+fa] = self.mu[fa+fb]
        return self.mu[fa+fb]

    def mutinfo_y(self,f):
        if f not in self.mu_y:
            self.mu_y[f],self.corr_y[f] = mutual_info(yy,self.xf[f])
        return self.mu_y[f]

    def corr(self,f):
        if f not in self.corr_y:
            self.mu_y[f],self.corr_y[f] = mutual_info(yy,self.xf[f])
        return self.corr_y[f]

    def alt_next_feature(self,f_sofar):
        ''' slow but more direct version '''
        muscore = dict()
        for f in set(self.xf.keys()) - set(f_sofar):
            muscore[f] = self.mutinfo_y(f)
            if f_sofar:
                mu_redundant = sum( [self.mutinfo(f,g) for g in f_sofar] ) / len(f_sofar)
                # Q: I wonder if mu_redundant = max( [...] ) would be a better chioce
                #mu_redundant = max( [self.mutinfo(f,g) for g in f_sofar] ) 
                muscore[f] -= mu_redundant            
    
        fbest = max(muscore, key=muscore.get)
        mubest = muscore[fbest]
        return fbest,mubest

    def next_feature(self,f_sofar):
        fbest = None
        mubest = -np.inf
        for f in self.flist:
            if f in f_sofar:
                continue
            mu = self.mutinfo_y(f)
            if mu < mubest:
                break
            if f_sofar:
                mu_redundant = sum( [self.mutinfo(f,g) for g in f_sofar] ) / len(f_sofar)
                mu -= mu_redundant
            if mu > mubest:
                fbest = f
                mubest = mu

        return fbest,mubest

    def get_features(self,n):
        f_sofar = []
        while len(f_sofar) < n:
            f,_ = self.next_feature(f_sofar)
            f_sofar.append(f)
        return f_sofar



    
    
