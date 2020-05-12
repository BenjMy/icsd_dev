# -*- coding: utf-8 -*-
"""
Created on Mon May 11 16:22:08 2020
@author: Benjamin
Estimation of initial model based on the physical assumption that a single source current
can describe the pattern of the masse anomaly
"""
from numpy import linalg as LA
import numpy as np
from scipy.stats import pearsonr


#%%


### Individual Misfit  
def normF1(A,b):
    """compute the norm between observation data and individual green functions"""
    F1=[]
    for i in range(np.shape(A)[1]):
        F1i = LA.norm((b-A[:,i]))
        F1.append(F1i)
    norm_F1 = (F1 - min(F1)) / (max(F1) - min(F1)) 
    # normalise such as the sum is equal to 1 for current conservation
    
    return norm_F1

def misfit_2_initialX0(A,b,iniM0):
    
    norm_F1 = normF1(A,b)
    
    """Transform the misfit to an initial solution M0"""
    x0F1=1./((norm_F1+1)*(norm_F1+1)) # Inverse misfit using a 1/x^2 transformation
    x0F1_sum= x0F1/sum(x0F1) # normalize such as sum equal to 1
    if iniM0=='cst':
        x0F1_sum=np.ones(x0F1_sum.shape)*0.1
    M0=x0F1_sum
    
    return norm_F1, M0


def productmoment(A,b):
    """ Compute the product moment correlation after Binley et al. 1999
    .. math:: 

        r_{k}= \frac{\sum_{i}(D_{I}-\overline{D})(F_{i}(I_{k})-\overline{F}(I_{k}))}{\sqrt{\sum_{i}(D_{I}-\overline{D})^{2}}\sum_{i}(F_{i}(I_{k})-\overline{F}(I_{k}))^{2}}
    where $D_{i}$ is the $i^{th}$ measured transfer resistance and $F_{i}(I_{k})$ is the $i^{th}$  transfer resistance computed to unit current at location k. 
    """
    # self.icsd_init()            
    # Estimate a covariance matrix, given data observation and weights and tranfert resitances measured.
    rpearson=[]
    for i in range(np.shape(A)[1]):
        corr, _ = pearsonr(b, A[:,i])
        rpearson.append(corr)
    M0=rpearson
    
    return M0
        
        
