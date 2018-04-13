# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 14:37:31 2018

@author: nhermans
"""
import Tools
import numpy as np

#functions
def fjc(f, par = Tools.default_pars()): #WLC extensible + WLC
    return 1-0.5*np.sqrt(par['kBT_pN_nm']/(f*par['Pss_nm']))

def wlc(f, par = Tools.default_pars()):
    """returns the WLC as a fraction of the DNA contour length"""
    return (1-0.5*np.sqrt(par['kBT_pN_nm']/(f*par['P_nm']))) + f/par['S_pN']

def wlc_fjc(f, SeqLen = 3000, par = Tools.default_pars()):
    return par['CLds_bp']*par['DNAds_nm']*wlc(f) + SeqLen*2*par['DNAss_nm']*fjc(f)

def fit_p(f, p, SeqLen = 3000, par = Tools.default_pars()): #WLC extensible + WLC
    return par['CLds_bp']*par['DNAds_nm']*wlc(f) + SeqLen*2*par['DNAss_nm']*(1-0.5*np.sqrt(par['kBT_pN_nm']/(f*p)))

def GC_dG(sequence, window):
    """ #dG calculator based on the nearest neighbour energies 
    Window: Specify the window to calculate the dG
    sequence: sequence string, upper case, only the sequence in the direction of unzipping"""
    
    #nearest neighbor parameters, value: kcal/mol
    NNparam =  {'AA': 1.0, 'TT': 1.0, 'AT': 0.88, 'TA': 0.58, 'CA': 1.45, 'TG': 1.45, 'GT': 1.44, 'AC': 1.44, 'CT': 1.28,
    'AG': 1.28, 'GA': 1.30, 'TC': 1.30, 'CG': 2.17, 'GC': 2.24, 'GG': 1.84, 'CC': 1.84}
    
    vals = []
    values = []
    for i in range(0, len(sequence)-2):
        s = sequence[i: i + 2] #dinucleotides
        dG = NNparam[s]
        vals.append(dG)
    for i in range(0,len(vals)): #lowpass filter
        values.append(np.average(vals[i:i+window]))
    return values