# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 13:30:25 2018

@author: nhermans
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 14:52:17 2018

@author: nhermans
"""
#from lmfit import Parameters

import numpy as np
from scipy import signal

def Define_Handles(Select=True, Pull=True, DelBreaks=True, MinForce=2, MinZ=0, MedFilt=False):
    """If analysis has to be done on only part of the data, these options can be used"""
    Handles = {}
    Handles['Select'] = Select
    Handles['Pulling'] = Pull
    Handles['DelBreaks'] = DelBreaks
    Handles['MinForce'] = MinForce
    Handles['MinZ'] = MinZ
    Handles['MedFilt'] = MedFilt 
    return Handles

def default_pars():
    """Default fitting parameters, returns a {dict} with 'key'= paramvalue
    par['CLds_bp'] = 700     #Length DNA handles in bp
    par['DNAds_nm'] = 0.34   #Rise per basepair (nm)
    par['CLds_nm'] = par['CLds_bp'] * par['DNAds_nm'] #Length in nm
    par['DNAss_nm'] = 0.67  #Interbasepair distance for single stranded DNA [nm/bp]
    par['P_nm'] = 50        #Persistencelength dsDNA
    par['Pss_nm'] = 0.67    #Persistencelength ssDNA
    par['S_pN'] = 1000      #Stretching modulus WLC
    par['z0_nm'] = 0        #offset
    par['kBT_pN_nm'] = 4.2  #Kb*T in pn/nm """
    par = {}
    par['CLds_bp'] = 700     #Length DNA handles in bp
    par['DNAds_nm'] = 0.34   #Rise per basepair (nm)
    par['CLds_nm'] = par['CLds_bp'] * par['DNAds_nm'] #Length in nm
    par['P_nm'] = 50        #Persistencelength dsDNA
    par['Pss_nm'] = 0.69    #Persistencelength ssDNA
    par['DNAss_nm'] = 0.67  #Interbasepair distance for single stranded DNA [nm/bp]
    par['S_pN'] = 1000      #Stretching modulus WLC
    par['z0_nm'] = 0        #offset
    par['kBT_pN_nm'] = 4.2  #Kb*T in pn/nm 
    par['MeasurementERR_nm'] = 5   
    return par

def read_data(Filename):
    """Open .dat/.fit files from magnetic tweezers"""
    f = open(Filename, 'r')
    #get headers
    headers = f.readlines()[0]
    headers = headers.split('\t')
    #get data
    f.seek(0) #seek to beginning of the file
    data = f.readlines()[1:]
    f.close()
    F = np.array([])
    Z = np.array([])
    T = np.array([])
    Z_Selected = np.array([])
    for idx,item in enumerate(data):                                            #Get all the data from the fitfile
        F = np.append(F,float(item.split()[headers.index('F (pN)')]))
        Z = np.append(Z,float(item.split()[headers.index('z (um)')])*1000)   
        T = np.append(T,float(item.split()[headers.index('t (s)')]))
        Z_Selected = np.append(Z_Selected,float(item.split()[headers.index('selected z (um)')])*1000)
    return F, Z, T, Z_Selected

def read_sequence(SequenceFile):
    """Reads in a textbased sequence file as a single uppercase string"""
    f = open(SequenceFile, 'r')
    sequence=f.read()
    sequence=sequence.rstrip()
    sequence=sequence.upper()
    f.close()
    return sequence
    
def read_log(Filename):
    """Open the corresponding .log files from magnetic tweezers"""
    f = open(Filename, 'r')
    lines = f.readlines()
    f.close()
    return lines

def log_pars(LogFile):
    """Reads in parameters from the logfile generate by the labview fitting program, returns a {dict} with 'key'= paramvalue"""
    par = {}
    par['L_bp'] = float(find_param(LogFile, 'L DNA (bp)'))
    par['P_nm'] = float(find_param(LogFile, 'p DNA (nm)'))
    par['S_pN'] = float(find_param(LogFile, 'S DNA (pN)'))
    par['degeneracy'] = 0
    par['z0_nm'] = 2
    par['k_pN_nm'] = float(find_param(LogFile, 'k folded (pN/nm)'))
    par['N_tot'] = float(find_param(LogFile, 'N nuc'))
    par['N4'] = float(find_param(LogFile, 'N unfolded [F0]'))
    par['NRL_bp'] = float(find_param(LogFile, 'NRL (bp)'))
    par['ZFiber_nm'] = float(find_param(LogFile, 'l folded (nm)'))
    par['G1_kT'] = 3
    par['G2_kT'] = 4
    par['DNAds_nm'] = 0.34 # rise per basepair (nm)
    par['kBT_pN_nm'] = 4.2 #pn/nm 
    par['Innerwrap_bp'] = 79 #number of basepairs in the inner turn wrap
    par['Fiber0_bp'] = par['L_bp']-(par['N_tot']*par['Innerwrap_bp'])  #Transition between fiber and beats on a string
    par['LFiber_bp'] = (par['N_tot']-par['N4'])*(par['NRL_bp']-par['Innerwrap_bp'])  #total number of bp in the fiber
    par['FiberStart_bp']  = par['Fiber0_bp']-par['LFiber_bp']
    par['MeasurementERR (nm)'] = 5
    return par

def find_param(Logfile, Param):
    """Find a parameter in the .log file"""
    for lines in Logfile:
        P =lines.split(' = ')
        if P[0]==Param:
            return P[1].strip('\n')
    print("<<<<<<<<<<", Param, "not found >>>>>>>>>>")
    return

def handle_data(F, Z, T, Z_Selected, Handles, Window=5):
    """Reads in parameters from the logfile generate by the labview fitting program"""
    if Handles['Select']:                                                       #If only the selected column is use do this
        F_Selected = np.delete(F, np.argwhere(np.isnan(Z_Selected)))
        T_Selected = np.delete(T, np.argwhere(np.isnan(Z_Selected)))
        Z_Selected = np.delete(Z, np.argwhere(np.isnan(Z_Selected))) 
        if len(Z_Selected)==0: 
            print('==> Nothing Selected!')
            return [], [], []
        else:
            F_Selected, Z_Selected, T_Selected = minforce(F_Selected, Z_Selected, T_Selected , Handles['MinForce'])
            return F_Selected, Z_Selected, T_Selected
    else:
        F_Selected = F
        Z_Selected = Z
        T_Selected = T
    
    if Handles['DelBreaks']: F_Selected ,Z_Selected, T_Selected = breaks(F_Selected, Z_Selected, T_Selected, 5000)
    if Handles['Pulling']: F_Selected, Z_Selected, T_Selected = removerelease(F_Selected, Z_Selected, T_Selected )
    if Handles['MinForce'] > 0: F_Selected, Z_Selected, T_Selected = minforce(F_Selected, Z_Selected, T_Selected , Handles['MinForce'])
    if Handles['MedFilt']: Z_Selected = signal.medfilt(Z_Selected, Window)
    return F_Selected, Z_Selected, T_Selected

def breaks(F, Z, T, Jump=1000):
    """Removes the data after a jump in z, presumably indicating the bead broke lose"""
    Test = Z[0]
    for i,x in enumerate(Z[1:]):
        if abs(x - Test) > Jump :
            F = F[:i]
            Z = Z[:i] 
            T = T[:i] 
            break
        Test = x
    return F, Z, T

def removerelease(F, Z, T):
    """Removes the release curve from the selected data"""
    test = 0
    Pullingtest = np.array([])
    for i,x in enumerate(F):
        if x < test:
            Pullingtest = np.append(Pullingtest,i)
        test = x
    F = np.delete(F, Pullingtest)
    Z = np.delete(Z, Pullingtest)
    T = np.delete(T, Pullingtest)
    return F, Z, T 

def minforce(F, Z,  T, Min_Force=2):
    """Removes the data below minimum force given"""
    Mask = F > Min_Force
    Z = Z[Mask]
    F = F[Mask]
    T = T[Mask]
    return F, Z, T

def maxextention(F, Z, T, Max_Extension): 
    """Removes the data above maximum extension given"""
    Mask = Z < Max_Extension
    Z = Z[Mask]
    F = F[Mask]
    T = T[Mask]
    return F, Z ,T


"""
#This function is not used atm
            
def write_data(Filename,Headers,Data):
    f = open(Filename, 'a')
#    import json
#    json.dump(str(Data),f)
    Headers='\t'.join(map(str,Headers))+'\n'
    f.write(Headers)
    Data='\t'.join(map(str,Data))+'\n'
    f.write(Data)
    f.close()
    return "resultsfile generated"
"""
