# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 17:58:45 2017
@author: nhermans & jheinsman
Files needed:
.fit force extension of unzipping curve (corrected for offset + drift)
sequence of unzipped DNA
"""

#Imported libraries

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import Tools
import Functions as func

#Import Data
file_location="G:\\test\\"
file_name="FC2_pBlue_05ul_eh_good_data_016_39"
SequenceFile="G:\\klaas\\Unzipping\\18S_seq.txt"
file_extension=".fit"
FilePath=file_location+file_name+file_extension

sequence=Tools.read_sequence(SequenceFile)
Pars = Tools.default_pars()
Force, Z, Time, Extension_nm = Tools.read_data(FilePath)
Handles = Tools.Define_Handles(Select=True, Pull=True, DelBreaks=True, MinForce=2.5, MinZ=0, MedFilt=False)
Force, Extension_nm, Time = Tools.handle_data(Force, Z, Time, Extension_nm, Handles, Pars)

GCcontent = func.GC_dG(sequence, window = 50) #Calls the function for calcultion of the GC content

#split data in unzipping and annealing
Force_df = Force - np.roll(Force,1)
timeup = Time[Force_df >= 0]
timedown = Time[Force_df < 0]
extensionup = Extension_nm[Force_df >= 0]
extensiondown = Extension_nm[Force_df < 0]
forceup = Force[Force_df >= 0]
forcedown = Force[Force_df < 0]

#Fit ssDNA persistence length
MaxForce=np.max(Force)
print(MaxForce)
if MaxForce < 25: MaxForce=25
Fit_Z = extensiondown[forcedown >= 20]
Fit_F = forcedown[forcedown >= 20]
Fit_Z = Fit_Z[Fit_F < MaxForce]
Fit_F = Fit_F[Fit_F < MaxForce]

popt = curve_fit(lambda f, p: func.fit_p(f, p, len(sequence)), Fit_F, Fit_Z, p0 = Pars['Pss_nm'])

print(popt)
Pars['Pss_nm'] = popt[0][0]

#Histogram data: change extension to bp
CLss_down_bp = (extensiondown-Pars['DNAds_nm']*Pars['CLds_bp']*func.wlc(forcedown, Pars))/(Pars['DNAss_nm']*2*func.fjc(forcedown, Pars))
CLss_up_bp = (extensionup-Pars['DNAds_nm']*Pars['CLds_bp']*func.wlc(forceup, Pars))/(Pars['DNAss_nm']*2*func.fjc(forceup, Pars))

#FE and histogram graph generator
#Legend creator for FE graph
a  = plt.scatter(extensionup, forceup, color = 'blue', marker = 'o', s = 1, label = 'Pull')
b  = plt.scatter(extensiondown, forcedown, color = 'green', marker = 'o', s = 1, label = 'Relax')
c, = plt.plot(func.wlc(Force, Pars)*Pars['DNAds_nm']*Pars['CLds_bp'],Force, color = 'black', linewidth=1.0, label = "WLC ")
d, = plt.plot(func.wlc(Force, Pars)*Pars['DNAds_nm']*Pars['CLds_bp']+func.fjc(Force, Pars)*Pars['DNAss_nm']*len(sequence)*2, Force, color = 'red',linewidth=1.0, label = 'DWLC')

#Graph plot parameters for FE graph
plt.tick_params(direction = 'in', top = 'on', right = 'on')
plt.xlabel('Extension [nm]')
plt.ylabel('Force [pN]')
plt.axis([0, 3000, 0, 25])
plt.legend(handles=[a, b, c, d])
plt.savefig('170707_63_FE_DWLC.pdf') #save file name for FE graph
plt.show()

#Contour length histogram generation
basepairs = list(range(0, len(sequence)-2)) #for plotting only

#Contour length histogram + GC-content creator
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.hist(CLss_up_bp, bins = int(len(sequence) / 20), color = 'blue', range = [0,len(sequence)] )
ax1.hist(CLss_down_bp, bins = int(len(sequence) / 20), color = 'green', range = [0,len(sequence)])
ax2.plot(basepairs,GCcontent, 'r-', label = 'dG based on sequence', linewidth = 0.5)

#Graph plot parameters for histogram
ax1.tick_params(top = 'on', direction = 'in')
ax2.tick_params(direction = 'in')
ax1.set_xlabel('# of unzipped basepairs')
ax1.set_xlim((0,len(sequence)))
ax1.set_yscale('symlog')
#ax1.set_ylim((1,10**4))
ax1.set_ylabel('Count', color='black')
ax2.set_ylabel('deltaG (kT)', color='black')
ax2.set_ylim(1,1.8)
plt.savefig('170707_63_GC_DWLC.pdf') #save file name for histogram + GC-content grapgh
plt.legend(loc=1)
plt.show()

FilePath=file_location+file_name+"_pull.dat"
file = open(FilePath, "w")

for i in range(0,len(timeup)):
    file.write("%s\t" %timeup[i])#,forceup[i],extensionup[i])
    file.write("%s\t" %forceup[i])
    file.write("%s\n" %extensionup[i])

file.close()

FilePath=file_location+file_name+"_release.dat"
file = open(FilePath, "w")

for i in range(0,len(timedown)):
    file.write("%s\t" %timedown[i])
    file.write("%s\t" %forcedown[i])
    file.write("%s\n" %extensiondown[i])

file.close()

FilePath=file_location+file_name+"_hist.dat"
file = open(FilePath, "w")

for i in range(0,len(CLss_up_bp)):
    a=0
    file.write("%s\t" % i)
    if CLss_up_bp[i]>0 and CLss_up_bp[i]<5000:
        file.write("%s\t" % CLss_up_bp[i])#,CLss_down_bp[i],GCcontent[i])
        a=1
    if i<len(CLss_down_bp) and CLss_down_bp[i]>0 and CLss_down_bp[i]<5000:
        if a==0: file.write("\t")
        file.write("%s\t" % CLss_down_bp[i])
        a=1
    if i<len(GCcontent):
        if a==0: file.write("\t \t")
        file.write("%s" % GCcontent[i])
        a=1
    if a==1:
        file.write("\n")
file.close()