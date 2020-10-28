# -*- coding: utf-8 -*-
"""
Code for Simpson's demodulation analysis of (in this case) \Chi(R) EXAFS data
"""
###############################################################################
# User's Setup

FILENAME='34570_Pd_Al2O3_H2xO2_300C_crop_rebin_crop_k3_FFT_k3_576-2879_MES.dat';

phistart=0;
phistop=350;
phistep=10;

K=1;            # User-defined K value set to 1
timestep=0.5;   # Time step (s) between acquisition of different spectra during the experiment
numframes=144;  # Number of frames in single cycle (144 for H2xO2 and 288 for CH4_CO2 CO)


###############################################################################

import numpy as np
import matplotlib.pyplot as plt
import csv, itertools

data_buffer = np.loadtxt(FILENAME)
energies=data_buffer[:,0];


def period_merger(allscans,numframes):
    merged_int=np.ones((np.ma.size(allscans,0),numframes))
    sz=np.ma.size(allscans,1)
    for index in range(0,numframes-1): 
        merged_int[:,index]=np.mean(allscans[:,np.arange(index,sz,numframes)],1)
    return merged_int

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx



import_data=period_merger(data_buffer[:,1:np.ma.size(data_buffer,1)], numframes); # All scans included

numscans=np.ma.size(import_data,1)    # Number of scans over the desired integration period T
deg=np.arange(phistart,phistop+phistep,phistep); phival=np.radians(deg)
timestamp=np.arange(1,numscans+1)*timestep
T=numscans*timestep

datatimesin=import_data*np.sin((K*2*np.pi/T)*timestamp);
datatimecos=import_data*np.cos((K*2*np.pi/T)*timestamp);

intI=np.ones((np.ma.size(import_data,0), len(phival)))
for index2 in range(0,len(phival)):
    for index1 in range(0,len(import_data)):
        buf=(datatimesin[index1,:]*np.cos(phival[index2])) + (datatimecos[index1,:]*np.sin(phival[index2])) 
        intI[index1,index2]=(2/T)*np.trapz(buf);
        #return intI


###############################################################################
# Graphical post-processing

fig = plt.figure()
plt.plot(energies, intI[:,0:-1], '.')
plt.xlabel('R ($\AA$)'); plt.ylabel('Intensity (Arb. Units)'); plt.show()
print('Now click where you would like to view the modulation analysis and then press return')
bufsel =plt.ginput(-1)
bufsel=np.asarray(bufsel); Xsel=bufsel[:,0]; Ysel=bufsel[:,1]


B=np.ones((len(phival),len(Xsel)+1)); B[:,0]=deg;
X_index=np.ones((len(Xsel)))
for ind in range(0,len(X_index)):
    X_index[ind] = find_nearest(energies, Xsel[ind])
    B[:,ind+1]=intI[int(X_index[ind]),:];
    #subplot(2,1,1)
    #plot(phival*(180/pi),intI(X_index(ind),:));

fig = plt.figure()
plt.subplot(2,1,1)
for indplot in range(1,np.ma.size(B,1)) :
    val= "%.3f" % energies[int(X_index[indplot-1])] # Shrink number to 3 digits only
    strlegend="R = {}".format(val)
    plt.plot(B[:,0], B[:,indplot], '.', label=strlegend)

TITLE="Demodulation angle for the {} user-defined R value(s)".format(len(Xsel))
plt.title(TITLE)
plt.xlabel('Demodulation Angle (deg.)'); plt.ylabel('Intensity (Arb. Units)');
plt.legend()
plt.show()

###############################################################################
# Threshold analysis
#thrs=input('Define Signal Threshold (from 0 to 1, corrisp to 0% to 100%): ')
thrs=0.33; # Only allowing values from 0 to 1 (namely 0% to 100%)

if thrs<0:
    thrs=0
elif thrs>1:
    thrs=1


C=[]
refI=thrs*np.max(abs(intI))

for ind in range(0,np.size(intI,1)): 
    indexAbove=np.where(abs(intI[:,ind]) > refI)
    plt.subplot(2,1,2)
    if len(indexAbove)>0:
       yabove=np.ones(np.size(indexAbove,1))*deg[ind]
       buff=[energies[indexAbove],yabove]; C=np.append(C,buff) 
       plt.plot(energies[indexAbove],yabove,'bo')
       
TITLE="Phase Value above Absolute Intensity threshold of {} %".format(thrs*100) 
plt.title(TITLE);
plt.xlim((np.min(energies), np.max(energies))); plt.ylim((-5, 365));
plt.xlabel('R ($\AA$)'); plt.ylabel('Phase Value (deg.)'); 
plt.show()       


###############################################################################
# I/O

ioval=input('Should I save the data? (y/n) ');
if ioval is 'y':
    FILELABEL=FILENAME + '_PhaseDemod.dat'
    with open(FILELABEL, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows( zip(energies,intI) )
    
    FILELABEL=FILENAME + '_IvsPhase.dat'
    with open(FILELABEL, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(B)
    
    FILELABEL=FILENAME + '_PhasevsR.dat'
    np.savetxt(FILELABEL,C)
        
    print('Data successfully saved, live long and prosper!')
elif ioval == 'n':
    print('Job done, live long and prosper!')
    
   


###############################################################################
                # Coded in a quick and dirty way by Matteo (x8498)