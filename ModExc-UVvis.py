"""
Code for Modulation Excitation analysis of UVvis data
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import math
import multiprocessing

###############################################################################
# Edit scan numbers for the analysis
#Scan_no = [72, 73, 74, 75, 76, 77]; # Edit scan numbers 
Scan_no=np.arange(72, 100+1,1);
Scan_files_in_dir = []; # No need to modify this one

e_col = 0; # Column index with energy
abs_col = 4; # Column index with absorbance

uv='even' # Change to odd if 'even' scan numbers corresponds to excitation with visible light
cycles_for_aver = 3;

phistart=0; # Phase values
phistop=350;
phistep=10;

K=1;            # User-defined K value set to 1
timestep=0.5;   # Time step (s) between acquisition of different spectra during the experiment

FILENAME='0.1 mM D6-D10 980 ms IT 1st Experiment'
###############################################################################
# Data loading 
 
for r, d, f in os.walk(os.getcwd()):
    for indexf in range(0,len(Scan_no)):
        for item in f:
            if (str(Scan_no[indexf])+'_') in item:
                Scan_files_in_dir.append(item)

buf=np.loadtxt(Scan_files_in_dir[0],skiprows=8);
uv_irr_e=np.zeros((np.size(buf,0),np.size(Scan_files_in_dir)))
uv_irr_abs=np.zeros((np.size(buf,0),np.size(Scan_files_in_dir)))
vis_irr_e=np.zeros((np.size(buf,0),np.size(Scan_files_in_dir)))
vis_irr_abs=np.zeros((np.size(buf,0),np.size(Scan_files_in_dir)))

for index in range(0,len(Scan_files_in_dir)):
    buf=np.loadtxt(Scan_files_in_dir[index],skiprows=8);

    if uv=='even':
        if Scan_no[index] % 2 == 0: 
            uv_irr_e[:,index] = buf[:,e_col];
            uv_irr_abs[:,index] = buf[:,abs_col];
        else:
            vis_irr_e[:,index] = buf[:,e_col];
            vis_irr_abs[:,index] = buf[:,abs_col];
            
    elif uv=='odd':
        if Scan_no[index] % 2 == 0:
            vis_irr_e[:,index] = buf[:,e_col];
            vis_irr_abs[:,index] = buf[:,abs_col];
        else:
            uv_irr_e[:,index] = buf[:,e_col];
            uv_irr_abs[:,index] = buf[:,abs_col];                
    else:
        print('Are you sure you gave me the right option for the excitation with UV light?')

uv_irr_e=np.delete(uv_irr_e,np.where(~uv_irr_e.any(axis=0))[0], axis=1)
uv_irr_abs=np.delete(uv_irr_abs,np.where(~uv_irr_abs.any(axis=0))[0], axis=1)
vis_irr_e=np.delete(vis_irr_e,np.where(~vis_irr_e.any(axis=0))[0], axis=1)
vis_irr_abs=np.delete(vis_irr_abs,np.where(~vis_irr_abs.any(axis=0))[0], axis=1)


###############################################################################
# Averaging over a user-given number of cycles

def average_cycles(allscans,numcycles):
    row_size=np.ma.size(allscans,0); 
    col_size=np.int_(np.ma.size(allscans,1)/numcycles); 
    if col_size >= np.int_(col_size): 
        col_size=math.ceil(col_size);
    merged_int=np.ones((row_size, col_size));
    #sz=np.ma.size(allscans,1)
    for index in range(0,col_size): 
        merged_int[:,index] = np.mean(allscans[:,np.arange(index*numcycles, ((index+1)*numcycles),1)],1)
    return merged_int

#cpu_no = int(os.environ['NUMBER_OF_PROCESSORS'])
#if cpu_no >= 4: 
#    PROCESSES=4
#    with multiprocessing.Pool(PROCESSES) as pool:
#        TASKS = [uv_irr_e, cycles_for_aver] + \
#                [uv_irr_abs, cycles_for_aver] + \
#                [vis_irr_e, cycles_for_aver] + \
#                [vis_irr_abs, cycles_for_aver]
#        results = (pool.map(average_cycles, t) for t in TASKS);   
#        print (results)
        #results = [pool.apply_async(average_cycles, t) for t in TASKS]
     
uv_irr_e_mer = average_cycles(uv_irr_e,cycles_for_aver);
uv_irr_abs_mer = average_cycles(uv_irr_abs,cycles_for_aver);
vis_irr_e_mer = average_cycles(vis_irr_e,cycles_for_aver);
vis_irr_abs_mer = average_cycles(vis_irr_abs,cycles_for_aver);

###############################################################################
# Graphical data consistency check
#plt.figure()
#for index in range(0,np.size(aaa,1)-1):
#    plt.plot(uv_irr_e_mer[:,index],uv_irr_abs_mer[:,index])
#    plt.plot(vis_irr_e_mer[:,index],vis_irr_abs_mer[:,index],'.')
#plt.show()
###############################################################################
# Phase modulation
deg=np.arange(phistart,phistop+phistep,phistep); phival=np.radians(deg);

def phase_integr_I(import_data, phival,K,cycles_for_aver,timestep):
    numscans=np.ma.size(import_data,1)
    T=numscans*cycles_for_aver*timestep;
    timestamp=np.arange(1,numscans+1)*timestep
    
    datatimesin=import_data*np.sin((K*2*np.pi/T)*timestamp);
    datatimecos=import_data*np.cos((K*2*np.pi/T)*timestamp);

    intI=np.ones((np.ma.size(import_data,0), len(phival)))
    for index2 in range(0,len(phival)):
        for index1 in range(0,len(import_data)):
            buf=(datatimesin[index1,:]*np.cos(phival[index2])) + (datatimecos[index1,:]*np.sin(phival[index2])) 
            intI[index1,index2]=(2/T)*np.trapz(buf);
    return intI

uv_intI = phase_integr_I(uv_irr_abs_mer, phival,K,cycles_for_aver,timestep);
vis_intI = phase_integr_I(vis_irr_abs_mer, phival,K,cycles_for_aver,timestep);

###############################################################################
# Graphical post-processing

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


################################# UV Excitation
fig = plt.figure()
plt.plot(uv_irr_e_mer[:,0], uv_intI[:,0:-1], '.');
plt.xlabel('Wavelength (nm)'); plt.ylabel('Absorbance Intensity (Arb. Units)'); 
TITLE="Spectra collected with UV excitation"; plt.title(TITLE)
plt.show()
print('Click where you would like to view the modulation analysis for UV data and then press return')
bufsel = plt.ginput(-1)
bufsel = np.asarray(bufsel); Xsel=bufsel[:,0]; Ysel=bufsel[:,1]


B_uv=np.ones((len(phival),len(Xsel)+1)); B_uv[:,0]=deg;
X_index=np.ones((len(Xsel)))
for ind in range(0,len(X_index)):
    X_index[ind] = find_nearest(uv_irr_e_mer[:,0], Xsel[ind])
    B_uv[:,ind+1]=uv_intI[int(X_index[ind]),:];

################################# Vis-light Excitation
fig = plt.figure()
plt.plot(vis_irr_e_mer[:,0], vis_intI[:,0:-1], '.');
plt.xlabel('Wavelength (nm)'); plt.ylabel('Absorbance Intensity (Arb. Units)'); 
TITLE="Spectra collected with Vis-light excitation"; plt.title(TITLE)
plt.show()
print('Click where you would like to view the modulation analysis for vis-light data and then press return')
bufsel = plt.ginput(-1)
bufsel = np.asarray(bufsel); Xsel2=bufsel[:,0]; Ysel2=bufsel[:,1]


B_vis=np.ones((len(phival),len(Xsel2)+1)); B_vis[:,0]=deg;
X_index2=np.ones((len(Xsel2)))
for ind in range(0,len(X_index2)):
    X_index2[ind] = find_nearest(vis_irr_e_mer[:,0], Xsel2[ind])
    B_vis[:,ind+1]=vis_intI[int(X_index2[ind]),:];


###############################################################################
# Phase Modulation graph 
fig2 = plt.figure()
plt.subplot(2,1,1)
for indplot in range(1,np.ma.size(B_uv,1)) :
    val= "%.2f" % uv_irr_e_mer[int(X_index[indplot-1]),0] # Shrink number to 3 digits only
    strlegend="Wavelength = {} nm".format(val)
    plt.plot(B_uv[:,0], B_uv[:,indplot], '.', label=strlegend)

TITLE="UV Irradiation - Mod. angle for the {} user-defined R value(s)".format(len(Xsel))
plt.title(TITLE)
plt.xlabel('Modulation Angle (deg.)'); plt.ylabel('Intensity (Arb. Units)');
plt.legend()

plt.subplot(2,1,2)
for indplot in range(1,np.ma.size(B_vis,1)) :
    val2= "%.2f" % vis_irr_e_mer[int(X_index2[indplot-1]),0] # Shrink number to 3 digits only
    strlegend="Wavelength = {} nm".format(val2)
    plt.plot(B_vis[:,0], B_vis[:,indplot], '.', label=strlegend)

TITLE="Vis-light Irradiation - Mod. angle for the {} user-defined R value(s)".format(len(Xsel2))
plt.title(TITLE)
plt.xlabel('Modulation Angle (deg.)'); plt.ylabel('Intensity (Arb. Units)');
plt.legend()

plt.show()
#plt.close('all')


###############################################################################
# I/O

answer = 'n' 
#answer = input('Should I save your normalised data in a file (yes/no): ')
if answer.lower().startswith("y"):
    FILELABEL=FILENAME + 'UV_PhaseDemod.dat'
    with open(FILELABEL, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows( zip(uv_irr_e_mer[:,0], uv_intI[:,0:-1]) )

    FILELABEL=FILENAME + 'VIS_PhaseDemod.dat'
    with open(FILELABEL, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows( zip(vis_irr_e_mer[:,0], vis_intI[:,0:-1]) )        
    
    FILELABEL=FILENAME + 'UV_IvsPhase.dat'
    with open(FILELABEL, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(B_uv)
    
    FILELABEL=FILENAME + 'VIS_IvsPhase.dat'
    with open(FILELABEL, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(B_vis)
        
    print('Data successfully saved, live long and prosper!')
if answer.lower().startswith("n"):
    print('Job done, live long and prosper!')
    
   


###############################################################################
                # Coded in a quick and dirty way by Matteo (x8498)