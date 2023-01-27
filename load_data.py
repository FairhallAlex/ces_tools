#%% Load CES Data

import numpy as np
import import_csv_data as icd
import os
from scipy.signal import butter, filtfilt, find_peaks, peak_widths, peak_prominences

rootFolder=r'C:\Users\goldsmith\Desktop\ces_tools\rootFolder'
powerFolders = [f.name for f in os.scandir(rootFolder) if f.is_dir()]
saveFiles=True
collectionFrequency=50000
nyq = 0.5*collectionFrequency

def dress_data(files):
    """ 
    Function for loading raw data in CSV format.
    Dresses data for normalization and peak picking.
    Discards last 10% of data to account for locking mechanism.

    Inputs: 
        files = list of data files in current folder
    Outputs:
        minTrans
        maxTrans
        maxRef
    
    i = index of file currently open
    a = index of analyte folder currently open
    x = index of power folder currently open
    filename not used in this loop
    data = raw data loaded from csv
    yDataTrans = APD signal for transmission measurements
    yDataRef = APD signal for reflection measurements
    invyDataTrans = inverted transmission data used for peak picking. 
    """
    minYTransFiles=[]
    maxYTransFiles=[]
    maxYRefFiles=[]
    for fileIdx in range(len(files)):
        print('On file '+str(fileIdx+1)+' of '+str(len(files))+'.')
        data=np.genfromtxt(files[fileIdx], delimiter=',', invalid_raise=False, skip_header=0)

        yDataTrans=data[0:int(len(data)*0.9),1] #transmission data
        b,c = butter(2, 1200/nyq)
        yDataTrans = filtfilt(b,c,yDataTrans)
        # Offset to make sure it goes above APD negative bias
        yDataTrans=yDataTrans+1
        invyDataTrans=-yDataTrans+(np.max(yDataTrans))  #Adding an offset to make normalisation possible, APD has negative minimum offset.

        yDataRef=data[0:int(len(data)*0.9),2] #reflection data
        b,c = butter(2, 1000/nyq)
        yDataRef = filtfilt(b,c,yDataRef)        
        ## Set reflection data minimum to 0
        if np.min(yDataRef)<0:
            yDataRef=yDataRef+(np.min(yDataRef)*-1)
        elif np.min(yDataRef)>0:
            yDataRef=yDataRef-np.min(yDataRef)

        minYTransFiles.append(np.min(yDataTrans))
        maxYTransFiles.append(np.max(invyDataTrans))
        maxYRefFiles.append(np.max(yDataRef))
    return np.min(minYTransFiles), np.max(maxYTransFiles), np.max(maxYRefFiles)

minYTransPowerFolders=[]
maxYTransPowerFolders=[]
maxYRefPowerFolders=[]

for powerIdx in range(len(powerFolders)):
    print(f'On power folder {powerIdx+1} of {len(powerFolders)}.')
    analyteFolders = [f.name for f in os.scandir(rootFolder+'/'+powerFolders[powerIdx]) if f.is_dir()]
    minTransAnalyteFolder=[]
    maxTransAnalyteFolder=[]
    maxRefAnalyteFolder=[]

    path = rootFolder+'/'+powerFolders[powerIdx]+'/'
    for analyteIdx in range(len(analyteFolders)):   
        print(f'\nOn analyte folder {analyteFolders[analyteIdx]}, ({analyteIdx+1} of {len(analyteFolders)}).') 
        files=icd.fn_import_csv(rootFolder+'/'+powerFolders[powerIdx]+'/'+analyteFolders[analyteIdx])
        minTrans,maxTrans,maxRef = dress_data(files)
        minTransAnalyteFolder.append(minTrans)
        maxTransAnalyteFolder.append(maxTrans)
        maxRefAnalyteFolder.append(maxRef)

    minYTransPowerFolders.append(np.min(minTransAnalyteFolder))
    maxYTransPowerFolders.append(np.max(maxTransAnalyteFolder))
    maxYRefPowerFolders.append(np.max(maxRefAnalyteFolder))

# %%
