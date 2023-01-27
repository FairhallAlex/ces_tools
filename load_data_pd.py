#%% Load CES Data

import numpy as np
import pandas as pd
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
    raw = raw data loaded from csv
    transmissionData = APD signal for transmission measurements
    filtTransData = dressed transmission data
    invyDataTrans = inverted transmission data used for peak picking
    reflectionData = APD signal for reflection measurements
    """
    minYTransFiles=[]
    maxYTransFiles=[]
    maxYRefFiles=[]
    for fileIdx in range(len(files)):
        print('On file '+str(fileIdx+1)+' of '+str(len(files))+'.')
        raw = pd.read_csv(files[fileIdx], names=['Time','Transmission','Reflection'])

        transmissionData = raw['Transmission'].loc[:0.9*len(raw['Transmission'])]
        b,c = butter(2, 1200/nyq)
        filtTransData = pd.Series(filtfilt(b,c,transmissionData))

        reflectionData = raw['Reflection'].loc[:0.9*len(raw['Reflection'])]
        b,c = butter(2, 1000/nyq)
        filtRefData = pd.Series(filtfilt(b,c,reflectionData))

        filtTransData.add(1) # Offset to make sure it goes above APD ground offset
        invyDataTrans=-filtTransData.add(filtTransData.max())  #Adding an offset to make normalisation possible, APD has negative minimum offset.

        ## Set reflection data minimum to 0
        if filtRefData.min()<0:
            filtRefData.add(filtRefData.min()*(-1))
        elif filtRefData.min()>0:
            filtRefData.sub(filtRefData.min())

        minYTransFiles.append(invyDataTrans.min())
        maxYTransFiles.append(invyDataTrans.max())
        maxYRefFiles.append(filtRefData.max())
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
def find_files(root,tree=None):
    list_exceptions = ['metadata.csv']
    print('\n'+root)
    _split = os.path.split(root)[1]
    if tree: folders = tree
    else: folders = {_split+' files':[]}

    for fol in os.listdir(root):
        obj = os.path.join(root,fol)
        if fol in list_exceptions: 
            pass
        elif os.path.isdir(obj):
            folders[fol]=find_files(obj)
        elif os.path.splitext(obj)[1].lower() == '.csv':
            folders[_split+' files'].append(obj)

    if folders[_split+' files'] == []:
        _=folders.pop(_split+' files')

    return folders

folderHierarchy = find_files(root=r'C:\Users\goldsmith\Desktop\ces_tools\rootFolder')

