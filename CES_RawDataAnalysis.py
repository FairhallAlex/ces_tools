#folder structure- main folder containing all data->Exps and controls as subfolders-> comparable datasets
#rootFolder= main folder with all data in it
#parentFolders=controls or data folders
#subfolders=different conditions within parentFolders folder.
#User defined variables

"""Cell 1-import libraries"""
import numpy as np
import matplotlib.pyplot as plt
import import_csv_data as icd
import os
from pylab import text
from scipy import signal
from scipy.signal import find_peaks, peak_widths, peak_prominences
import csv
from matplotlib.ticker import FormatStrFormatter


rootFolder=''

saveFiles=True
relPeakAmpCoefTrans=0.5 #Value between 0-1
relPeakAmpCoefRef=0.5 #Value between 0-1
collectionFrequency=50000
peakHeight=0.5
figureFolder=''
timePeakFilter=250

#Code starts here, don't mess with it.
parentFolders = [f.name for f in os.scandir(rootFolder) if f.is_dir()]

#Trans= transmission, Ref= reflection
filename_list=[]
IntMean=[]
peakAmpTransMean=[]
peakAmpTransStd=[]
peakAmpRefMean=[]
peakAmpRefStd=[]
NPeaksTrans=[]
NPeaksRef=[]
minYTransParentFolders=[]
maxYTransParentFolders=[]
minYRefParentFolders=[]
maxYRefParentFolders=[]

plt.rcParams.update({'figure.max_open_warning': 0})

#Data analysis bit
print('Part one-Normalization')
for x in range(len(parentFolders)):
    subfolders = [f.name for f in os.scandir(rootFolder+'/'+parentFolders[x]) if f.is_dir()]
    minTransSubfolder=[]
    maxTransSubfolder=[]
    minRefSubfolder=[]
    maxRefSubfolder=[]
    for a in range(len(subfolders)):
        for fname in os.listdir(rootFolder+'/'+parentFolders[x]+'/'+subfolders[a]):
            if not any(fname.endswith('.csv') for fname in os.listdir(rootFolder+'/'+parentFolders[x]+'/'+subfolders[a]))and a<len(subfolders):
                a=a+1
            elif not any(fname.endswith('.csv') for fname in os.listdir(rootFolder+'/'+parentFolders[x]+'/'+subfolders[a])) and a==len(subfolders):
                break
        files=icd.fn_import_csv(rootFolder+'/'+parentFolders[x]+'/'+subfolders[a])
        if saveFiles==True:
            if not os.path.exists(rootFolder+'/'+parentFolders[x]+'/'+subfolders[a]+'/Plots'+figureFolder):
                os.makedirs(rootFolder+'/'+parentFolders[x]+'/'+subfolders[a]+'/Plots'+figureFolder)
        minYTransFiles=[]
        minYRefFiles=[]
        maxYTransFiles=[]
        maxYRefFiles=[]
        #meanAmp=[]
        #meanWidth=[]
        #stdAmp=[]
        #stdWidth=[]
        #nPeaks=[]
        
        """
        This first for loop is to calculate the minimum and maximum values in every dataset.
        These values are then used for a normalization step, which enables relative comparison of peak ampltiudes across comparable datsets 
        """
        for i in range(len(files)):
            print('On file '+str(i+1)+' of '+str(len(files))+' of subfolder '+str(a+1)+' of '+str(len(subfolders))+' and of parentfolder '+str(x+1)+' of '+str(len(parentFolders)))
            filename_plusext = os.path.basename(files[i])
            filename=os.path.splitext(filename_plusext)[0]
            filename_list.append(os.path.splitext(filename_plusext)[0])
            data=np.genfromtxt(files[i], delimiter=',', dtype=np.str, skip_header=0)
            data= data.astype(np.float)
            xData=data[0:int(len(data)*0.9),0] 
            yDataTrans=data[0:int(len(data)*0.9),1] #transmission data
            nyq = 0.5 * (collectionFrequency)
            b,c = signal.butter(2, 1200/nyq)
            yDataTrans = signal.filtfilt(b,c,yDataTrans)
            yDataRef=data[0:int(len(data)*0.9),2] #reflection data
            b,c = signal.butter(2, 1000/nyq)
            yDataRef = signal.filtfilt(b,c,yDataRef)
            yDataTrans=yDataTrans+1
            invyDataTrans=-yDataTrans+(np.max(yDataTrans))  #Adding an offset to make normalisation possible, APD has negative minimum offset.
            if np.min(yDataRef)<0:
                yDataRef=yDataRef+(np.min(yDataRef)*-1)
            elif np.min(yDataRef)>0:
                yDataRef=yDataRef-np.min(yDataRef)
            resolution=xData[1]-xData[0] #temporal interval in seconds
            #yData=yData-np.min(yData)
            minYTransFiles.append(np.min(yDataTrans))
            maxYTransFiles.append(np.max(invyDataTrans))
            #minYRefFiles.append(np.min(yDataRef))
            maxYRefFiles.append(np.max(yDataRef))

        minTransSubfolder.append(np.min(minYTransFiles))
        maxTransSubfolder.append(np.max(maxYTransFiles))
        #minRefSubfolder.append(np.min(minYRefFiles))
        maxRefSubfolder.append(np.max(maxYRefFiles))

    minYTransParentFolders.append(np.min(minTransSubfolder))
    maxYTransParentFolders.append(np.max(maxTransSubfolder))
    #minYRefParentFolders.append(np.min(minRefSubfolder))
    maxYRefParentFolders.append(np.max(maxRefSubfolder))

"""
This is where the processing starts
"""    
print('Part two-Analysis and plotting')
for x in range(len(parentFolders)):
    subfolders = [f.name for f in os.scandir(rootFolder+'/'+parentFolders[x]) if f.is_dir()]
    for a in range(len(subfolders)):
        for fname in os.listdir(rootFolder+'/'+parentFolders[x]+'/'+subfolders[a]):
            if not any(fname.endswith('.csv') for fname in os.listdir(rootFolder+'/'+parentFolders[x]+'/'+subfolders[a]))and a<len(subfolders):
                a=a+1
            elif not any(fname.endswith('.csv') for fname in os.listdir(rootFolder+'/'+parentFolders[x]+'/'+subfolders[a])) and a==len(subfolders):
                break
                
        files=icd.fn_import_csv(rootFolder+'/'+parentFolders[x]+'/'+subfolders[a])
        peakAmpTransAll=[]
        peakWidthTransAll=[]
        yDataTransAll=[]
        peakAmpRefAll=[]
        peakWidthRefAll=[]
        yDataRefAll=[]
        
        for i in range(len(files)):
            print('On file '+str(i+1)+' of '+str(len(files))+' of subfolder '+str(a+1)+' of '+str(len(subfolders))+' and of parentfolder '+str(x+1)+' of '+str(len(parentFolders)))
            filename_plusext = os.path.basename(files[i])
            filename=os.path.splitext(filename_plusext)[0]
            filename_list.append(os.path.splitext(filename_plusext)[0])
            data=np.genfromtxt(files[i], delimiter=',', dtype=np.str, skip_header=0)
            data= data.astype(np.float)
            xData=data[0:int(len(data)*0.9),0] 
            yDataTrans=data[0:int(len(data)*0.9),1] #transmission data, multiply aprot(27-11-22) and cmyc (26-11-22) data by 0.87 to correct for transmission APD gain change 
            nyq = 0.5 * (collectionFrequency) #Nyquist sampling
            b,c = signal.butter(2, 1200/nyq) #2500 is default
            yDataTrans = signal.filtfilt(b,c,yDataTrans)
            yDataRef=data[0:int(len(data)*0.9),2] #reflection data
            b,c = signal.butter(2, 1000/nyq)
            yDataRef = signal.filtfilt(b,c,yDataRef)
            if np.min(yDataRef)<0:
                yDataRef=yDataRef+(np.min(yDataRef)*-1)
            elif np.min(yDataRef)>0:
                yDataRef=yDataRef-np.min(yDataRef)
            resolution=xData[1]-xData[0]

            yDataTransAll.extend(yDataTrans)
            yDataRefAll.extend(yDataRef)
            
            #Subplot 0, plot raw un-normalised data
            fig,axes=plt.subplots(2,3,figsize=(23, 10))
            fig.suptitle(str(filename))
            axes[0,0].plot(xData, yDataTrans, color='teal')
            axes[0,0].set_ylim((np.min(yDataTrans)-0.1),(np.max(yDataTrans)*1.1))
            axes[0,0].set_xlabel('Time (s)')
            axes[0,0].set_ylabel('APD voltage (V)')
            axes[0,0].set_title('Raw Data')
            axes[0,0].legend(['Transmission'])

            axes[1,0].plot(xData, yDataRef, color='teal')
            axes[1,0].set_ylim((np.min(yDataRef)*0.95),(np.max(yDataRef)*1.1))
            axes[1,0].set_xlabel('Time (s)')
            axes[1,0].set_ylabel('APD voltage (V)')
            axes[1,0].legend(['Reflection'])

            yDataTrans=yDataTrans+1
            globalNormYTrans=yDataTrans/minYTransParentFolders[x]
            axes[0,1].plot(xData, globalNormYTrans, color='teal')
            axes[0,1].set_ylim(0.8, (np.max(globalNormYTrans))*1.01)
            axes[0,1].set_xlabel('Time (s)')
            axes[0,1].set_ylabel('Normalized intensity (a.u.)')
            axes[0,1].set_title('Normalized Data')

            globalNormYRef=yDataRef/maxYRefParentFolders[x]
            axes[1,1].plot(xData, globalNormYRef, color='teal')
            axes[1,1].set_ylim(np.min(globalNormYRef)*0.95, (np.max(globalNormYRef))*1.1)
            axes[1,1].set_xlabel('Time (s)')
            axes[1,1].set_ylabel('Normalized intensity (a.u.)')
            axes[1,1].set_ylim(0, 1.01)
            
            #Pick transmission dips on globally normalised data
            invyDataTrans=-yDataTrans+(np.max(yDataTrans))
            invyDataTransNorm=invyDataTrans/maxYTransParentFolders[x]
            #need this later for extracting 
            relPeakAmpTrans=relPeakAmpCoefTrans
            peaksTrans,_ = find_peaks(invyDataTransNorm,prominence=relPeakAmpTrans,distance=timePeakFilter, width=0, rel_height=peakHeight)

            #Pick reflection peaks on globally normalised data
            relPeakAmpRef=relPeakAmpCoefRef
            peaksRef,_ = find_peaks(globalNormYRef,prominence=relPeakAmpRef,distance=timePeakFilter, width=0, rel_height=peakHeight)
        
            #plot labelled peaks and thresholds
            axes[0,2].plot(xData,yDataTrans,color='cornflowerblue')
            axes[0,2].plot(xData[peaksTrans],(yDataTrans[peaksTrans]),'x', color='orange')
            axes[0,2].set_ylabel('APD voltage (V)')
            axes[0,2].set_xlabel('Time (s)')
            axes[0,2].set_ylim((np.min(yDataTrans)-0.1),(np.max(yDataTrans)*1.1))
            #axes[0,1].text(np.min(xData)*1.05,np.max(yDataTrans)*1.05,str(np.round(len(peakAmpTrans)/np.max(xData),3))+'peaks/s')
            axes[0,2].set_title('Picked peaks')

            axes[1,2].plot(xData,yDataRef,color='cornflowerblue')
            axes[1,2].plot(xData[peaksRef],(yDataRef[peaksRef]),'x', color='orange')
            axes[1,2].set_ylabel('APD voltage (V)')
            axes[1,2].set_xlabel('Time (s)')
            axes[1,2].set_ylim((np.min(yDataRef)*0.95),(np.max(yDataRef)*1.1))
            #axes[1,1].text(np.min(xData)*1.05,np.max(yDataRef)*0.95,str(np.round(len(peakAmpRef)/np.max(xData),3))+'peaks/s')

            plt.rcParams['font.sans-serif'] = ['Arial']
            plt.rcParams.update({'font.size': 16})
            
            if saveFiles==True:
                plt.savefig(rootFolder+'/'+parentFolders[x]+'/'+subfolders[a]+'/Plots'+figureFolder +'/'+ filename +'_Raw_Norm_PickedPeaks.pdf', dpi=200,bbox_inches='tight')
            
            if len(peaksTrans)>2 and len(peaksRef)<2:
                cellWidthsTrans=peak_widths(invyDataTrans,peaksTrans)[0]
                peakWidthsTrans=cellWidthsTrans*resolution
                peakAmpTrans=peak_prominences(invyDataTrans,peaksTrans)[0]
                peakAmpTransAll.extend(peakAmpTrans)
                peakWidthTransAll.extend(peakWidthsTrans)

            elif len(peaksTrans)>2 and len(peaksRef)>2:
            
                cellWidthsTrans=peak_widths(invyDataTrans,peaksTrans)[0]
                peakWidthsTrans=cellWidthsTrans*resolution
                peakAmpTrans=peak_prominences(invyDataTrans,peaksTrans)[0]
                
                cellWidthsRef=peak_widths(yDataRef,peaksRef)[0]
                peakWidthsRef=cellWidthsRef*resolution
                peakAmpRef=peak_prominences(yDataRef,peaksRef)[0]

                peakAmpTransAll.extend(peakAmpTrans)
                peakWidthTransAll.extend(peakWidthsTrans)
                peakAmpRefAll.extend(peakAmpRef)
                peakWidthRefAll.extend(peakWidthsRef)


        """plot histogram of peak ampltitudes across all comparable data sets"""
        if len(peakAmpTransAll)>2 and len(peakAmpRefAll)<2:
            minwidth=(np.min(peakWidthTransAll))*0.8
            maxwidth=(np.max(peakWidthTransAll))*1.2
            minAmp=(np.min(peakAmpTransAll))*0.8
            maxAmp=(np.max(peakAmpTransAll))*1.2
            
            fig,axes=plt.subplots(1,3,figsize=(20, 5))

            "transmission prominence histogram"
            #fig.suptitle(str(parentFolders[x]+subfolders[a]+'_prominence and width histograms'),x=0.5, y=1.00005)
            fig.suptitle(str(parentFolders[x]+subfolders[a]+'_prominence and width histograms'),x=0.5, y=1)
            numBins = np.linspace(minAmp,maxAmp,80)
            counts,histbins=np.histogram(peakAmpTransAll, bins=numBins)
            #counts = counts / np.max(counts)
            axes[0].hist(histbins[:-1], histbins, weights=counts,color='darkorange')
            axes[0].set_ylabel('Counts')
            axes[0].set_xlabel('Peak prominence (V)')
            axes[0].legend(['Transmitted peak prominence'])
            axes[0].set_title(str(subfolders[a]+'_Prominence Histogram'))

            """plot histogram of peak widths for all comparable datasets """ 
            "Transmission width histogram"
            #peakWidthTransAll=peakWidthTransAll*1000
            numBins = np.linspace(minwidth,maxwidth,80)
            counts,histbins=np.histogram(peakWidthTransAll, bins=numBins)
            #counts = counts / np.max(counts)
            axes[1].hist(histbins[:-1], histbins, weights=counts,color='darkgreen')
            axes[1].set_ylabel('Counts')
            axes[1].set_xlabel('Peak width (s)')
            axes[1].xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
            axes[1].legend(['Transmitted peak width'])
            #plt.yscale('log')
            #axes[0,1].set_ylim(np.min(peakWidthTransAll)-(np.min(peakWidthTransAll)*0.95),np.max(peakWidthTransAll)*1.05)
            axes[1].set_title(str(subfolders[a]+'_Width Histograms'))

            axes[2].scatter(peakWidthTransAll, peakAmpTransAll)
            axes[2].set_ylabel('Peak prominence (V)')
            axes[2].set_xlabel('Peak width (s)')
            axes[2].xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
            axes[2].set_title(str(subfolders[a]+'AmpvsWidth'))
            axes[2].set_ylim([np.min(peakAmpTransAll)*0.9,np.max(peakAmpTransAll)*1.1])
            axes[2].set_xlim([np.min(peakWidthTransAll)*0.5,np.max(peakWidthTransAll)*1.1])
            
            plt.rcParams['font.sans-serif'] = ['Arial']
            plt.rcParams.update({'font.size': 16})
            if saveFiles==True:
                plt.savefig(rootFolder+'/'+parentFolders[x]+'/'+subfolders[a]+'/Plots'+figureFolder+ '/' + filename +'_Width_Prominence_Histograms_Trans.pdf', dpi=200,bbox_inches='tight')
            
            with open(rootFolder+'/'+parentFolders[x]+'/'+subfolders[a]+'_MeanParams_HighInt.txt', 'w',encoding='UTF8', newline='') as f:  
                f.writelines(['Mean trans Width= '+str(round(np.mean(peakWidthTransAll),5))+' ± '+str(round(np.std(peakWidthTransAll),5))+u'\u00B1'+' Mean trans Amplitude= '+str(round(np.mean(peakAmpTransAll),5))+' ± '+str(round(np.std(peakAmpTransAll),5))+u'\u00B1'+' N_peaks/file= '+str(len(peakAmpTransAll)/len(files))])

        elif len(peakAmpTransAll)>2 and len(peakAmpRefAll)>2:

            minwidth=(np.min(peakWidthTransAll)*0.01,np.min(peakWidthRefAll)*0.5)
            maxwidth=(np.max(peakWidthTransAll)*1.2,np.max(peakWidthRefAll)*1.2)
            minAmp=(np.min(peakAmpTransAll)*0.8,np.min(peakAmpRefAll)*0.8)
            maxAmp=(np.max(peakAmpTransAll)*1.2,np.max(peakAmpRefAll)*1.2)

            fig,axes=plt.subplots(2,3,figsize=(25, 10))

            "transmission prominence histogram"
            #fig.suptitle(str(parentFolders[x]+subfolders[a]+'_prominence and width histograms'),x=0.5, y=1.00005)
            fig.suptitle(str(parentFolders[x]+subfolders[a]+'_prominence and width histograms'),x=0.5, y=1)
            #numBins = np.linspace((np.min(peakAmpTransAll)),(np.max(peakAmpTransAll)), int(np.sqrt(len(peakAmpTransAll))*2))
            numBins = np.linspace(minAmp[0],maxAmp[0],80)
            TransAmpcounts,TransAmphistbins=np.histogram(peakAmpTransAll, bins=numBins)
            #TransAmpcounts = TransAmpcounts / np.max(TransAmpcounts)
            axes[0,0].hist(TransAmphistbins[:-1], TransAmphistbins, weights= TransAmpcounts,color='darkorange')
            axes[0,0].set_ylabel('Counts')
            axes[0,0].set_xlabel('Peak prominence (V)')
            axes[0,0].legend(['Transmitted peak prominence'])
            axes[0,0].set_title(str(subfolders[a]+'_Prominence Histogram'))
            axes[0,0].set_xlim(minAmp[0],maxAmp[0])

            "reflection prominence histogram"
            #numBins = np.linspace((np.min(peakAmpRefAll)),(np.max(peakAmpRefAll)), int(np.sqrt(len(peakAmpRefAll))*2))
            numBins = np.linspace(minAmp[1],maxAmp[1],80)
            RefAmpcounts,RefAmphistbins=np.histogram(peakAmpRefAll, bins=numBins)
            #RefAmpcounts = RefAmpcounts / np.max(RefAmpcounts)
            axes[1,0].hist(RefAmphistbins[:-1], RefAmphistbins, weights=RefAmpcounts,color='darkorange')
            axes[1,0].set_ylabel('Counts')
            axes[1,0].set_xlabel('Peak prominence (V)')
            axes[1,0].legend(['Reflected peak prominence'])
            axes[1,0].set_xlim(minAmp[1],maxAmp[1])
            #axes[1,0].set_title(str(parentFolders[x]+subfolders[a]+'_prominenceHistogram'),x=0.5, y=1.1)
            
            """plot histogram of peak widths for all comparable datasets """ 
            "Transmission width histogram"
            #peakWidthTransAll=peakWidthTransAll*1000
            #numBins = np.linspace((np.min(peakWidthTransAll)),(np.max(peakWidthTransAll)), int(np.sqrt(len(peakWidthTransAll))*2))
            numBins = np.linspace(minwidth[0],maxwidth[0],80)
            TransWidthcounts,TransWidthhistbins=np.histogram(peakWidthTransAll, bins=numBins)
            #TransWidthcounts = TransWidthcounts / np.max(TransWidthcounts)
            axes[0,1].hist(TransWidthhistbins[:-1], TransWidthhistbins, weights=TransWidthcounts,color='darkgreen')
            axes[0,1].set_ylabel('Counts')
            axes[0,1].set_xlabel('Peak width (s)')
            axes[0,1].xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
            axes[0,1].legend(['Transmitted peak widths'])
            #axes[0,1].set_ylim(np.min(peakWidthTransAll)-(np.min(peakWidthTransAll)*0.95),np.max(peakWidthTransAll)*1.05)
            axes[0,1].set_title(str(subfolders[a]+'_Width Histograms'))
            axes[0,1].set_xlim(minwidth[0],maxwidth[0])

            "Reflection width histogram"
            #peakWidthRefAll=peakWidthRefAll*1000
            #numBins = np.linspace((np.min(peakWidthRefAll)),(np.max(peakWidthRefAll)), int(np.sqrt(len(peakWidthRefAll))*2))
            numBins = np.linspace(minwidth[1],maxwidth[1],80)
            RefWidthcounts,RefWidthhistbins=np.histogram(peakWidthRefAll, bins=numBins)
            #RefWidthcounts = RefWidthcounts / np.max(RefWidthcounts)
            axes[1,1].hist(RefWidthhistbins[:-1], RefWidthhistbins, weights=RefWidthcounts,color='darkgreen')
            axes[1,1].set_ylabel('Counts')
            axes[1,1].set_xlabel('Peak width (s)')
            axes[1,1].legend(['Reflected peak widths'])
            axes[1,1].xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
            #plt.yscale('log')
            #axes[1,1].set_ylim(np.min(peakWidthRefAll)-(np.min(peakWidthRefAll)*0.95),np.max(peakWidthRefAll)*1.05)
            axes[1,1].set_xlim(minwidth[1],maxwidth[1])

            axes[0,2].scatter(peakWidthTransAll, peakAmpTransAll,s=2)
            axes[0,2].set_ylabel('Peak prominence (V)')
            axes[0,2].set_xlabel('Peak width (s)')
            axes[0,2].xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
            axes[0,2].set_title(str(subfolders[a]+'_AmpvsWidth'))
            axes[0,2].set_ylim([np.min(peakAmpTransAll)*0.9,np.max(peakAmpTransAll)*1.1])
            #axes[0,2].set_xlim([np.min(peakWidthTransAll)*0.5,np.max(peakWidthTransAll)*1.1])
            axes[0,2].set_xlim(minwidth[0],maxwidth[0])
            axes[0,2].set_ylim(minAmp[0],maxAmp[0])
            
            axes[1,2].scatter(peakWidthRefAll, peakAmpRefAll,s=2)
            axes[1,2].set_ylabel('Peak prominence (V)')
            axes[1,2].set_xlabel('Peak width (s)')
            axes[1,2].xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
            axes[1,2].set_ylim([np.min(peakAmpRefAll)*0.8,np.max(peakAmpRefAll)*1.1])
            #axes[1,2].set_xlim([np.min(peakWidthRefAll)*0.7,np.max(peakWidthRefAll)*1.1])
            axes[1,2].set_xlim(minwidth[1],maxwidth[1])
            axes[1,2].set_ylim(minAmp[1],maxAmp[1])
        
            plt.rcParams['font.sans-serif'] = ['Arial']
            plt.rcParams.update({'font.size': 16})
            if saveFiles==True:
                plt.savefig(rootFolder+'/'+parentFolders[x]+'/'+subfolders[a]+'/Plots'+figureFolder+ '/' + filename +'_Width_Prominence_Histograms_TransRef.pdf', dpi=200,bbox_inches='tight')

            writeFile = open(rootFolder+'/'+parentFolders[x]+'/'+subfolders[a]+'/Plots'+figureFolder+'/'+filename+'_HistogramData.csv', "w",newline='')
            writer = csv.writer(writeFile,delimiter=',')
            writer.writerow(['PeakWidthBinsTrans','PeakWidthCountsTrans','PeakWidthBinsRef','PeakWidthCountsRef','PeakPromBinsTrans','PeakPromCountsTrans','PeakPromBinsRef','PeakPromCountsRef'])
            rows=zip(TransWidthhistbins,TransWidthcounts,RefWidthhistbins,RefWidthcounts,TransAmphistbins,TransAmpcounts,RefAmphistbins,RefAmpcounts)
            for row in rows:
                writer.writerow(row)
            #writer.writerow([TransWidthhistbins,TransWidthcounts,RefWidthhistbins,RefWidthcounts,TransAmphistbins,TransAmpcounts,RefAmphistbins,RefAmpcounts,peakWidthTransAll,peakAmpTransAll,peakWidthRefAll,peakAmpRefAll])
            writeFile.close()

            writeFile = open(rootFolder+'/'+parentFolders[x]+'/'+subfolders[a]+'/Plots'+figureFolder+'/'+filename+'PeakWidths_PeakProminences.csv', "w",newline='')
            writer = csv.writer(writeFile,delimiter=',')
            writer.writerow(['RawWidthsTrans','RawPromTrans','RawWidthsRef','RawPromRef'])
            rows=zip(peakWidthTransAll,peakAmpTransAll,peakWidthRefAll,peakAmpRefAll)
            for row in rows:
                writer.writerow(row)
            #writer.writerow([TransWidthhistbins,TransWidthcounts,RefWidthhistbins,RefWidthcounts,TransAmphistbins,TransAmpcounts,RefAmphistbins,RefAmpcounts,peakWidthTransAll,peakAmpTransAll,peakWidthRefAll,peakAmpRefAll])
            writeFile.close()

            with open(rootFolder+'/'+parentFolders[x]+'/'+subfolders[a]+figureFolder+'_MeanParams.txt', 'w',encoding='UTF8', newline='') as f:  
                f.writelines(['Mean trans Width= '+str(round(np.mean(peakWidthTransAll),5))+' ± '+str(round(np.std(peakWidthTransAll),5))+u'\u00B1'+' Mean trans Amplitude= '+str(round(np.mean(peakAmpTransAll),5))+' ± '+str(round(np.std(peakAmpTransAll),5))+u'\u00B1'+' N_trans_peaks/file= '+str(len(peakAmpTransAll)/len(files))+' Mean ref Width= '+str(round(np.mean(peakWidthRefAll),5))+' ± '+str(round(np.std(peakWidthRefAll),5))+u'\u00B1'+' Mean Ref Amplitude= '+str(round(np.mean(peakAmpRefAll),5))+' ± '+str(round(np.std(peakAmpRefAll),5))+u'\u00B1'+' N_ref_peaks/file= '+str(len(peakAmpRefAll)/len(files))])
    
    plt.close('all')

print("yay we're finished")