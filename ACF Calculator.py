# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 12:19:34 2023

@author: Daniel

This code takes an input of a CSV including time and intensity columns, and outputs a CSV of the autocorrelation function
"""

import statsmodels.api as sm
import csv


#File paths are predefined, and must manually be switched between reflection and transmission, and the four proteins
getFolder=r"D:\Stitched Files\Raw Stitch"
getFile='\\Stitched Data Rf.csv'
#getFile='\\Stitched Data Tr.csv'
saveFolder=r"D:\Raw Stitched\Stitched Data"
#subfolder='\\Aprotinin'
#subfolder='\\CarbAn'
#subfolder='\\cMyc'
subfolder='\\Strep'
Results='\\Results'
Raw='\\Raw'
ACF='\\ACF'

#The data are imported as time and intensity
time=[]
intens=[]
file = open(getFolder+subfolder+getFile) 
csvreader = csv.reader(file)
for row in csvreader:
    try:
        time.append(float(row[0]))
        intens.append(float(row[1]))
    except IndexError:
        pass
file.close()

#The autocorrelation function, g, is generated and shorted to only 1000000 points. 
g = sm.tsa.acf(intens,nlags=len(time))
gtemp=[]
for j in range(0,1000000):
    gtemp.append(g[j])
g=gtemp
#A list of autocorrelation lag times is generated. 
gtime=[]
for j in range(0,1000000):
    gtime.append(time[j]-time[0])

#The autocorrelation data are stored
writeFile = open(saveFolder+subfolder+Results+ACF+'\ACF Rf.csv', "w",newline='')
#writeFile = open(saveFolder+subfolder+Results+ACF+'\ACF Tr.csv', "w",newline='')
writer = csv.writer(writeFile,delimiter=',')
writer.writerow(['G Time','G'])
for n in range(0,len(gtime)):
    writer.writerow([gtime[n],g[n]])
writeFile.close()