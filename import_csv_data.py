#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 29 19:54:58 2020

@author: lisa-marianeedham
"""
import os
import numpy as np

def fn_import_csv(location):
    
    """
    Imports .txt files and creates a list. 
    Inputs: file location as path
    
    """
    import os
    import numpy as np 
    files=[]
    for file in os.listdir(location):
        if file.endswith(".txt"):
            files.append(file)
        elif file.endswith(".csv"):
            files.append(file)
    for i in range(len(files)):
        files[i]=os.path.join(location, files[i])

    return files


def fn_load_array(file_list):
    
    """
   Loads csv data into array
    
    """

    
    for i in range(len(file_list)):
        
        data = np.genfromtxt(file_list[i], delimiter='\t', dtype=np.str, skip_header=1)
        
        return
    
def fn_write_text_file(folder,filename,subData,data):
     
    path=folder
    completeName = os.path.join(path, filename+".txt")                  
    with open(completeName,mode="w+") as file:
        for subData in data:
            file.write(str(linewidth)+"\n") 