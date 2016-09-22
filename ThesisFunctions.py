# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 19:44:43 2016

@author: Rogers Lu
"""

##This program is used to extract the features

def ExtractFeatures(index,scanImageData):
    ## define some necessary variables
    #  define the offset 
    #  read the offset settings and window settings from those two files
    offset = []
    window = []
    f_offset = open ('offsetSettings.txt', 'r')
    f_window = open ('windowSettings.txt', 'r') 
    linesOffset = f_offset.readlines()
    linesWindow = f_window.readlines()    
    window = linesWindow[0].split()
    for line0 in linesOffset :
        line0 = line0.split()
        offset.append(line0)

    ## Change the list of string into list of integer
    rows1 = len(offset)
    cols1 = len(offset[0]) 
    for row1 in range(rows1):
        for col1 in range(cols1):
            offset[row1][col1] = int(offset[row1][col1])
            
    for i in range(0,len(window)):
        window[i]=int(window[i])

    ## Calculate the offset index
    index_offset = []
    for offsetTemp in offset:
        index_offset.append([a+b for a,b in zip(index,offsetTemp)])

    ## Read the scan image
    #scanImageFile = scanImagePath
    #scanImage = nib.load(scanImageFile)
    #scanImageData = scanImage.get_data()

    ##Loop through the window area and calculate the mean intensity value within this window
    #instaniate a feature list
    feature = []
    patchSize = window[0] * window[1] * window[2]
    for indexT in index_offset:
        s = 0
        for k in range(window[0]):
            for m in range(window[1]):
                for n in range(window[2]):
                    a = indexT[0]+k
                    b = indexT[1]+m
                    c = indexT[2]+n
                    s = s + scanImageData[a,b,c]
        feature.append(s/patchSize)
    
    return feature
    f_offset.close()
    f_window.close()
    

## This function is used to extract 60 features containing the difference between two windows

def ExtractFeatures60(index,scanImageData):
    ## define some necessary variables
    #  define the offset 
    #  read the offset settings and window settings from those two files
    offset1 = []
    offset2 = []
    window = []
    f_offset1 = open ('offsetSettings.txt', 'r')
    f_offset2 = open ('offsetSettings2.txt', 'r')    
    f_window = open ('windowSettings.txt', 'r') 
    linesOffset1 = f_offset1.readlines()
    linesOffset2 = f_offset2.readlines()
    linesWindow = f_window.readlines()    
    window = linesWindow[0].split()
    for line1 in linesOffset1 :
        line1 = line1.split()
        offset1.append(line1)
        
    for line2 in linesOffset2 :
        line2 = line2.split()
        offset2.append(line2)

    ## Change the list of string into list of integer
    rows1 = len(offset1)
    cols1 = len(offset1[0]) 
    rows2 = len(offset2)
    cols2 = len(offset2[0]) 
    
    for row1 in range(rows1):
        for col1 in range(cols1):
            offset1[row1][col1] = int(offset1[row1][col1])
   
    for row2 in range(rows2):
        for col2 in range(cols2):
            offset2[row2][col2] = int(offset2[row2][col2])
            
    for i in range(0,len(window)):
        window[i]=int(window[i])

    ## Calculate the offset index
    index_offset1 = []
    for offsetTemp1 in offset1:
        index_offset1.append([a+b for a,b in zip(index,offsetTemp1)])
    
    index_offset2 = []
    for offsetTemp2 in offset2:
        index_offset2.append([a2+b2 for a2,b2 in zip(index,offsetTemp2)])



    ## Read the scan image
    #scanImageFile = scanImagePath
    #scanImage = nib.load(scanImageFile)
    #scanImageData = scanImage.get_data()

    ##Loop through the window area and calculate the mean intensity value within this window
    #instaniate a feature list
    indexT1 = []
    indexT2 = []
    num_offset = len(index_offset1)
    feature = []
    patchSize = window[0] * window[1] * window[2]
    for num in range(num_offset):
        indexT1 = index_offset1[num]
        indexT2 = index_offset2[num]
        s1 = 0
        s2 = 0
        for k in range(window[0]):
            for m in range(window[1]):
                for n in range(window[2]):
                    a1 = indexT1[0]+k
                    b1 = indexT1[1]+m
                    c1 = indexT1[2]+n
                    a2 = indexT2[0]+k
                    b2 = indexT2[1]+m
                    c2 = indexT2[2]+n
                    s1 = s1 + scanImageData[a1,b1,c1]
                    s2 = s2 + scanImageData[a2,b2,c2]
        feature.append(s1/patchSize)
        feature.append((s1-s2)/patchSize)
    
    return feature
    f_offset1.close()
    f_offset2.close()
    f_window.close()
    
    
## This module is used to determine if the voxel is a boundary voxel
## If yes, then return true
def detectBoundary(index,segmImageData):
    label = 1
    if segmImageData[index[0],index[1],index[2]] == 1:
        for i in range(-1,1):
            for j in range(-1,1):
                for k in range(-1,1):
                    p1 = index[0] + i
                    p2 = index[1] + j
                    p3 = index[2] + k
                    label = label * segmImageData[p1,p2,p3]
                    if label == 0:
                        return True
                        break
                break
            break    
    if label == 1:
        return False
            
## This module is used to return all the classifying voxels
## The input voxel is boundary voxel
## The input 'index' is a tuple since it needs to be immutable
## The returned voxels should only be negative neighborhood voxels and the center voxel
def returnClassifyingVoxels(index,segmImageData):
    neighborVoxels = set()
    neighborVoxels.add(index)
    tup = ()
    for i in range(-1,1):
        for j in range(-1,1):
            for k in range(-1,1):
                p1 = index[0] + i
                p2 = index[1] + j
                p3 = index[2] + k
                tup = (p1,p2,p3)
                if segmImageData[p1,p2,p3] ==0:
                     neighborVoxels.add(tup) 
    return neighborVoxels                
    
def detectRepetition():
    f_index = open ('Index34.txt', 'r')
    indexNew = set()
    linesIndex = f_index.readlines()
    for line0 in linesIndex :
        line0 = line0.split()
        tup = (line0)
        indexNew.add(tup)
    return len(indexNew)