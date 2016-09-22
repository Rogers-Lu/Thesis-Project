# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 20:30:06 2016

@author: Administrator
"""

import nibabel as nib

index = [267,250,105]
offset1 = []
offset2 = []
window = []
scanImageFile = 'scan_34.nii.gz'
scanImage = nib.load(scanImageFile)
scanImageData = scanImage.get_data()

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
f_offset1.close()
f_offset2.close()
f_window.close()
    