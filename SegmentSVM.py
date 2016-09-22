# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 15:05:12 2016

@author: Rogers Lu
"""

## This program is used to classify the boundary voxels using SVM model

import ThesisFunctions as tf
import os
import numpy as np
import nibabel as nib
from sklearn import preprocessing
from sklearn.metrics import accuracy_score
from sklearn.externals import joblib

## Read the built model
clf = joblib.load('modelSVMDif.pkl')

## Read the image file(in a loop)
#dataPath = 'E:\MasterStudy\ThesisProject\Data\thresholded\auto'
#autoImageFile = os.path.join(dataPath, 'auto_34.nii.gz')
autoImageFile = 'auto_42.nii.gz'
scanImageFile = 'scan_42.nii.gz'
autoImage = nib.load(autoImageFile)
scanImage = nib.load(scanImageFile)
segmImageData = autoImage.get_data()
scanImageData = scanImage.get_data()
## copy a segmImageData
newSegmImageData = segmImageData

## Before calling the feature extraction function, boundary voxels need to be deteted
## Use 'Set' to store the index of the perspective target voxels and to avoid repetition
dim1 = len (segmImageData[:,0,0])
dim2 = len (segmImageData[0,:,0])
dim3 = len (segmImageData[0,0,:])
classifyingIndex = set()
neighborVoxels = set()
tupIndex = ()
feature = []

## Extract the those perspevtive classifying voxels first
## After the loop, classifyingIndex should contain all the perspective voxels
for d1 in range(dim1):
    for d2 in range(dim2):
        for d3 in range(dim3):
            #index0 = [d1,d2,d3]
            tupIndex = (d1,d2,d3)            
            if tf.detectBoundary(tupIndex,segmImageData):
                neighborVoxels = tf.returnClassifyingVoxels(tupIndex,segmImageData)
                for voxel in neighborVoxels:
                    if voxel not in classifyingIndex:                
                        classifyingIndex.add(voxel)  
                        
## Get all the true label of those classifying voxels
manualImageFile = 'SF42_mws.nii.gz'
manuaImage = nib.load(manualImageFile)
manualImageData = manuaImage.get_data()
manualLabels = []

## Call the extracting feature funciton
for indexP in classifyingIndex:
    feature.append(tf.ExtractFeatures60(indexP,scanImageData))
    manualLabels.append(manualImageData[indexP[0],indexP[1],indexP[2]])
## Classify this voxel
# First we scale the feature
scaler = preprocessing.StandardScaler().fit(feature)
featureScaled = scaler.transform(feature)
predicted = clf.predict(featureScaled)
predicted = predicted.tolist()
for i in range(0,len(predicted)):
     predicted[i] = int(predicted[i])
accuracy = accuracy_score (manualLabels, predicted)

#predicted = predicted.reshape(1,-1)
#b1 = voxel[0]
#b2 = voxel[1]
#b3 = voxel[2]
#print(predicted[0])
#newSegmImageData[b1,b2,b3] = predicted[0,0]

## Write a new image according to newly classified labels
## Use the classified labels to update the value of the segmentation image

outputImage = autoImage
outputImageData = outputImage.get_data()
count = 0
for indexP in classifyingIndex:
    outputImageData[indexP] = predicted[count]
    count = count + 1

new_image = nib.Nifti1Image(outputImageData, outputImage.affine)
nib.save(new_image, 'newSegm42Dif.nii.gz')
