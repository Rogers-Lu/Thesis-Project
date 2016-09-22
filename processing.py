# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 14:33:49 2016

@author: Rogers LU
"""
# This is the program of SVM applied in my thesis 
# The main goal of this program is to tune the parameters and get a best model 
# of this particular problem

from sklearn.externals import joblib
#from sklearn.metrics import confusion_matrix
#from sklearn import metrics
from sklearn.metrics import accuracy_score
from sklearn import svm
from sklearn import preprocessing
#from sklearn.cross_validation import train_test_split

trainingFile = open ('trainingD.txt', 'r')
testFile = open ('testD.txt', 'r')

linesTraining = trainingFile.readlines()
linesTest = testFile.readlines()
trainingFeature =[]
trainingLabel =[]
testFeature = []
testLabel = []

for line1 in linesTraining :
    line1 = line1.split(',')
    length = len(line1)
    training = line1[0:(length - 1)]
    labelTraining = line1[-1]  
    trainingFeature.append(training)
    trainingLabel.append(labelTraining)
  
for line2 in linesTest :
    line2 = line2.split(',')
    length = len(line2)
    test = line2[0:(length - 1)]
    labelTest = line2[-1]  
    testFeature.append(test)
    testLabel.append(labelTest)

#Scaling the features
scaler = preprocessing.StandardScaler().fit(trainingFeature)
trainingFeatureScaled = scaler.transform(trainingFeature)
#trainingFeature = preprocessing.scale(trainingFeature)
#testFeature = preprocessing.scale(testFeature)
testFeatureScaled = scaler.transform(testFeature)


#Use exhaustive search to find the optimal parameters
#C = [8, 9, 10, 11, 12]
#gamma = [0.005,0.01,0.02, 0.03,0.04]

## This is the optimal parameters I can get
ca = 8
g = 0.01
#param_grid = [{'C': [1, 10, 100, 1000], 'gamma': [0.1,0.01,0.001, 0.0001,0.00001], 'kernel': ['rbf']}]

accuracy = []
#for ca in C:
#    for g in gamma:
clf = svm.SVC(C = ca,kernel='rbf', gamma = g, cache_size = 3000)
#clf=svm.SVC() # This is the default setting;
clf.fit(trainingFeatureScaled, trainingLabel)
predictedTest = clf.predict(testFeatureScaled)
accuracy.append(accuracy_score (testLabel, predictedTest))

#Save the built model
joblib.dump(clf, 'modelSVMDif.pkl')
trainingFile.close()
testFile.close()