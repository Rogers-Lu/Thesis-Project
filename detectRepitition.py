# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 16:56:46 2016

@author: Rogers Lu
"""

f_index = open ('Index34.txt', 'r')
indexNew = set()
index = []
linesIndex = f_index.readlines()
for line00 in linesIndex :
    line00 = line00.split()
    index.append(line00)
    tup = (line00[0],line00[1],line00[2])
    indexNew.add(tup)
length1 = len(indexNew)
length2 = len(index)