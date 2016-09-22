#!/usr/bin/python
# -*- coding: utf-8 -*-
import os

executable = "D://EvaluateSegmentation-master//build//Release//EvaluateSegmentation.exe"
folder_manuals = "E://MasterStudy//ThesisProject//Data//thresholded//manual//"
folder_auto = "E://MasterStudy//ThesisProject//Data//thresholded//auto//"
output_filename = "out.csv"

Cases = range(1,51)

output_file = open(output_filename, "w")
output_file.write("Case,DICE\n")

for case in Cases:
    if case == "23":
        caseStr_manual = "23a"
        caseStr_auto = "23a"
    else:
        caseStr_manual = str(case).zfill(2)
        caseStr_auto = str(case)

    manual = folder_manuals + "SF" + caseStr_manual + "_mws.nii.gz"
    auto = folder_auto + "auto_" + caseStr_auto + ".nii.gz"

    if os.path.exists(manual):
        cmd = executable + " " + manual + " " + auto + " -use DICE"
        output = os.popen(cmd).read()
        dice = ""
        for line in output.split("\n"):
            line_split = line.split()
            if len(line_split) > 2 and line_split[0] == "DICE":
                dice = line_split[2]
                break
        print "DICE for case", caseStr_auto, ":", dice
        output_file.write(caseStr_auto + "," + dice + "\n")

output_file.close()
