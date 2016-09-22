 #!/usr/bin/python
# -*- coding: utf-8 -*-
import os

executable1 = "E://MasterStudy//ThesisProject//ExperimentResults//ExtractFeatures.exe"
executable2 = "E://MasterStudy//ThesisProject//ExperimentResults//TrainClassifier.exe"
executable3 = "E://MasterStudy//ThesisProject//ExperimentResults//Segment.exe"
executable4 = "D://EvaluateSegmentation-master//build//Release//EvaluateSegmentation.exe"
folder_training = "E://MasterStudy//ThesisProject//ExperimentResults//training//"
folder_test = "E://MasterStudy//ThesisProject//ExperimentResults//test//"
folder_model = "E://MasterStudy//ThesisProject//ExperimentResults//model//"
folder_scan = "E://MasterStudy//ThesisProject//Data//data1//scans//"
path = "E://MasterStudy//ThesisProject//ExperimentResults//"
folder_segment = "E://MasterStudy//ThesisProject//ExperimentResults//segmentation//"
folder_auto = "E://MasterStudy//ThesisProject//Data//thresholded//auto//"
folder_manuals = "E://MasterStudy//ThesisProject//Data//thresholded//manual//"
output_filename1 = "E://MasterStudy//ThesisProject//ExperimentResults//out_accuracy.csv"
output_filename2 = "E://MasterStudy//ThesisProject//ExperimentResults//out.csv"


Cases = range(1,2)

output_file1 = open(output_filename1, "w")
output_file1.write("Case,Accuracy\n")
output_file2 = open(output_filename2, "w")
output_file2.write("Case,Dice\n")

#This is used to execute ExtractFeatures.exe
image = path + "inputListofImage.txt"
mask = path + "inputListofMask.txt"
offset = path + "offsetSettings.txt"
window = path + "windowSettings.txt"

#This is used to indicate the path of the generated classification model
model = folder_model + "RandomForestModel.xml"

#This is used to indicate the folder name of the generated training data and test data
training = folder_training + "training" + ".txt"
test = folder_test + "test" + ".txt"

#Extract the features
if os.path.exists(path):
        cmd1 = executable1 + " " + image + " " + mask + " " + offset + " " + window + " " + training + " " + test
        output1 = os.popen(cmd1).read()

#Train the model
if os.path.exists(training):
        cmd2 = executable2 + " " + training + " " + test + " " + model
        output2 = os.popen(cmd2).read()
        accuracy = ""
        for line in output2.split("\n"):
            line_split = line.split()
            if len(line_split) > 2 and line_split[0] == "Correct":
                accuracy = line_split[3]
                break
        print "Accuracy for case", ":", accuracy
        output_file1.write( accuracy + "\n")
output_file1.close()

#Update the mask using the built model and 
#This is used to calculate the dice score
for case in Cases:
        if case == "23":
                caseStr_manual = "23a"
                caseStr_auto = "23a"
        else:
                caseStr_manual = str(case).zfill(2)
                caseStr_auto = str(case)
                
        manual = folder_manuals + "SF" + caseStr_manual + "_mws.nii.gz"
        auto = folder_auto + "auto_" + caseStr_auto + ".nii.gz"
        scan = folder_scan + "scan_" + caseStr_auto + ".nii.gz"
        segmentation = folder_segment + "Segmentation" + caseStr_auto + ".nii.gz"
        if os.path.exists(path):
                cmd3 = executable3 + " " + scan + " " + auto + " " + model + " " + offset + " " + window + " " + segmentation
                output3 = os.popen(cmd3).read()
                cmd4 = executable4 + " " + manual + " " + segmentation + " -use DICE"
                output4 = os.popen(cmd4).read()
                dice = ""
                for line in output4.split("\n"):
                        line_split = line.split()
                        if len(line_split) > 2 and line_split[0] == "DICE":
                                dice = line_split[2]
                                break
                print "DICE for case", caseStr_auto, ":", dice
                output_file2.write(caseStr_auto + "," + dice + "\n")

output_file2.close()


