#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 12:27:11 2022
Code to obtain all the commands to extract the T2/R2 * of our GaussNewton treated datas. You just get the list, you have to put " cd $HOME/brainvisa-5.0.4 " command in the terminal, next copy and pasteall the results of this algorithm
@author: Orfanidis Alexandre
"""
import nibabel as ni
import glob #glob allow to ignore some parts of a file name
import os
import numpy as np
import subprocess

def aimsROIfeatures(n_sub,n_visits,inputPath,outputPath):
    #init = 'cd $HOME/brainvisa-5.0.4' 
    #os.chdir('/home/user7t/brainvisa-5.0.4')
    for sub in range (n_sub):
        inputPath2=inputPath + 'R2star-sub-' + '%03d' % (sub+1) 
        maskpath="/mnt/senior/derivatives/VolBrain_7T/VolBrain_7T_T1w_UNIDEN_eroded2vxl_and_fusion" +'/sub-' + '%03d' % (sub+1) + '/'
        #print(inputPath2) #follow the algorithm easily
        outputPath2=outputPath + 'sub-' + '%03d' % (sub+1) 
        #print(outputPath2)
        for visits in range(n_visits):
                 inputPath3 = inputPath2 + '_ses-7TV' + "%02d" % (visits+1) + '.nii.gz' 
                 maskpath2 = maskpath + 'ses-7TV' + "%02d" % (visits+1) + '/' + 'sub-' + '%03d' % (sub+1)  + '_ses-7TV' + "%02d" % (visits+1) + '_space-T2starmap_desc-VolBrain_crisp_lab_mask.nii.gz'
                 outputPath3=outputPath2 + '_ses-7TV' + "%02d" % (visits+1) + '_R2Quant_lab.csv'
                 #print(inputPath3)
                 #print(maskpath2)
                 #print(outputPath3)
                 if os.path.exists(inputPath3):
                     #AimsRoiFeatures(s=inputPath3, o=outputPath)
                     command = 'AimsRoiFeatures -i ' + maskpath2 + ' -s ' + inputPath3 + ' -f csv' + '  -o ' + outputPath3
                     print(command) #final print result
                     #os.chdir('HOME/brainvisa')
                     #os.system(command)
                    # _T2Quant_crisp_lab.csv
                     #${subz}${sub}_ses-7T${vzero}${visite}_T2Quant_lab.csv
                        #glob.glob(inputPath + '*echo-' +"%02d" % (echo+1) +'_part-mag_MEGRE.nii.gz') #create the path where we'll fetch the file
                         #.save(img, outputPath+'sub-' + '%03d' % (sub+185) + '_ses-7TV'+ "%02d" % (visits+1)+".nii.gz") #save the img        
               
aimsROIfeatures(185,10,"/home/user7t/aorfanidis/senior/AOrfanidis/Gauss-Newton_results/","/home/user7t/aorfanidis/senior/CSV/AimsRoiFeatures_fusion/")

#cd $HOME/brainvisa 
#AimsRoiFeatures -i /mnt/senior/derivatives/VolBrain_7T/VolBrain_7T_T1w_UNIDEN_eroded2vxl_and_fusion/sub-181/ses-7TV05/sub-181_ses-7TV05_space-T2starmap_desc-VolBrain_crisp_lab_mask.nii.gz -s /home/user7t/aorfanidis/senior/Résultats_Gauss-Newton/T2star-sub-181_ses-7TV05.nii.gz -f csv  -o /home/user7t/aorfanidis/senior/CSV/AimsRoiFeatures_fusion/sub-181_ses-7TV05_T2Quant_lab.csv
