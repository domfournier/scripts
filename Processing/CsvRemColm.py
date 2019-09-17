import re
import os
import numpy as np
"""
CsvRemColm

Function to remove columns from a csv file.

INPUT
sep: Column seperator
col_1 ... : Columns number to be removed

Writen: Feb 6, 2017
Author: dominiquef@mirageoscience.com

"""

## Input
work_dir = 'C:\Users\shannonf\Desktop\Current Gocad\OD\Drilling\OD_drilling\New folder'
input_file = '9 OD_Sample.csv'
out_file = '9 OD_Sample_EDT.csv'

colm = np.asarray([2, 5, 6, 7, 8, 9, 10, 11, 12, -5, -4, -3, -2, -1, 0])


sep = ','

## SCRIPT STARTS HERE # #
fidin = open(work_dir + os.path.sep + input_file, 'r')
fidout = open(work_dir + os.path.sep + out_file, 'w')

line = fidin.readline()
line = np.array(line.split(sep))

# Change index to python 0 and remove return
colm -= 1
colm = np.r_[colm]
out = np.ones(len(line), dtype=bool)
out[colm] = 0

count = 0
while len(line)>1:
#while count < 2: 
   
    line = line[out]
    
    if not re.match('\n',line[-1]):
        line[-1] = line[-1] + '\n'

    fidout.write(",".join(line))
    line = fidin.readline()
    line = np.array(line.split(sep))
    count+=1

fidout.close()
fidin.close()