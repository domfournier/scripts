# -*- coding: utf-8 -*-
"""
Run H3DTD forward models on MIRA cluster

Run through a list of Observation, mesh, model files and call H3DTD on the
cluster. The program than copies and rename the pred0 file to final directory.

Created on Sat Nov 12 10:27:25 2016

@author: dominiquef
"""

import os
import glob
import numpy as np
import re

# # INPUT # #

# # MAIN # #
root_dir = os.getcwd()
inp_dir = root_dir + '/INP_DIR/'
out_dir = root_dir + '/OUT/'


inv_dir = os.listdir(inp_dir)
# Look at input directory and load in all models
for curr_dir in inv_dir:

    if os.path.isdir(inp_dir + curr_dir):

        os.chdir(inp_dir + curr_dir)

        # Look for model files
        mfiles = glob.glob("*.sus")
        if mfiles:
            count = 0
            for m in mfiles:
                iter = int(re.search('\d+', m).group())
                if iter > count:
                    keep = m
                    count = iter

            tID = int(re.search('\d+', curr_dir).group())

            os.system('cp ' + m + ' ' + out_dir + 'Model'+str(tID)+'.sus')
            os.system('cp Tile.msh ' + out_dir + 'Mesh'+str(tID)+'.msh')



    # Move and rename pred0 and delete inversion files
    # os.system('cp recv_h3dtd.txt '+ out_dir + '/FWR_Tile_'+str(tID)+'.dat')
    # os.system('rm recv_h3dtd.txt')
    # os.system('rm h3dtd.log')
    #os.system('rm h3dtdinv.out')
    #os.system('rm h3dtdinv_stat.txt')
