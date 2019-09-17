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
import sys

# # INPUT # #
input_dir = 'INP_DIR/'
out_dir = 'OUT_DIR'
nodes = ['g1', 'g2']
ppn = ['12', '12']

nQ = 2

# # MAIN # #
os.chdir(input_dir)
inv_dir = os.listdir('.')
print(inv_dir)
count = 0
# Loop through the tiles and Q a job
for ii in inv_dir:

    count += 1
    # Cycle through the nodes
    idn = np.remainder(count, len(nodes))

    if not os.path.isdir(ii):
        continue

    os.chdir(ii)

    fid = open('inv.inp', 'w')

    fid.write('0   ! irest\n')
    fid.write('1  ! mode\n')
    fid.write('1  0.025   ! par tolc\n')
    fid.write('Tile_data.dat ! observations file\n')
    fid.write('maginv3d.mtx\n')
    fid.write('null  ! initial model\n')
    fid.write('null  ! reference model\n')
    fid.write('null  !active cell file\n')
    fid.write('0 1  ! lower, upper bounds\n')
    fid.write('4.5e-5 1 1 0.0625 ! Le, Ln, Lz\n')
    fid.write('SMOOTH_MOD_DIF\n')
    fid.write('null  ! weighting file\n')
    fid.write('0\n')
    fid.close()

    # Write a pbs file and Q
    fid = open('Tile' + str(count) + '.pbs', 'w')
    fid.write('#PBS -l nodes=' + nodes[idn] + ':ppn=' + ppn[idn] + ':gold\n')
    fid.write('cd $PBS_O_WORKDIR\n')
    fid.write('magsen3d sens.inp\n')
    fid.write('maginv3d_newbounds inv.inp\n')
    fid.close()

    os.system('qsub ' + 'Tile' + str(count) + '.pbs')

    os.chdir('../')

    # Move and rename pred0 and delete inversion files
    # os.system('cp recv_h3dtd.txt '+ out_dir + '/FWR_Tile_'+str(tID)+'.dat')
    # os.system('rm recv_h3dtd.txt')
    # os.system('rm h3dtd.log')
    #os.system('rm h3dtdinv.out')
    #os.system('rm h3dtdinv_stat.txt')
