#!/usr/bin/env python                                            
from __future__ import division
from __future__ import print_function

import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import json
import sys
import os

# moviepy  -  pip install moviepy
# had to hack anaconda2/envs/clouds/lib/python2.7/site-packages/imageio/plugins/pillowmulti.py
# to handle 'wu' quantizer

from netCDF4 import Dataset # pip install netCDF4

from moviepy.video.io.bindings import mplfig_to_npimage
import moviepy.editor as mpy
exp_dir = 'long3'
if len(sys.argv) > 1:
    exp_dir = sys.argv[1]

in_file = exp_dir+'/spifs.nc'

params = {'legend.fontsize': 'small',
          'axes.labelsize': 'small',
          'axes.titlesize':'small',
          'xtick.labelsize':'7',
          'ytick.labelsize':'7'}
plt.rcParams.update(params)

    
# create output directory
#try:
#    os.mkdir(out_dir)
#except:
#    pass


rootgrp = Dataset(in_file, "r")

def makeFig():
    #fig, ax = plt.subplots(1,figsize=(5,3), facecolor='white')
    
    fig, (sf, sf2, sf3) = plt.subplots(1, 3, sharey=True, facecolor='white', figsize=(8,6))

    sf.set_ylim([0, 3000])
    sf.set_xlim([240, 300])
    l1o, = sf.plot([0], [0], 'o-',  color='blue', mec='none', label='OpenIFS')
    l1d, = sf.plot([0], [0], '--', color='red', label='Dales')
    l1d_2, = sf.plot([0], [0], ':', color='green', label='Dales own', lw=4)
    sf.set_ylabel('Height (m)')
    sf.set_xlabel('T (K)')

    
    sf2.set_xlim([0, 1e-2])
    l2o, = sf2.plot([0], [0], 'o-',  color='blue', mec='none')
    l2d, = sf2.plot([0], [0], '--', color='red')
    sf2.set_xlabel('Qt')
    
    
    sf3.set_xlim([240, 305])
    l3o, = sf3.plot([0],[0], '-',  color='blue')
    l3d, = sf3.plot([0], [0], '--', color='red')
    sf3.set_xlabel('THL (K)')

    
    time = 't ' 
    plt.suptitle(time)

    return fig, l1o, l1d, l1d_2, l2o, l2d, l3o, l3d, 

fps = 24
g = '385'


def makeFrame(t):
    global fig, g
    i = round(fps*t)
  #  print(i)
    
    #l1.set_ydata( zz(2*np.pi*t/duration))  # <= Update the curve

    # Temperature
    l1o.set_xdata(rootgrp[g+'/T'][i][:])
    l1o.set_ydata(rootgrp[g+'/Zf'][i][:])
    l1d.set_xdata(rootgrp[g+'/t'][i][:])
    #print(rootgrp['/zf'][:])
    l1d.set_ydata(rootgrp['/zf'][:])

    l1d_2.set_xdata(rootgrp[g+'/t_'][i][:])
    l1d_2.set_ydata(rootgrp['/zf'][:])
    
    # q_t
    l2o.set_xdata(rootgrp[g+'/SH'][i][:])
    l2o.set_ydata(rootgrp[g+'/Zf'][i][:])
    l2d.set_xdata(rootgrp[g+'/qt'][i][:])
    l2d.set_ydata(rootgrp['/zf'][:])

    # Theta_l
    l3o.set_xdata(rootgrp[g+'/THL'][i][:])
    l3o.set_ydata(rootgrp[g+'/Zf'][i][:])   #note: interpolated - uses dales' height levels
    l3d.set_xdata(rootgrp[g+'/thl'][i][:])
    l3d.set_ydata(rootgrp['/zf'][:])
    time = 't = %0.0f s'%rootgrp['Time'][i] + ' ps = %0.0f'%rootgrp[g+'/Pf'][i][:] [-1]
    fig.suptitle(time)
    
    return mplfig_to_npimage(fig) # RGB image of the figure

fig, l1o, l1d, l1d_2, l2o, l2d, l3o, l3d = makeFig()


nframes=len(rootgrp['Time'])-1

for g in rootgrp.groups:
    animation = mpy.VideoClip(makeFrame, duration=nframes/fps)
    animation.write_videofile(exp_dir + "/profile-%s.mp4"%g, fps=fps)

    


