#!/usr/bin/env python                                            
from __future__ import division
from __future__ import print_function

import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib.cm as cm
import json
import sys
import os

# moviepy  -  pip install moviepy
# with current imageio(2.1.1)
#   had to hack anaconda2/envs/clouds/lib/python2.7/site-packages/imageio/plugins/pillowmulti.py
#   to handle 'wu' quantizer
#   gif output is still broken sometimes
# option: use imageio 1.6: pip install imageio=1.6

from moviepy.video.io.bindings import mplfig_to_npimage
import moviepy.editor as mpy
from PIL import Image  # using pillow, PIL fork.

from netCDF4 import Dataset # pip install netCDF4

exp_dir = 'long3'
if len(sys.argv) > 1:
    exp_dir = sys.argv[1]
    

out_dir = exp_dir+'/out/'
in_file = exp_dir+'/spifs.nc'

lwp_dir = exp_dir+'/lwp/'


params = {'legend.fontsize': 'small',
          'axes.labelsize': 'small',
          'axes.titlesize':'small',
          'xtick.labelsize':'7',
          'ytick.labelsize':'7'}
plt.rcParams.update(params)

# create output directories
for d in (out_dir, lwp_dir):
    try:
        os.mkdir(d)
    except:
        pass
    



# ---- Color & image functions ----

# return an array normalized to the interval 0...1
def normalize(q, low=None, high=None):
    if low == None:
        low = np.min(q)
    if high == None:
        high = np.max(q)
    f = 1.0 / (high - low)
    norm = np.clip((q - low) * f, 0.0, .9999)
    return norm

# Apply a gamma transform, returns q**gamma
# If input is in [0...1] output is also in [0...1]
def gamma(q, gamma):
    return q**gamma

def scale(q, scale=1.0, gamma=1.0):
    return numpy.clip((q*scale)**gamma, 0, .99)

# linear grayscale color map
# uses array broadcasting - add a new axis to q, then multiply along that axis by color
def cm_gray(q):
    color = np.ones((3))*255 # white
    #color = np.array((128,255,16))
    return np.uint8(q[:, :,  np.newaxis] * color)

# apply matplotlib colormap: [0...1] -> RGB bytes
# http://stackoverflow.com/questions/10965417/how-to-convert-numpy-array-to-pil-image-applying-matplotlib-colormap

def data_to_PIL_image(a, cmap=cm.gist_earth):
    img = Image.fromarray(cmap(a, bytes=True))
    return img

# -------------------------------------- #




# matplotlib plot of LWP field.
# Save as an image file if the name is given, otherwise show().
def plot_lwp(lwp, fig_name=None):
    fig,ax = plt.subplots(1, 1, facecolor='white', figsize=(8,6))
    ax.imshow(lwp, cmap ='Blues_r')

    if fig_name:
        fig.savefig(fig_name)
    else:
        plt.show()

# generate a pillow RGB image of the LWP field
#
# TODO: sensible handling of resizing with aspect 
#       maybe input a scale factor instead of a pixel size?
def image_lwp(lwp, fig_name=None, size=None):

    lwp = scale(lwp, 2, 1.5)
    img = data_to_PIL_image(lwp, cmap=cm.Blues_r)
    #img.show(command='viewer.py')
    
    if size:
        img = img.resize((size,size),Image.BILINEAR)
    
    if fig_name:        
        img.save(fig_name)

    return img
        
# load lwp from a file, either .npy or .npz
def load_lwp(file_name):
    #print(file_name)
    if file_name[-1] == 'y':
        lwp = numpy.load(file_name)       # when loading .npy
    else:
        lwp = numpy.load(file_name)['lwp'] # when loading .npz
    return lwp   

lwp_max = 0
rootgrp = Dataset(in_file, "r")


def plot_stats():
    pass


 

#define a makeframe funcion, with a specific g (and fps, cmap)
def frame_maker(g, fps=24, cmap=cm.Blues_r, size=None):

    # paints an LWP field as an RGB picture
    # returns an RGB numpy array, with values in the range 0...255
    # usable with moviepy.VideoClip
    #
    # rootgrp, the group name g, and cmap are global variables 
    def makeframe(t):
        i = round(fps*t) # time to index
        #time = rootgrp[g+'/time']
        lwp = rootgrp[g+'/lwp'][i]

        img = cmap(scale(lwp, scale=2, gamma=.5), bytes=True) # returns RGBA data
        img = img[:,:,0:3] # discard alpha
        #return img
        pic = Image.fromarray(img) # convert to PIL image
        if size:
            pic = pic.resize((size,size),Image.BILINEAR)
        img = numpy.array(pic) # convert PIL image to numpy array
        return img
    return makeframe

r = 4 # number of images per row
clips = [] # list of lists - each list is one row
row = []

fps=24
steps = len(rootgrp['/Time']) - 1 #-1 to avoid rounding errors
dur = steps * 900 / 30.0 / fps

for g in rootgrp.groups:
    clip = mpy.VideoClip(frame_maker(g, size=512), duration=dur)
    row.append(clip)
    if len(row) == r: # start a new row
        clips.append(row)
        row = []

if len(row): # an unfinished row?
    clips.append(row)
        
final_clip = mpy.clips_array(clips)
final_clip.write_videofile(exp_dir+'/clouds.mp4', fps=fps)


