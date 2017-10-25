import logging
import math
import numpy
import matplotlib.pyplot as plt
from dalesdata import dataslice

# Matplotlib plotting backend for dalesview.

log = logging.getLogger(__name__)

class mpl4dales(object):

    def __init__(self):
        pass

    def get_figkwargs(self, **kwargs):
        figkwargs = {}
        if "figsize" in kwargs: figkwargs['figsize'] = kwargs['figsize']
        if "dpi" in kwargs: figkwargs['dpi'] = kwargs['dpi']
        if "facecolor" in kwargs: figkwargs['facecolor'] = kwargs['facecolor']
        if "edgecolor" in kwargs: figkwargs['edgecolor'] = kwargs['edgecolor']
        if "frameon" in kwargs: figkwargs['frameon'] = kwargs['frameon']
        if "subplotpars" in kwargs: figkwargs['subplotpars'] = kwargs['subplotpars']
        if "tight_layout" in kwargs: figkwargs['tight_layout'] = kwargs['tight_layout']
        if len(figkwargs)==0:
            figkwargs = None
        return figkwargs

    # Plotting interface method, entrypoint of the class
    def plot(self,data,**kwargs):
        if "show" not in kwargs: 
            kwargs["show"] = True
        if(isinstance(data,dataslice.profile)):
            numtimes = len(data.times)
            threshold2d = kwargs.get("threshold2d",10)
            if(numtimes <= threshold2d):
                return self.plot_profiles(data,**kwargs)
            else:
                return self.plot_profiles_2d(data,**kwargs)
        elif(isinstance(data,dataslice.timeseries)):
            return self.plot_series(data,**kwargs)


    # Plots 1D profiles, vertically. Can handle multiple (up to 10) timesteps
    # of a variable profile.
    def plot_profiles(self,profile,**kwargs):
        if(len(profile.shape) == 1):
            vals,heights = numpy.array(profile[:]),numpy.array(profile.heights)
            vmask = numpy.not_equal(vals,profile.missvals.get(profile.variable,None))
            hmask = numpy.not_equal(heights,profile.missvals.get("height",None))
            mask = numpy.logical_and(vmask,hmask)
            if 'sub_col' in kwargs and 'axs' in kwargs:
                axs = kwargs['axs']
                sub_col = kwargs['sub_col']
                axs[sub_col].plot(vals[mask],heights[mask],color = 'b',linewidth = 2.)
            else:
                plt.plot(vals[mask],heights[mask],color = 'b',linewidth = 2.)
        elif(len(profile.shape) == 2):
            numplots = len(profile.times)
            f = 1. / (numplots)
            lines = []
            heights = numpy.array(profile.heights)
            hmask = numpy.not_equal(heights,profile.missvals.get("height",None))
            for i in range(numplots):
                time = profile.times[i]
                if(time == profile.missvals["time"]): continue
                vals = numpy.array(profile[i,:])
                vmask = numpy.not_equal(vals,profile.missvals.get(profile.variable,None))
                mask = numpy.logical_and(vmask,hmask)
                timelabel = "time " + str(profile.times[i]) + " s"
                if 'sub_col' in kwargs and 'axs' in kwargs:
                    axs = kwargs['axs']
                    sub_col = kwargs['sub_col']
                    lines.append(axs[sub_col].plot(vals[vmask],profile.heights,color = 'b',alpha = math.exp(-2*i/numplots),label = timelabel,linewidth = 2.)[0])
                else:
                    lines.append(plt.plot(vals[vmask],profile.heights,color = 'b',alpha = math.exp(-2*i/numplots),label = timelabel,linewidth = 2.)[0])
            if(numplots < 10):
                plt.legend(handles = lines)
        else:
            raise Exception("Plotting profiles with shape %s is not supported" % str(profile.shape))
        plt.ylabel("height [m]")
        xlabel = profile.variable
        if(profile.unit):
            xlabel += " [" + profile.unit + "]"
        plt.xlabel(xlabel)
        figkwargs = self.get_figkwargs(**kwargs)
        if figkwargs is not None: plt.figure(**figkwargs)
        if(kwargs["show"]==True): plt.show()


    # Plots 2D profile fields.
    def plot_profiles_2d(self,profile,**kwargs):
        #print("Creating 2D-plot...")
        if 'sub_col' in kwargs and 'axs' in kwargs:
            axs = kwargs['axs']
            sub_col = kwargs['sub_col']
            contours = axs[sub_col].contourf(numpy.array(profile.times)/3600.,profile.heights,profile[:,:].transpose(),20,cmap = plt.cm.jet)
        else:
            contours = plt.contourf(numpy.array(profile.times)/3600.,profile.heights,profile[:,:].transpose(),20,cmap = plt.cm.jet)
        plt.ylabel("height [m]")
        plt.xlabel("time [h]")
        clabel = profile.variable
        cb = plt.colorbar(contours)
        if(profile.unit):
            clabel += " [" + profile.unit + "]"
        cb.set_label(clabel)
        figkwargs = self.get_figkwargs(**kwargs)
        if figkwargs is not None: plt.figure(**figkwargs)
        if(kwargs["show"]==True): plt.show()


    # Plots time series
    def plot_series(self,timeseries,**kwargs):
        #print("Creating timeseries-plot")
        if(len(timeseries.shape) == 1):
            vals,times = numpy.array(timeseries[:]),numpy.array(timeseries.times)
            vmask = numpy.not_equal(vals,timeseries.missvals.get(timeseries.variable,None))
            tmask = numpy.not_equal(times,timeseries.missvals.get("time",None))
            mask = numpy.logical_and(vmask,tmask)
            if 'sub_col' in kwargs and 'axs' in kwargs:
                axs = kwargs['axs']
                sub_col = kwargs['sub_col']
                axs[sub_col].plot(times[mask]/3600.,vals[mask],color = 'r',linewidth = 2)
            else:
                plt.plot(times[mask]/3600.,vals[mask],color = 'r',linewidth = 2)
            plt.xlabel("time [h]")
            ylabel = timeseries.variable
            if(timeseries.unit):
                ylabel += " [" + timeseries.unit + "]"
            plt.ylabel(ylabel)
        else:
            raise Exception("Plotting time series with shape %s is not supported" % str(timeseries.shape))
        figkwargs = self.get_figkwargs(**kwargs)
        if figkwargs is not None: plt.figure(**figkwargs)       
        if(kwargs["show"]==True): plt.show()


    # Plots multiple data sets in subplots.
    def multiplot(self,datalist,**kwargs):
        numdata = float(len(datalist))
        ncols = int(math.ceil(math.sqrt(numdata/2)))
        nrows = math.ceil(numdata/float(ncols))#int(math.ceil(math.sqrt(numdata/3)))
        k = 0
        print("number of cols {}, number of rows {}, total of {} plots".format(ncols, nrows, int(numdata)))
        
        figkwargs = self.get_figkwargs(**kwargs)
        #if figkwargs is not None: 
        #fig = plt.figure(**figkwargs)
        #else:
        fig, axs = plt.subplots(nrows, ncols, **figkwargs)
        #fig = plt.figure() 
        for row in axs:
            for col in range(1, ncols):    
                if (k < len(datalist)): 
                    kwargs['sub_col'] = col
                    kwargs['axs'] = row
                    kwargs["show"] = True
                    self.plot(datalist[k], **kwargs)
                    k = k + 1
        #for i in range(nrows):
        #    for j in range(ncols):
        #        if(k < len(datalist)):                    
        #            plt.subplot(nrows, ncols, k + 1)
        #            #ax = fig.add_subplot(nrows, ncols, k+1)
        #            kwargs["show"]= False
        #            self.plot(datalist[k], **kwargs)
        #        k += 1
        plt.subplots_adjust(left = 0.05,right = 0.95,bottom = 0.02,top = 0.98,hspace = 0.4,wspace = 0.5)
        plt.show()
        return None
