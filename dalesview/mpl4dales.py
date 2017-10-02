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

    # Plotting interface method, entrypoint of the class
    def plot(self,data,**kwargs):
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
        if(kwargs.get("show",True)): plt.show()


    # Plots 2D profile fields.
    def plot_profiles_2d(self,profile,**kwargs):
        contours = plt.contourf(numpy.array(profile.times)/3600.,profile.heights,profile[:,:].transpose(),20,cmap = plt.cm.jet)
        plt.ylabel("height [m]")
        plt.xlabel("time [h]")
        clabel = profile.variable
        cb = plt.colorbar(contours)
        if(profile.unit):
            clabel += " [" + profile.unit + "]"
        cb.set_label(clabel)
        if(kwargs.get("show",True)): plt.show()


    # Plots time series
    def plot_series(self,timeseries,**kwargs):
        if(len(timeseries.shape) == 1):
            vals,times = numpy.array(timeseries[:]),numpy.array(timeseries.times)
            vmask = numpy.not_equal(vals,timeseries.missvals.get(timeseries.variable,None))
            tmask = numpy.not_equal(times,timeseries.missvals.get("time",None))
            mask = numpy.logical_and(vmask,tmask)
            plt.plot(times[mask]/3600.,vals[mask],color = 'r',linewidth = 2)
            plt.xlabel("time [h]")
            ylabel = timeseries.variable
            if(timeseries.unit):
                ylabel += " [" + timeseries.unit + "]"
            plt.ylabel(ylabel)
        else:
            raise Exception("Plotting time series with shape %s is not supported" % str(timeseries.shape))
        if(kwargs.get("show",True)): plt.show()


    # Plots multiple data sets in subplots.
    def multiplot(self,datalist,**kwargs):
        numdata = float(len(datalist))
        ncols = int(math.ceil(math.sqrt(3*numdata/2)))
        nrows = int(math.ceil(math.sqrt(2*numdata/3)))
        k = 0
        fig = plt.figure()
        for i in range(nrows):
            for j in range(ncols):
                if(k < len(datalist)):
                    plt.subplot(nrows,ncols,k + 1)
                    kwargs["show"] = False
                    self.plot(datalist[k],**kwargs)
                k += 1
        plt.subplots_adjust(left = 0.02,right = 0.98,bottom = 0.02,top = 0.98,hspace = 0.5,wspace = 0.4)
        plt.show()
        return None
