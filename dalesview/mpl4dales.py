import logging
import math
import numpy
import matplotlib.pyplot as plt
from dalesdata import dataslice

# Matplotlib plotting backend for dalesview.

log = logging.getLogger(__name__)


class Mpl4Dales(object):

    def __init__(self):
        pass

    # Plotting interface method, entry point of the class
    @staticmethod
    def plot(data, **kwargs):
        if isinstance(data, dataslice.Profile):
            num_times = len(data.times)
            threshold2d = kwargs.get("threshold2d", 10)
            if num_times <= threshold2d:
                return Mpl4Dales.plot_profiles(data, **kwargs)
            else:
                return Mpl4Dales.plot_profiles_2d(data, **kwargs)
        elif isinstance(data, dataslice.TimeSeries):
            return Mpl4Dales.plot_series(data, **kwargs)

    # Plots 1D profiles, vertically. Can handle multiple (up to 10) time steps
    # of a variable profile.
    @staticmethod
    def plot_profiles(profile, **kwargs):
        if len(profile.shape) == 1:
            vals, heights = numpy.array(profile[:]), numpy.array(profile.heights)
            vmask = numpy.not_equal(vals, profile.miss_vals.get(profile.variable, None))
            hmask = numpy.not_equal(heights, profile.miss_vals.get("height", None))
            mask = numpy.logical_and(vmask, hmask)
            plt.plot(vals[mask], heights[mask], color='b', linewidth=2.)
        elif len(profile.shape) == 2:
            num_plots = len(profile.times)
            lines = []
            heights = numpy.array(profile.heights)
            hmask = numpy.not_equal(heights, profile.miss_vals.get("height", None))
            for i in range(num_plots):
                time = profile.times[i]
                if time == profile.miss_vals["time"]:
                    continue
                vals = numpy.array(profile[i, :])
                vmask = numpy.not_equal(vals, profile.miss_vals.get(profile.variable, None))
                mask = numpy.logical_and(vmask, hmask)
                time_label = "time " + str(profile.times[i]) + " s"
                lines.append(plt.plot(vals[mask], profile.heights, color='b', alpha=math.exp(-2 * i / num_plots),
                                      label=time_label, linewidth=2.)[0])
            if num_plots < 10:
                plt.legend(handles=lines)
        else:
            raise Exception("Plotting profiles with shape %s is not supported" % str(profile.shape))
        plt.ylabel("height [m]")
        x_label = profile.variable
        if profile.unit:
            x_label += " [" + profile.unit + "]"
        plt.xlabel(x_label)
        if kwargs.get("show", True):
            plt.show()

    # Plots 2D profile fields.
    @staticmethod
    def plot_profiles_2d(profile, **kwargs):
        contours = plt.contourf(numpy.array(profile.times) / 3600., profile.heights, profile[:, :].transpose(), 20,
                                cmap=plt.cm.jet)
        plt.ylabel("height [m]")
        plt.xlabel("time [h]")
        c_label = profile.variable
        cb = plt.colorbar(contours)
        if profile.unit:
            c_label += " [" + profile.unit + "]"
        cb.set_label(c_label)
        if kwargs.get("show", True):
            plt.show()

    # Plots time series
    @staticmethod
    def plot_series(timeseries, **kwargs):
        if len(timeseries.shape) == 1:
            vals, times = numpy.array(timeseries[:]), numpy.array(timeseries.times)
            vmask = numpy.not_equal(vals, timeseries.miss_vals.get(timeseries.variable, None))
            tmask = numpy.not_equal(times, timeseries.miss_vals.get("time", None))
            mask = numpy.logical_and(vmask, tmask)
            plt.plot(times[mask] / 3600., vals[mask], color='r', linewidth=2)
            plt.xlabel("time [h]")
            y_label = timeseries.variable
            if timeseries.unit:
                y_label += " [" + timeseries.unit + "]"
            plt.ylabel(y_label)
        else:
            raise Exception("Plotting time series with shape %s is not supported" % str(timeseries.shape))
        if kwargs.get("show", True):
            plt.show()

    # Plots multiple data sets in subplots.
    @staticmethod
    def multiplot(data_list, **kwargs):
        num_data = float(len(data_list))
        ncols = int(math.ceil(math.sqrt(3 * num_data / 2)))
        nrows = int(math.ceil(math.sqrt(2 * num_data / 3)))
        k = 0
        for i in range(nrows):
            for j in range(ncols):
                if k < len(data_list):
                    plt.subplot(nrows, ncols, k + 1)
                    kwargs["show"] = False
                    Mpl4Dales.plot(data_list[k], **kwargs)
                k += 1
        plt.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.98, hspace=0.5, wspace=0.4)
        plt.show()
        return None
