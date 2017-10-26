from netCDF4 import Dataset
from netCDF4 import MFDataset
import _pickle as pickle
from matplotlib import pyplot
from matplotlib import contour
import numpy
import math
import pandas

'''Script to generate diagnostic comparison plots using measurement data as a reference. '''

## DONE select max height index for OIFS
## DONE select height indices belonging to pre-specified heights
## select timelevels for comparisons
## create movie of profiles: QT, QL, U, V


def n_plots(ax, Xaxis, Ydata,xlab, ylab, labels):
    for i in range(0, len(Xaxis)):
        ax.plot(Xaxis[i], Ydata[i], label = labels[i])
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    lgd, ax.legend(loc=0, bbox_to_anchor=(1.1, 1), prop={'size': 14})
    return ax, lgd
    
    

def two_y_scales(ax1, xaxis1, xaxis2, data1, data2, c1, c2, ylab1, ylab2, xlab, labels):
    """

    Parameters
    ----------
    ax : axis
        Axis to put two scales on

    x-axis1 : array-like
        x-axis values for data set 1 
    
    x-axis2 : array-like
        x-axis values for data set 2

    data1: array-like
        Data for left hand scale

    data2 : array-like
        Data for right hand scale

    c1 : color
        Color for line 1

    c2 : color
        Color for line 2

    Returns
    -------
    ax : axis
        Original axis
    ax2 : axis
        New twin axis
    """
    ax2 = ax1.twinx()

    plt1, = ax1.plot(xaxis1, data1, '-+', color=c1, label = labels[0], )
    ax1.set_xlabel(xlab)
    ax1.set_ylabel(ylab1)
    plt2, = ax2.plot(xaxis2, data2, '-+', color=c2, label = labels[1])
    ax2.set_ylabel(ylab2)
    
    lns = [plt1, plt2]
    labs = [l.get_label() for l in lns]
    ax1.grid('on')
    lgd = ax1.legend(lns, labs, bbox_to_anchor=(1.1, 1), prop={'size': 13})
    
    return ax1, ax2, lgd

def two_x_scales(ax1, xaxis1, xaxis2, data1, data2, c1, c2, xlab1, xlab2, ylab, labels):
    """

    Parameters
    ----------
    ax : axis
        Axis to put two scales on

    x-axis1 : array-like
        x-axis values for data set 1 
    
    x-axis2 : array-like
        x-axis values for data set 2

    data1: array-like
        Data for left hand scale

    data2 : array-like
        Data for right hand scale

    c1 : color
        Color for line 1

    c2 : color
        Color for line 2

    Returns
    -------
    ax : axis
        Original axis
    ax2 : axis
        New twin axis
    """
    ax2 = ax1.twiny()

    plt1, = x1.plot(xaxis1, data1, '-+', color=c1, label = labels[0], )
    ax1.set_ylabel(ylab)
    ax1.set_xlabel(xlab1)
    plt2, = ax2.plot(xaxis2, data2, '-+', color=c2, label = labels[1])
    ax2.set_xlabel(xlab2)
    
    lns = [ax1, ax2]
    labs = [l.get_label() for l in lns]
    ax1.grid('on')
    lgd = ax1.legend(lns, labs, bbox_to_anchor=(1.1, 1), prop={'size': 13})    

    return ax1, ax2, lgd

# Fetching the data
## DALES/OIFS
base_location = "C:/Users\Bram van Es/Dropbox/eScience/data_analysis/" # C:/Users\Bram van Es/ /home/bramiozo/
netcdf = Dataset(base_location+"data/cabau_750steps/spifs/spifs_750_cabau.nc", "r") 
ds = netcdf[list(netcdf.groups.keys())[0]]

## MEASUREMENTS
meascdf = Dataset(base_location+"data/cabau_750steps/dales/profiles.001.nc", "r")
ds_meas = meascdf[list()]

measurements_meta = Dataset(base_location+"data/cabau_750steps/meas/cesar_tower_meteo_lb1_t10_v1.1_201204.nc", "r") 
measurements = []
measurements.append({'datestr': '13-04-2012','ds':Dataset(base_location+"data/cabau_750steps/meas/scm_in.RACMO_Cabauw_2012041312.nc", "r")})
measurements.append({'datestr': '14-04-2012','ds':Dataset(base_location+"data/cabau_750steps/meas/scm_in.RACMO_Cabauw_2012041412.nc", "r")})
measurements.append({'datestr': '15-04-2012','ds':Dataset(base_location+"data/cabau_750steps/meas/scm_in.RACMO_Cabauw_2012041512.nc", "r")})
measurements.append({'datestr': '16-04-2012','ds':Dataset(base_location+"data/cabau_750steps/meas/scm_in.RACMO_Cabauw_2012041612.nc", "r")})
measurements.append({'datestr': '17-04-2012','ds':Dataset(base_location+"data/cabau_750steps/meas/scm_in.RACMO_Cabauw_2012041712.nc", "r")})
measurements.append({'datestr': '18-04-2012','ds':Dataset(base_location+"data/cabau_750steps/meas/scm_in.RACMO_Cabauw_2012041812.nc", "r")})
measurements.append({'datestr': '19-04-2012','ds':Dataset(base_location+"data/cabau_750steps/meas/scm_in.RACMO_Cabauw_2012041912.nc", "r")})


measurements_rds = []
measurements_rds.append({'datestr': '13-04-2012','ds':Dataset(base_location+"data/cabau_750steps/meas/rds_DeBilt_20120413.nc", "r")})
measurements_rds.append({'datestr': '14-04-2012','ds':Dataset(base_location+"data/cabau_750steps/meas/rds_DeBilt_20120414.nc", "r")})
measurements_rds.append({'datestr': '15-04-2012','ds':Dataset(base_location+"data/cabau_750steps/meas/rds_DeBilt_20120415.nc", "r")})
measurements_rds.append({'datestr': '16-04-2012','ds':Dataset(base_location+"data/cabau_750steps/meas/rds_DeBilt_20120416.nc", "r")})
measurements_rds.append({'datestr': '17-04-2012','ds':Dataset(base_location+"data/cabau_750steps/meas/rds_DeBilt_20120417.nc", "r")})
measurements_rds.append({'datestr': '18-04-2012','ds':Dataset(base_location+"data/cabau_750steps/meas/rds_DeBilt_20120418.nc", "r")})
measurements_rds.append({'datestr': '19-04-2012','ds':Dataset(base_location+"data/cabau_750steps/meas/rds_DeBilt_20120419.nc", "r")})

time_series = ds['time'][:]/3600/24
Time_series = netcdf['Time'][:]/3600/24
height_series = netcdf['zf'][:]
Height_series = ds['Zf']

meas_heights = [100, 200, 300, 1000, 2000, 3000, 4000] # heights used for comparison
ref_heights = [numpy.argmin(numpy.abs(h-height_series)) for h in meas_heights]
Ref_heights = [numpy.argmin(numpy.abs(h-Height_series[1])) for h in meas_heights]

meas_times = numpy.linspace(0, 7, num=15) # in half days 
ref_times = [numpy.argmin(numpy.abs(t - Time_series)) for t in meas_times]

E = len(Height_series[1])-1
B = numpy.argmin(numpy.abs(numpy.max(height_series)-Height_series[1]))
numSteps = E - B +1
oifs_range = list(numpy.linspace(E, B, numSteps).astype("int"))

diff_heights = numpy.diff(numpy.hstack((0, height_series)))
diff_Heights = []
for i in range(0, len(Time_series)):
    diff_Heights.append(-numpy.diff(numpy.hstack((Height_series[i], 0))))
    

###
profile_agg_OIFS = ['U', 'V', 'THL', 'f_SH', 'f_U', 'f_V', 'QT', 'QL', 'f_T', 'f_A', 
                                'Pf', 'Ph', 'QI', 'Tv', 'T', 'Zf', 'Zh', 'SH']
profile_agg_DALES = ['u', 'v', 'thl', 't_', 't', 'qt', 'ql', 'f_thl', 
                                'f_qt', 'f_u', 'f_v', 'presf']
xy_DALES = ['lwp', 'rwp', 'twp']

## plot tuples
profile_comparison = [('u', 'U'), ('v','V'), ('thl', 'THL'), ('qt','QT'), ('ql','QL'), 
                     ('t','T'), ('f_u', 'f_U'), ('f_v', 'f_V'), ('presf','Pf')]

profile_comparison_DALES_dual= [('f_thl', 'thl'), ('f_qt', 'qt')] # two vertical axes
profile_comparison_OIFS_mono = [('Pf', 'Ph', 'Psurf')]
profile_comparison_OIFS_dual = [('f_A', 'A')] # two vertical axes
##
### generating aggregates
##
agg_time_series = {}
for varname in profile_agg_DALES:
    total_dales = []
    for i in range(1, len(Time_series)):
        total_dales.append({'timestep': i, 'total':sum(ds[varname][i][:])/len(ds[varname][i][:])})
    agg_time_series[varname]= {'ts': pandas.DataFrame(data=total_dales)}

for varname in profile_agg_OIFS:
    total_OIFS = []
    for i in range(1, len(Time_series)):
        total_OIFS.append({'timestep': i, 'total':sum(ds[varname][i][oifs_range])/len(ds[varname][i][oifs_range])})
    agg_time_series[varname]= {'ts': pandas.DataFrame(data=total_OIFS)}

## Custom aggregation for LWP and TWP
####
total_lw = []
total_LW = []
R = 8.3144621
for i in range(1, len(Time_series)):
    total_lw.append({'timestep': i, 'total':sum(ds['ql'][i][:]*(ds['presf'][i][:]/ds['t'][i][:]/R)*diff_heights[:])});
agg_time_series['lwp'] = {'ts': pandas.DataFrame(data=total_lw)}

for i in range(1, len(Time_series)):
    total_LW_2.append({'timestep': i, 'total':sum(ds['QL'][i][oifs_range]*(ds['Pf'][i][oifs_range]/ds['T'][i][oifs_range]/R)*diff_Heights[i][oifs_range])});
agg_time_series['LWP'] = {'ts': pandas.DataFrame(data=total_LW)}

####
total_tw = []
total_TW = []
for i in range(1, len(Time_series)):
    total_tw.append({'timestep': i, 'total':sum(ds['qt'][i][:]*(ds['presf'][i][:]/ds['t'][i][:]/R)*diff_heights[:])});
agg_time_series['twp'] = {'ts': pandas.DataFrame(data=total_tw)}

for i in range(1, len(Time_series)):
    total_TW.append({'timestep': i, 'total':sum(ds['QT'][i][oifs_range]*(ds['Pf'][i][oifs_range]/ds['T'][i][oifs_range]/R)*diff_Heights[i][oifs_range])});
agg_time_series['TWP'] = {'ts': pandas.DataFrame(data=total_TW)}
####
#### 

## PLOTS
# time evolution of aggregate profile values
fig_size = (20,10)
dpi_ = 300
output_folder = base_location+'/output'

for plot_tuple in profile_comparison:
    figure = pyplot.figure(figsize=fig_size)
    pyplot.plot(Time_series[1:], agg_time_series[plot_tuple[0]]['ts']['total'], label = 'DALES')
    pyplot.plot(Time_series[1:], agg_time_series[plot_tuple[1]]['ts']['total'], label = 'OIFS')
    pyplot.xlabel("time(days)")
    pyplot.ylabel(plot_tuple[0])
    pyplot.title('Evolution of profile averages for {} over time'.format(plot_tuple[1]))
    pyplot.legend()
    pyplot.savefig(output_folder+'/DALES_OIFS_{}_total_timeseries.png'.format(plot_tuple[1]), dpi = dpi_ )
    pyplot.show()
    pyplot.close()

# DALES dual plots
for plot_tuple in profile_comparison_DALES_dual:
    fig, ax = pyplot.subplots(figsize=fig_size)
    ax1, ax2, lgd = two_y_scales(ax, Time_series[1:], Time_series[1:], 
                            agg_time_series[plot_tuple[0]]['ts']['total'],
                            agg_time_series[plot_tuple[1]]['ts']['total'], 
                              'r', 'b', plot_tuple[0], plot_tuple[1], 
                            'time(days)', labels=[plot_tuple[0], plot_tuple[1]])
    pyplot.title('Dales evolution of profile averages for {} over time'.format(plot_tuple[1]))
    pyplot.savefig("output/DALES_dual_{}_total_timeseries.png".format(plot_tuple[1]), 
                   dpi = dpi_, bbox_extra_artists=(lgd,), bbox_inches='tight')
    pyplot.show()
    pyplot.close()

# OIFS dual plots
for plot_tuple in profile_comparison_OIFS_dual:
    fig, ax = pyplot.subplots(figsize=fig_size)
    ax1, ax2, lgd = two_y_scales(ax, Time_series[1:], Time_series[1:], 
                            agg_time_series[plot_tuple[0]]['ts']['total'],
                            agg_time_series[plot_tuple[1]]['ts']['total'], 
                              'r', 'b', plot_tuple[0], plot_tuple[1], 
                            'time(days)', labels=[plot_tuple[0], plot_tuple[1]])
    pyplot.title('Dales evolution of profile averages for {} over time'.format(plot_tuple[1]))
    pyplot.savefig("output/OIFS_dual_{}_total_timeseries.png".format(plot_tuple[1]), 
                   dpi = dpi_, bbox_extra_artists=(lgd,), bbox_inches='tight')
    pyplot.show()
    pyplot.close()


# time evolution of profile values at specific heights, two y-axes for different units
# DALES

# DALES 
for var in profile_agg_DALES:
    fig, ax = pyplot.subplots(figsize=fig_size)
    Xdata = []
    Ydata = []
    labels = []
    print(var)
    for k in ref_heights:
        Xdata.append(Time_series[1:])
        Ydata.append(ds2[var][1:,k])
        labels.append(str(height_series[k])+"(m)")
    
    ax1, lgd = n_plots(ax, Xdata, Ydata, "time(days)", var, labels)
    pyplot.title('DALES, Height-based comparison of profile averages over time')
    pyplot.savefig("output/DALES_{}_levels_timeseries.png".format(var), 
                   dpi = dpi_, bbox_extra_artists=(lgd,), bbox_inches='tight')
    pyplot.show()
    pyplot.close()

# OIFS 
for var in profile_agg_OIFS:
    fig, ax = pyplot.subplots(figsize=fig_size)
    Xdata = []
    Ydata = []
    labels = []
    for k in Ref_heights:
        Xdata.append(Time_series[1:])
        Ydata.append(ds2[var][1:,k])
        labels.append(str(Height_series[1, k])+"(m)")
    
    ax1, lgd = n_plots(ax, Xdata, Ydata, "time(days)", var, labels)
    pyplot.title('OIFS, Height-based comparison of profile averages over time')
    pyplot.savefig("output/OIFS_{}_levels_timeseries.png".format(var), 
                   dpi = dpi_, bbox_extra_artists=(lgd,), bbox_inches='tight')
    pyplot.show()
    pyplot.close()


# height profiles at specific times, two x-axes for different units

pyplot.savefig(outfile+'.png', dpi=600, bbox_extra_artists=(lgd,), bbox_inches='tight')
pyplot.savefig(outfile+'.svg')
pyplot.savefig(outfile+'.eps')