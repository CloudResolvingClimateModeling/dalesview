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

    plt1, = ax1.plot(xaxis1, data1, '-+', color=c1, label = labels[0])
    ax1.set_ylabel(ylab)
    ax1.set_xlabel(xlab1)
    plt2, = ax2.plot(xaxis2, data2, '-+', color=c2, label = labels[1])
    ax2.set_xlabel(xlab2)
    
    lns = [plt1, plt2]
    labs = [l.get_label() for l in lns]
    ax1.grid('on')
    lgd = ax1.legend(lns, labs, bbox_to_anchor=(1.1, 1), prop={'size': 13})    

    return ax1, ax2, lgd

## SETTABLE
###########################################

# PUT ALL SETTABLES IN A .JSON file, load as dictionary..

# Fetching the data
## DALES/OIFS
base_location = "C:/Users\Bram van Es/Dropbox/eScience/data_analysis/" # C:/Users\Bram van Es/ /home/bramiozo/
#netcdf = Dataset(base_location+"data/cabau_750steps/spifs/spifs_750_cabau.nc", "r") 
netcdf = Dataset(base_location+"data/cabau_fixed_forcings/spifs.nc", "r") 
ds = netcdf[list(netcdf.groups.keys())[0]]

## MEASUREMENTS
measurements_meta = Dataset(base_location+"data/cabau_750steps/meas/cesar_tower_meteo_lb1_t10_v1.1_201204.nc", "r") 
measurements = []
measurements.append({'datestr': '13-04-2012 12:00','ds':Dataset(base_location+"data/cabau_750steps/meas/scm_in.RACMO_Cabauw_2012041312.nc", "r")})
measurements.append({'datestr': '14-04-2012 12:00','ds':Dataset(base_location+"data/cabau_750steps/meas/scm_in.RACMO_Cabauw_2012041412.nc", "r")})
measurements.append({'datestr': '15-04-2012 12:00','ds':Dataset(base_location+"data/cabau_750steps/meas/scm_in.RACMO_Cabauw_2012041512.nc", "r")})
measurements.append({'datestr': '16-04-2012 12:00','ds':Dataset(base_location+"data/cabau_750steps/meas/scm_in.RACMO_Cabauw_2012041612.nc", "r")})
measurements.append({'datestr': '17-04-2012 12:00','ds':Dataset(base_location+"data/cabau_750steps/meas/scm_in.RACMO_Cabauw_2012041712.nc", "r")})
measurements.append({'datestr': '18-04-2012 12:00','ds':Dataset(base_location+"data/cabau_750steps/meas/scm_in.RACMO_Cabauw_2012041812.nc", "r")})
measurements.append({'datestr': '19-04-2012 12:00','ds':Dataset(base_location+"data/cabau_750steps/meas/scm_in.RACMO_Cabauw_2012041912.nc", "r")})

measurements_rds = []
measurements_rds.append({'datestr': '13-04-2012 12:00','ds':Dataset(base_location+"data/cabau_750steps/meas/rds_DeBilt_20120413.nc", "r")})
measurements_rds.append({'datestr': '14-04-2012 12:00','ds':Dataset(base_location+"data/cabau_750steps/meas/rds_DeBilt_20120414.nc", "r")})
measurements_rds.append({'datestr': '15-04-2012 12:00','ds':Dataset(base_location+"data/cabau_750steps/meas/rds_DeBilt_20120415.nc", "r")})
measurements_rds.append({'datestr': '16-04-2012 12:00','ds':Dataset(base_location+"data/cabau_750steps/meas/rds_DeBilt_20120416.nc", "r")})
measurements_rds.append({'datestr': '17-04-2012 12:00','ds':Dataset(base_location+"data/cabau_750steps/meas/rds_DeBilt_20120417.nc", "r")})
measurements_rds.append({'datestr': '18-04-2012 12:00','ds':Dataset(base_location+"data/cabau_750steps/meas/rds_DeBilt_20120418.nc", "r")})
measurements_rds.append({'datestr': '19-04-2012 12:00','ds':Dataset(base_location+"data/cabau_750steps/meas/rds_DeBilt_20120419.nc", "r")})

time_series = ds['time'][:]/3600/24
Time_series = netcdf['Time'][:]/3600/24
height_series = netcdf['zf'][:]
Height_series = ds['Zf']


fig_size = (20,10)
dpi_ = 300
output_folder = base_location+'/output'

meas_heights = [100, 200, 300, 1000, 2000, 3000, 4000] # heights used for comparison
meas_times = numpy.linspace(0, 7, num=8) # in half days 

###
profile_agg_OIFS = ['U', 'V', 'THL', 'f_SH', 'f_U', 'f_V', 'QT', 'QL', 'f_T', 'f_A', 'Pf', 'Ph', 'QI', 'Tv', 'T', 'Zf', 'Zh', 'SH']
profile_agg_DALES = ['u', 'v', 'thl', 't_', 't', 'qt', 'ql', 'f_thl', 'f_qt', 'f_u', 'f_v', 'presf']
xy_DALES = ['lwp', 'rwp', 'twp']

## plot tuples
profile_comparison =  [('u', 'U'), ('v','V'), ('thl', 'THL'), ('qt','QT'), ('ql','QL'), 
                            ('t','T'), ('f_u', 'f_U'), ('f_v', 'f_V'), ('presf','Pf'), 
                            ('lwp','LWP'), ('twp', 'TWP')]
profile_comparison_DALES_mono = ["(t, t_)"]
profile_comparison_DALES_dual = [('f_thl', 'thl'), ('f_qt', 'qt'), ('qt', 'thl')] # two vertical axes
profile_comparison_OIFS_mono  = [('Pf', 'Ph', 'Psurf')]
profile_comparison_OIFS_dual  = [('f_A', 'A'), ('QT', 'THL')] # two vertical axes

profile_height_DALES = [('t', 'qt'), ('qt', 'thl'), ('f_thl', 'f_qt'), ('ql', 'qt')] # 
profile_height_OIFS = [('T', 'QT'), ('QT', 'THL'), ('QL', 'QT'), ('f_SH', 'SH'), ('A', 'SH')] #

profile_height_comparison = [('t', 'T'), ('qt', 'QT'), ('ql', 'QL')]
profile_height_comparison_water = [('lwp', 'LWP'), ('twp', 'TWP')]


####################

meas_times_night = numpy.linspace(0., 7, num=8) # in half days 
meas_times_afternoon = numpy.linspace(0.5, 7, num=8) # in half days
ref_times_night = [numpy.argmin(numpy.abs(t - Time_series)) for t in meas_times_night]
ref_times_afternoon = [numpy.argmin(numpy.abs(t - Time_series)) for t in meas_times_afternoon]

ref_times_night[0]=1 if ref_times_night[0]==0 else ref_times_night[0]
ref_times_afternoon[0]=1 if ref_times_afternoon[0]==0 else ref_times_afternoon[0]



ref_heights = [numpy.argmin(numpy.abs(h-height_series)) for h in meas_heights]
Ref_heights = [numpy.argmin(numpy.abs(h-Height_series[1])) for h in meas_heights]

ref_times = [numpy.argmin(numpy.abs(t - Time_series)) for t in meas_times]

E = len(Height_series[1])-1
B = numpy.argmin(numpy.abs(numpy.max(height_series)-Height_series[1]))
numSteps = E - B +1
oifs_range = list(numpy.linspace(E, B, numSteps).astype("int"))

diff_heights = numpy.diff(numpy.hstack((0, height_series)))
diff_Heights = []
for i in range(0, len(Time_series)):
    diff_Heights.append(-numpy.diff(numpy.hstack((Height_series[i], 0))))
    

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

lwp_profiles = ds['ql'][:, :]*(ds['presf'][:, :]/ds['t_'][:][:]/R)*diff_heights[:]
twp_profiles = ds['qt'][:, :]*(ds['presf'][:, :]/ds['t_'][:][:]/R)*diff_heights[:]
LWP_profiles = ds['QL'][:, :]*(ds['Pf'][:, :]/ds['T'][:][:]/R)*diff_Heights[:][:]
TWP_profiles = ds['QL'][:, :]*(ds['Pf'][:, :]/ds['T'][:][:]/R)*diff_Heights[:][:]
####
#### 

## PLOTS
####################################################################
# time evolution of aggregate profile values #######################
####################################################################


for plot_tuple in profile_comparison:
    figure = pyplot.figure(figsize=fig_size)
    pyplot.plot(Time_series[1:], agg_time_series[plot_tuple[0]]['ts']['total'], label = 'DALES')
    pyplot.plot(Time_series[1:], agg_time_series[plot_tuple[1]]['ts']['total'], label = 'OIFS')
    pyplot.xlabel("time(days)")
    pyplot.ylabel(plot_tuple[0])
    pyplot.title('Evolution of profile averages for {} over time'.format(plot_tuple[1]))
    pyplot.legend()
    pyplot.savefig(output_folder+'/averages_over_time/DALES_OIFS_{}_total_timeseries.png'.format(plot_tuple[1]), dpi = dpi_ )
    pyplot.show()
    pyplot.close()


# DALES mono plots

for plot_tuple in profile_comparison_DALES_mono:
    fig, ax = pyplot.subplots(figsize=fig_size)
    Xdata = []
    Ydata = []
    labels = []
    for var in plot_tuple:
        Xdata.append(Time_series[1:])
        Ydata.append(agg_time_series[var]['ts']['total'])
        labels.append(var)
    
    ax1, lgd = n_plots(ax, Xdata, Ydata, "time(days)", "/".join(labels), labels)
    pyplot.title('DALES, Comparison of profile averages for {} over time'.format(",".join(labels)))
    pyplot.savefig(output_folder+"/averages_over_time/DALES_{}_total_timeseries.png".format(plot_tuple[1]), 
                   dpi = dpi_, bbox_extra_artists=(lgd,), bbox_inches='tight')
    pyplot.show()
    pyplot.close()


# OIFS mono plots
for plot_tuple in profile_comparison_OIFS_mono:
    fig, ax = pyplot.subplots(figsize=fig_size)
    Xdata = []
    Ydata = []
    labels = []
    for var in plot_tuple:
        Xdata.append(Time_series[1:])
        Ydata.append(agg_time_series[var]['ts']['total'])
        labels.append(var)
    
    ax1, lgd = n_plots(ax, Xdata, Ydata, "time(days)", "/".join(labels), labels)
    pyplot.title('OIFS, Comparison of profile averages for {} over time'.format(",".join(labels)))
    pyplot.savefig(output_folder+"/averages_over_time/OIFS_{}_total_timeseries.png".format(plot_tuple[1]), 
                   dpi = dpi_, bbox_extra_artists=(lgd,), bbox_inches='tight')
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
    pyplot.savefig(output_folder+"/averages_over_time/DALES_dual_{}_total_timeseries.png".format(plot_tuple[1]), 
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
    pyplot.savefig(output_folder+"/averages_over_time/OIFS_dual_{}_total_timeseries.png".format(plot_tuple[1]), 
                   dpi = dpi_, bbox_extra_artists=(lgd,), bbox_inches='tight')
    pyplot.show()
    pyplot.close()


#################################################################################################
### time evolution of profile values at specific heights, two y-axes for different units ########
#################################################################################################

# DALES 
for var in profile_agg_DALES:
    fig, ax = pyplot.subplots(figsize=fig_size)
    Xdata = []
    Ydata = []
    labels = []
    print(var)
    for k in ref_heights:
        Xdata.append(Time_series[1:-1])
        Ydata.append(ds[var][1:-1,k])
        labels.append(str(height_series[k])+"(m)")
    
    ax1, lgd = n_plots(ax, Xdata, Ydata, "time(days)", var, labels)
    pyplot.title('DALES, Height-based comparison of profile averages over time')
    pyplot.savefig(output_folder+"/height_levels/DALES_{}_levels_timeseries.png".format(var), 
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
        Xdata.append(Time_series[1:-1])
        Ydata.append(ds[var][1:-1,k])
        labels.append(str(Height_series[1, k])+"(m)")
    
    ax1, lgd = n_plots(ax, Xdata, Ydata, "time(days)", var, labels)
    pyplot.title('OIFS, Height-based comparison of profile averages over time')
    pyplot.savefig(output_folder+"/height_levels/OIFS_{}_levels_timeseries.png".format(var), 
                   dpi = dpi_, bbox_extra_artists=(lgd,), bbox_inches='tight')
    pyplot.show()
    pyplot.close()

#####################################################################################
## height profiles at specific times, two x-axes for different units  ###############
#####################################################################################
# DALES height profiles at reference times
for timelevel in ref_times:
    timestring = 't = %0.2f days'%(Time_series[timelevel]) #(netcdf_2['Time'][timelevel]/3600/24)
    for plot_tuple in profile_height_DALES:
        fig, ax = pyplot.subplots(figsize=fig_size)
        var1 = plot_tuple[0]
        var2 = plot_tuple[1]
        ax1, ax2, lgd = two_x_scales(ax, 
                              ds[var1][timelevel, :], 
                              ds[var2][timelevel, :],
                              height_series,
                              height_series,
                              'r', 'b', 
                              var1, var2,  
                              'height(m)', labels=[var1, var2])
        pyplot.title("DALES, Average {} and {} at {}".format(var1, var2, timestring), y=1.14)
        pyplot.savefig(output_folder+"/height_profiles/DALES_{}-{}_height_profile_{}_days.png".format(var1, var2, Time_series[timelevel]), 
                       dpi = dpi_, bbox_extra_artists=(lgd,), bbox_inches='tight')
        #pyplot.show()
        pyplot.close()

# OIFS height profiles at reference times
for timelevel in ref_times:
    timestring = 't = %0.2f days'%(Time_series[timelevel]) #(netcdf_2['Time'][timelevel]/3600/24)
    for plot_tuple in profile_height_OIFS:
        fig, ax = pyplot.subplots(figsize=fig_size)
        var1 = plot_tuple[0]
        var2 = plot_tuple[1]
        ax1, ax2, lgd = two_x_scales(ax, 
                              ds[var1][timelevel, oifs_range], 
                              ds[var2][timelevel, oifs_range],
                              Height_series[timelevel, oifs_range],
                              Height_series[timelevel, oifs_range],
                              'r', 'b', 
                              var1, var2,  
                              'height(m)', labels=[var1, var2])
        pyplot.title("OIFS, Average {} and {} at {}".format(var1, var2, timestring), y=1.14)
        pyplot.savefig(output_folder+"/height_profiles/OIFS_{}-{}_height_profile_{}_days.png".format(var1, var2, Time_series[timelevel]), 
                       dpi = dpi_, bbox_extra_artists=(lgd,), bbox_inches='tight')
        #pyplot.show()
        pyplot.close()

# DALES OIFS height profiles at reference times
for timelevel in ref_times:
    timestring = 't = %0.2f days'%(Time_series[timelevel]) #(netcdf_2['Time'][timelevel]/3600/24)
    for plot_tuple in profile_height_comparison:
        figure = pyplot.figure(figsize=fig_size)
        pyplot.plot(ds[plot_tuple[0]][timelevel, :], height_series, label = 'DALES')
        pyplot.plot(ds[plot_tuple[1]][timelevel, oifs_range], Height_series[timelevel, oifs_range], label = 'OIFS')
        pyplot.xlabel(plot_tuple[0])
        pyplot.ylabel("height(m)")
        pyplot.title('Height profile of average {} at {}'.format(plot_tuple[1], timestring))
        pyplot.legend()
        pyplot.savefig(output_folder+"/height_profiles/DALES_OIFS_{}_height_profile_{}_days.png".format(plot_tuple[1], 
                                                                                Time_series[timelevel]), dpi = dpi_)
        #pyplot.show()
        pyplot.close()

# DALES OIFS height profiles for water paths at reference times
for timelevel in ref_times:
    timestring = 't = %0.2f days'%(Time_series[timelevel]) #(netcdf_2['Time'][timelevel]/3600/24)
    figure = pyplot.figure(figsize=fig_size)
    pyplot.plot(lwp_profiles.data[timelevel, :], height_series, label = 'DALES')
    pyplot.plot(LWP_profiles.data[timelevel, oifs_range], Height_series[timelevel, oifs_range], label = 'OIFS')
    pyplot.xlabel(plot_tuple[0])
    pyplot.ylabel("height(m)")
    pyplot.title('Height profile of average LWP at {}'.format(timestring))
    pyplot.legend()
    pyplot.savefig(output_folder+"/height_profiles/DALES_OIFS_LWP_height_profile_{}_days.png".format(Time_series[timelevel]), dpi = dpi_)
    #pyplot.show()
    pyplot.close()

# DALES OIFS height profiles for water paths at reference times
for timelevel in ref_times:
    timestring = 't = %0.2f days'%(Time_series[timelevel]) #(netcdf_2['Time'][timelevel]/3600/24)
    figure = pyplot.figure(figsize=fig_size)
    pyplot.plot(twp_profiles[timelevel, :], height_series, label = 'DALES')
    pyplot.plot(TWP_profiles[timelevel, oifs_range], Height_series[timelevel, oifs_range], label = 'OIFS')
    pyplot.xlabel(plot_tuple[0])
    pyplot.ylabel("height(m)")
    pyplot.title('Height profile of average TWP at {}'.format(timestring))
    pyplot.legend()
    pyplot.savefig(output_folder+"/height_profiles/DALES_OIFS_TWP_height_profile_{}_days.png".format(Time_series[timelevel]), dpi = dpi_)
    #pyplot.show()
    pyplot.close()

########################################################################################
### Comparison plots with measurement data, height profiles ############################
########################################################################################


E = len(measurements[0]['ds']['height_f'][0])-1
B = numpy.argmin(numpy.abs(numpy.max(height_series)-measurements[0]['ds']['height_f'][0]))
numSteps = E - B +1
meas_h_range = list(numpy.linspace(E, B, numSteps).astype("int"))
####
B = numpy.argmin(numpy.abs(numpy.max(height_series)-measurements_rds[0]['ds']['height'][0]))
numSteps = B +1
meas_h_range_rds_afternoon = list(numpy.linspace(0, B, numSteps).astype("int"))
####
B = numpy.argmin(numpy.abs(numpy.max(height_series)-measurements_rds[0]['ds']['height'][1]))
numSteps = B +1
meas_h_range_rds_night = list(numpy.linspace(0, B, numSteps).astype("int"))

## T
# rds: T, absolute velocity, theta_l, Q_L,
### NIGHT
############
for idx, meas in enumerate(measurements_rds):
    pyplot.figure(figsize=fig_size)
    pyplot.plot(meas['ds']['T'][0,meas_h_range_rds_night], meas['ds']['height'][0,meas_h_range_rds_night], label = 'measurement')
    pyplot.plot(ds['t_'][ref_times_night[idx],:], height_series, label = 'DALES')
    pyplot.plot(ds['T'][ref_times_night[idx],oifs_range], Height_series[ref_times_night[idx],oifs_range], '-o', label = 'OIFS')
    pyplot.xlabel("T(K)")
    pyplot.ylabel("height(m)")
    pyplot.title("Temperature versus height, {} 00:00AM".format(meas['datestr']))
    pyplot.legend()
    pyplot.savefig(output_folder+"/meas_comparison/MEAS_Comparison_T_height_profile_{}_night.png".format(meas['datestr']), dpi = dpi_)
    #pyplot.show()

# rds: T, absolute velocity, theta_l, normal: cloud fraction, q, ql, qi, u, v
for idx, meas in enumerate(measurements_rds):
    pyplot.figure(figsize=fig_size)
    pyplot.plot(meas['ds']['THETA'][0,meas_h_range_rds_night], meas['ds']['height'][0,meas_h_range_rds_night], label = 'measurement')
    pyplot.plot(ds['thl'][ref_times_night[idx],:], height_series, label = 'DALES')
    pyplot.plot(ds['THL'][ref_times_night[idx],oifs_range], Height_series[ref_times_night[idx],oifs_range], '-o', label = 'OIFS')
    pyplot.xlabel("T(K)")
    pyplot.ylabel("height(m)")
    pyplot.title("Theta_l versus height, {} 00:00AM".format(meas['datestr']))
    pyplot.legend()
    pyplot.savefig(output_folder+"/meas_comparison/MEAS_Comparison_THETA_height_profile_{}_night.png".format(meas['datestr']), dpi = dpi_)
    #pyplot.show()
    
for idx, meas in enumerate(measurements_rds):
    pyplot.figure(figsize=fig_size)
    pyplot.plot(meas['ds']['QV'][0,meas_h_range_rds_night], meas['ds']['height'][0,meas_h_range_rds_night], label = 'measurement')
    pyplot.plot(ds['qt'][ref_times_night[idx],:] - ds['ql'][ref_times_night[idx],:] , height_series, label = 'DALES')
    pyplot.plot(ds['SH'][ref_times_night[idx],oifs_range], Height_series[ref_times_night[idx],oifs_range], '-o', label = 'OIFS')
    pyplot.xlabel("QV(K)")
    pyplot.ylabel("height(m)")
    pyplot.title("QV versus height, {} 00:00AM".format(meas['datestr']))
    pyplot.legend()
    pyplot.savefig(output_folder+"/meas_comparison/MEAS_Comparison_QV_height_profile_{}_night.png".format(meas['datestr']), dpi = dpi_)
    #pyplot.show()
        
    

### MORNING
##############
for idx, meas in enumerate(measurements_rds):
    pyplot.figure(figsize=fig_size)
    pyplot.plot(meas['ds']['T'][1,meas_h_range_rds_afternoon], meas['ds']['height'][1,meas_h_range_rds_afternoon], label = 'measurement')
    pyplot.plot(ds['t_'][ref_times_afternoon[idx],:], height_series, label = 'DALES')
    pyplot.plot(ds['T'][ref_times_afternoon[idx],oifs_range], Height_series[ref_times_afternoon[idx],oifs_range], '-o', label = 'OIFS')
    pyplot.xlabel("T(K)")
    pyplot.ylabel("height(m)")
    pyplot.title("Temperature versus height, {} 12:00AM".format(meas['datestr']))
    pyplot.legend()
    pyplot.savefig(output_folder+"/meas_comparison/MEAS_Comparison_T_height_profile_{}_afternoon.png".format(meas['datestr']), dpi = dpi_)
    #pyplot.show()

# rds: T, absolute velocity, theta_l, normal: cloud fraction, q, ql, qi, u, v
for idx, meas in enumerate(measurements_rds):
    pyplot.figure(figsize=fig_size)
    pyplot.plot(meas['ds']['THETA'][1,meas_h_range_rds_afternoon], meas['ds']['height'][1,meas_h_range_rds_afternoon], label = 'measurement')
    pyplot.plot(ds['thl'][ref_times_afternoon[idx],:], height_series, label = 'DALES')
    pyplot.plot(ds['THL'][ref_times_afternoon[idx],oifs_range], Height_series[ref_times_afternoon[idx],oifs_range], '-o', label = 'OIFS')
    pyplot.xlabel("T(K)")
    pyplot.ylabel("height(m)")
    pyplot.title("Theta_l versus height, {} 12:00AM".format(meas['datestr']))
    pyplot.legend()
    pyplot.savefig(output_folder+"/meas_comparison/MEAS_Comparison_THETA_height_profile_{}_afternoon.png".format(meas['datestr']), dpi = dpi_)
    #pyplot.show()
    
for idx, meas in enumerate(measurements_rds):
    pyplot.figure(figsize=fig_size)
    pyplot.plot(meas['ds']['QV'][1,meas_h_range_rds_afternoon], meas['ds']['height'][1,meas_h_range_rds_afternoon], label = 'measurement')
    pyplot.plot(ds['qt'][ref_times_afternoon[idx],:] - ds['ql'][ref_times_afternoon[idx],:] , height_series, label = 'DALES')
    pyplot.plot(ds['SH'][ref_times_afternoon[idx],oifs_range], Height_series[ref_times_afternoon[idx],oifs_range], '-o', label = 'OIFS')
    pyplot.xlabel("QV(K)")
    pyplot.ylabel("height(m)")
    pyplot.title("QV versus height, {} 12:00AM".format(meas['datestr']))
    pyplot.legend()
    pyplot.savefig(output_folder+"/meas_comparison/MEAS_Comparison_QV_height_profile_{}_afternoon.png".format(meas['datestr']), dpi = dpi_)
    #pyplot.show()


##
## normal: cloud fraction, q, ql, qi, u, v
for idx, meas in enumerate(measurements):
    pyplot.figure(figsize=fig_size)
    pyplot.plot(numpy.sqrt(meas['ds']['u'][0, meas_h_range]**2+meas['ds']['v'][0, meas_h_range]**2), meas['ds']['height_f'][0, meas_h_range], label = 'measurement')
    pyplot.plot(numpy.sqrt(ds['u'][ref_times[idx],:]**2+ds['v'][ref_times[idx],:]**2), height_series, label = 'DALES')
    pyplot.plot(numpy.sqrt(ds['U'][ref_times[idx],oifs_range]**2+ds['V'][ref_times[idx],oifs_range]**2), Height_series[ref_times[idx],oifs_range], '-o', label = 'OIFS')
    pyplot.xlabel("Vabs(m/s)")
    pyplot.ylabel("height(m)")
    pyplot.title("Absolute velocity versus height, {} 12:00AM".format(meas['datestr']))
    pyplot.legend()
    pyplot.savefig(output_folder+"/meas_comparison/MEAS_Comparison_Vabs_height_profile_{}_afternoon.png".format(meas['datestr']), dpi = dpi_)
    #pyplot.show()
    
# rds: U velocity
for idx, meas in enumerate(measurements):
    pyplot.figure(figsize=fig_size)
    pyplot.plot(meas['ds']['u'][0, meas_h_range], meas['ds']['height_f'][0, meas_h_range], label = 'measurement')
    pyplot.plot(ds['u'][ref_times[idx],:], height_series, label = 'DALES')
    pyplot.plot(ds['U'][ref_times[idx],oifs_range], Height_series[ref_times[idx],oifs_range], '-o', label = 'OIFS')
    pyplot.xlabel("U(m/s)")
    pyplot.ylabel("height(m)")
    pyplot.title("Velocity u versus height, {} 12:00AM".format(meas['datestr']))
    pyplot.legend()
    pyplot.savefig(output_folder+"/meas_comparison/MEAS_Comparison_U_height_profile_{}_afternoon.png".format(meas['datestr']), dpi = dpi_)
    #pyplot.show()

# rds: V velocity
for idx, meas in enumerate(measurements):
    pyplot.figure(figsize=fig_size)
    pyplot.plot(meas['ds']['v'][0, meas_h_range], meas['ds']['height_f'][0, meas_h_range], label = 'measurement')
    pyplot.plot(ds['v'][ref_times[idx],:], height_series, label = 'DALES')
    pyplot.plot(ds['V'][ref_times[idx],oifs_range], Height_series[ref_times[idx],oifs_range], '-o', label = 'OIFS')
    pyplot.xlabel("V(m/s)")
    pyplot.ylabel("height(m)")
    pyplot.title("Velocity u versus height, {} 12:00AM".format(meas['datestr']))
    pyplot.legend()
    pyplot.savefig(output_folder+"/meas_comparison/MEAS_Comparison_V_height_profile_{}_afternoon.png".format(meas['datestr']), dpi = dpi_)
    #pyplot.show()

# rds: QL
for idx, meas in enumerate(measurements):
    pyplot.figure(figsize=fig_size)
    pyplot.plot(meas['ds']['ql'][0, meas_h_range], meas['ds']['height_f'][0, meas_h_range], label = 'measurement')
    pyplot.plot(ds['ql'][ref_times[idx],:], height_series, label = 'DALES')
    pyplot.plot(ds['QL'][ref_times[idx],oifs_range], Height_series[ref_times[idx],oifs_range], '-o', label = 'OIFS')
    pyplot.xlabel("QL")
    pyplot.ylabel("height(m)")
    pyplot.title("QL versus height, {} 12:00AM".format(meas['datestr']))
    pyplot.legend()
    pyplot.savefig(output_folder+"/meas_comparison/MEAS_Comparison_QL_height_profile_{}_afternoon.png".format(meas['datestr']), dpi = dpi_)
    #pyplot.show()

# rds: QT
for idx, meas in enumerate(measurements):
    pyplot.figure(figsize=fig_size)
    pyplot.plot(meas['ds']['q'][0, meas_h_range], meas['ds']['height_f'][0, meas_h_range], label = 'measurement')
    pyplot.plot(ds['qt'][ref_times[idx],:], height_series, label = 'DALES')
    pyplot.plot(ds['QT'][ref_times[idx],oifs_range], Height_series[ref_times[idx],oifs_range], '-o', label = 'OIFS')
    pyplot.xlabel("QT")
    pyplot.ylabel("height(m)")
    pyplot.title("QT versus height, {} 12:00AM".format(meas['datestr']))
    pyplot.legend()
    pyplot.savefig(output_folder+"/meas_comparison/MEAS_Comparison_QT_height_profile_{}_afternoon.png".format(meas['datestr']), dpi = dpi_)
    #pyplot.show()


# rds: Cloud fraction
for idx, meas in enumerate(measurements):
    pyplot.figure(figsize=fig_size)
    pyplot.plot(meas['ds']['cloud_fraction'][0, meas_h_range], meas['ds']['height_f'][0, meas_h_range], label = 'measurement')
    #pyplot.plot(ds['qt'][ref_times[idx],:], height_series, label = 'DALES')
    pyplot.plot(ds['A'][ref_times[idx],oifs_range], Height_series[ref_times[idx],oifs_range], '-o', label = 'DALES')
    pyplot.xlabel("A")
    pyplot.ylabel("height(m)")
    pyplot.title("Cloud fraction versus height, {} 12:00AM".format(meas['datestr']))
    pyplot.legend()
    pyplot.savefig(output_folder+"/meas_comparison/MEAS_Comparison_A_height_profile_{}_afternoon.png".format(meas['datestr']), dpi = dpi_)
    #pyplot.show()