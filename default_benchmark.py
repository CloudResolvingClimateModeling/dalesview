from netCDF4 import Dataset
from netCDF4 import MFDataset
import _pickle as pickle
from matplotlib import pyplot
from matplotlib import contour
import numpy
import math
import pandas

'''Script to generate diagnostic comparison plots using measurement data as a reference. '''

## select max height index for OIFS
## select height indices belonging to pre-specified heights
## select timelevels for comparisons
## create movie of profiles: QT, QL, U, V


def n_plots(ax, Xaxis, Ydata,xlab, ylab, labels):
    for i in range(0, len(Xaxis)):
        ax.plot(Xaxis[i], Ydata[i], label = labels[i])
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.legend(loc=0, bbox_to_anchor=(1.25, 1), prop={'size': 14})
    return ax
    
    

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

    ax1.plot(xaxis1, data1, '-+', color=c1, label = labels[0], )
    ax1.set_xlabel(xlab)
    ax1.set_ylabel(ylab1)
    ax2.plot(xaxis2, data2, '-+', color=c2, label = labels[1])
    ax2.set_ylabel(ylab2)
    
    lns = [ax1, ax2]
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc=1, prop={'size': 9})
    
    return ax1, ax2

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

    ax1.plot(xaxis1, data1, '-+', color=c1, label = labels[0], )
    ax1.set_ylabel(ylab)
    ax1.set_xlabel(xlab1)
    ax2.plot(xaxis2, data2, '-+', color=c2, label = labels[1])
    ax2.set_xlabel(xlab2)
    
    lns = [ax1, ax2]
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc=1, prop={'size': 9})
    
    return ax1, ax2


## DALES/OIFS
netcdf = Dataset("C:/Users\Bram van Es/Dropbox/eScience/data_analysis/data/cabau_750steps/spifs/spifs_750_cabau.nc", "r") 
ds = netcdf[list(netcdf.groups.keys())[0]]

## MEASUREMENTS
meascdf = Dataset("C:/Users\Bram van Es/Dropbox/eScience/data_analysis/data/cabau_750steps/dales/profiles.001.nc", "r")
ds_meas = meascdf[list()]

measurements_meta = Dataset("C:/Users\Bram van Es/Dropbox/eScience/data_analysis/data/cabau_750steps/meas/cesar_tower_meteo_lb1_t10_v1.1_201204.nc", "r") 
measurements = []
measurements.append(Dataset("C:/Users\Bram van Es/Dropbox/eScience/data_analysis/data/cabau_750steps/meas/scm_in.RACMO_Cabauw_2012041312.nc", "r"))
measurements.append(Dataset("C:/Users\Bram van Es/Dropbox/eScience/data_analysis/data/cabau_750steps/meas/scm_in.RACMO_Cabauw_2012041412.nc", "r"))
measurements.append(Dataset("C:/Users\Bram van Es/Dropbox/eScience/data_analysis/data/cabau_750steps/meas/scm_in.RACMO_Cabauw_2012041512.nc", "r"))
measurements.append(Dataset("C:/Users\Bram van Es/Dropbox/eScience/data_analysis/data/cabau_750steps/meas/scm_in.RACMO_Cabauw_2012041612.nc", "r"))
measurements.append(Dataset("C:/Users\Bram van Es/Dropbox/eScience/data_analysis/data/cabau_750steps/meas/scm_in.RACMO_Cabauw_2012041712.nc", "r"))
measurements.append(Dataset("C:/Users\Bram van Es/Dropbox/eScience/data_analysis/data/cabau_750steps/meas/scm_in.RACMO_Cabauw_2012041812.nc", "r"))

measurements_rds = []
measurements_rds.append(Dataset("C:/Users\Bram van Es/Dropbox/eScience/data_analysis/data/cabau_750steps/meas/rds_DeBilt_20120413.nc", "r"))
measurements_rds.append(Dataset("C:/Users\Bram van Es/Dropbox/eScience/data_analysis/data/cabau_750steps/meas/rds_DeBilt_20120414.nc", "r"))
measurements_rds.append(Dataset("C:/Users\Bram van Es/Dropbox/eScience/data_analysis/data/cabau_750steps/meas/rds_DeBilt_20120415.nc", "r"))
measurements_rds.append(Dataset("C:/Users\Bram van Es/Dropbox/eScience/data_analysis/data/cabau_750steps/meas/rds_DeBilt_20120416.nc", "r"))
measurements_rds.append(Dataset("C:/Users\Bram van Es/Dropbox/eScience/data_analysis/data/cabau_750steps/meas/rds_DeBilt_20120417.nc", "r"))
measurements_rds.append(Dataset("C:/Users\Bram van Es/Dropbox/eScience/data_analysis/data/cabau_750steps/meas/rds_DeBilt_20120418.nc", "r"))

time_series = ds['time'][:]/3600/24
Time_series = netcdf['Time'][:]/3600/24
height_series = netcdf['zf'][:]
Height_series = ds['Zf']

E = len(Height_series[1])-1
B = 67
numSteps = E - B +1
oifs_range = list(numpy.linspace(E, B, numSteps).astype("int"))

diff_heights = numpy.diff(numpy.hstack((0, height_series)))
diff_Heights = []
for i in range(0, len(Time_series)):
    diff_Heights.append(-numpy.diff(numpy.hstack((Height_series[i], 0))))
    


### generating aggregates

total_t = []
total_T = []
for i in range(1, len(Time_series)):
    total_t.append({'timestep': i, 'total':sum(netcdf['38062/t'][i][:])/len(netcdf['38062/t'][i][:])});
t_total_df = pandas.DataFrame(data=total_t)
for i in range(1, len(Time_series)):
    total_T.append({'timestep': i, 'total':sum(netcdf['38062/T'][i][oifs_range])/len(netcdf['38062/T'][i][oifs_range])});
T_total_df = pandas.DataFrame(data=total_T)
####
total_thl = []
total_THL = []
for i in range(1, len(Time_series)):
    total_thl.append({'timestep': i, 'total':sum(netcdf['38062/thl'][i][:])/len(netcdf['38062/thl'][i][:])});
thl_total_df = pandas.DataFrame(data=total_thl)
for i in range(1, len(Time_series)):
    total_THL.append({'timestep': i, 'total':sum(netcdf['38062/THL'][i][oifs_range])/len(netcdf['38062/THL'][i][oifs_range])});
THL_total_df = pandas.DataFrame(data=total_THL)
####
total_qt = []
total_QT = []
for i in range(1, len(Time_series)):
    total_qt.append({'timestep': i, 'total':sum(netcdf['38062/qt'][i][:])/len(netcdf['38062/qt'][i][:])});
qt_total_df = pandas.DataFrame(data=total_qt)
for i in range(1, len(Time_series)):
    total_QT.append({'timestep': i, 'total':sum(netcdf['38062/QT'][i][oifs_range])/len(netcdf['38062/QT'][i][oifs_range])});
QT_total_df = pandas.DataFrame(data=total_QT)
####
total_ql = []
total_QL = []
for i in range(1, len(Time_series)):
    total_ql.append({'timestep': i, 'total':sum(netcdf['38062/ql'][i][:])});
ql_total_df = pandas.DataFrame(data=total_ql)
for i in range(1, len(Time_series)):
    total_QL.append({'timestep': i, 'total':sum(netcdf['38062/QL'][i][oifs_range])});
QL_total_df = pandas.DataFrame(data=total_QL)
####
total_v = []
total_V = []
for i in range(1, len(Time_series)):
    total_v.append({'timestep': i, 'total':sum(netcdf['38062/v'][i][:])/len(netcdf['38062/v'][i][:])});
v_total_df = pandas.DataFrame(data=total_v)
for i in range(1, len(Time_series)):
    total_V.append({'timestep': i, 'total':sum(netcdf['38062/V'][i][oifs_range])/len(netcdf['38062/V'][i][oifs_range])});
V_total_df = pandas.DataFrame(data=total_V)
#####
total_u = []
total_U = []
for i in range(1, len(Time_series)):
    total_u.append({'timestep': i, 'total':sum(netcdf['38062/u'][i][:])/len(netcdf['38062/u'][i][:])});
u_total_df = pandas.DataFrame(data=total_u)
for i in range(1, len(Time_series)):
    total_U.append({'timestep': i, 'total':sum(netcdf['38062/U'][i][oifs_range])/len(netcdf['38062/U'][i][oifs_range])});
U_total_df = pandas.DataFrame(data=total_U)
####
total_fthl = []
total_fT = []
for i in range(1, len(Time_series)):
    total_fthl.append({'timestep': i, 'total':sum(netcdf['38062/f_thl'][i][:])/len(netcdf['38062/f_thl'][i][:])});
fthl_total_df = pandas.DataFrame(data=total_fthl)
for i in range(1, len(Time_series)):
    total_fT.append({'timestep': i, 'total':sum(netcdf['38062/f_T'][i][oifs_range])/len(netcdf['38062/f_T'][i][oifs_range])});
fT_total_df = pandas.DataFrame(data=total_fT)
####
total_fu = []
total_fU = []
for i in range(1, len(Time_series)):
    total_fu.append({'timestep': i, 'total':sum(netcdf['38062/f_u'][i][:])/len(netcdf['38062/f_u'][i][:])});
fu_total_df = pandas.DataFrame(data=total_fu)
for i in range(1, len(Time_series)):
    total_fU.append({'timestep': i, 'total':sum(netcdf['38062/f_U'][i][oifs_range])/len(netcdf['38062/f_U'][i][oifs_range])});
fU_total_df = pandas.DataFrame(data=total_fU)
####
total_fv = []
total_fV = []
for i in range(1, len(Time_series)):
    total_fv.append({'timestep': i, 'total':sum(netcdf['38062/f_v'][i][:])/len(netcdf['38062/f_v'][i][:])});
fv_total_df = pandas.DataFrame(data=total_fv)
for i in range(1, len(Time_series)):
    total_fV.append({'timestep': i, 'total':sum(netcdf['38062/f_V'][i][oifs_range])/len(netcdf['38062/f_V'][i][oifs_range])});
fV_total_df = pandas.DataFrame(data=total_fV)
####
total_fSH = []
for i in range(1, len(Time_series)):
    total_fSH.append({'timestep': i, 'total':sum(netcdf['38062/f_SH'][i][oifs_range])/len(netcdf['38062/f_SH'][i][oifs_range])});
total_fSH_df = pandas.DataFrame(data=total_fSH)

total_fqt = []
for i in range(1, len(Time_series)):
    total_fqt.append({'timestep': i, 'total':sum(netcdf['38062/f_qt'][i])/len(netcdf['38062/f_qt'][i])});
total_fqt_df = pandas.DataFrame(data=total_fqt)
####
total_A = []
for i in range(1, len(Time_series)):
    total_A.append({'timestep':i, 'Atot':sum(netcdf['38062/A'][i])})
total_A_df =pandas.DataFrame(data=total_A)
####
total_lw = []
total_LW = []
R = 8.3144621
for i in range(1, len(Time_series)):
    total_lw.append({'timestep': i, 'total':sum(netcdf['38062/ql'][i][:]*(netcdf['38062/presf'][i][:]/netcdf['38062/t'][i][:]/R)*diff_heights[:])});
lw_total_df = pandas.DataFrame(data=total_lw)

for i in range(1, len(Time_series)):
    total_LW_2.append({'timestep': i, 'total':sum(netcdf['38062/QL'][i][oifs_range]*(netcdf['38062/Pf'][i][oifs_range]/netcdf['38062/T'][i][oifs_range]/R)*diff_Heights[i][oifs_range])});
LW_total_df = pandas.DataFrame(data=total_LW)
####
total_tw = []
total_TW = []
for i in range(1, len(Time_series)):
    total_tw.append({'timestep': i, 'total':sum(netcdf['38062/qt'][i][:]*(netcdf['38062/presf'][i][:]/netcdf['38062/t'][i][:]/R)*diff_heights[:])});
tw_total_df = pandas.DataFrame(data=total_tw)

for i in range(1, len(Time_series)):
    total_TW.append({'timestep': i, 'total':sum(netcdf['38062/QT'][i][oifs_range]*(netcdf['38062/Pf'][i][oifs_range]/netcdf['38062/T'][i][oifs_range]/R)*diff_Heights[i][oifs_range])});
TW_total_df = pandas.DataFrame(data=total_TW)
####
ref_heights = [0, 4, 8, 12, 40, 80, 120, 159]
height_series[ref_heights]
Ref_heights = [90, 87, 85, 84, 78, 73, 70, 67]
Height_series[2][Ref_heights]



## PLOTS


pyplot.savefig(outfile+'.png', dpi=600)
pyplot.savefig(outfile+'.svg')
pyplot.savefig(outfile+'.eps')