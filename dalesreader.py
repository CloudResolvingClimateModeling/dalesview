import logging
import re
import numpy
import netCDF4

# File readers for DALES experiments. We support single-timestep prifile ASCII
# file and netcdf files for time-dependent output profiles or time series.

log = logging.getLogger(__name__)

# Factory method for cerating file readers.
def make_filereader(filename):
    result = None
    if(filename.endswith(".nc")):
        result = netcdf_reader(filename)
    else:
        result = ascii_reader(filename)
    log.info("Detected %d variables in file %s" % (len(result.variables),filename))
    return result


# Base class for the DALES file readers.
class dales_reader(object):

    def __init__(self,filepath):
        self.filepath = filepath
        self.variables = []
        self.dimensions = []
        self.units = []
        self.data = []

    # Returns the heights read from the file. If the file only contains time series,
    # returns a list with a single zero height.
    def get_heights(self,var):
        return []

    # Returns the times read from the file in seconds since the starttime.
    # If the file only contains time-independent profiles, returns a list with
    # a single zero time step.
    def get_times(self,var):
        return []

    # Returns the unit of the given variable.
    def get_unit(self,var):
        if(var in self.variables):
            index = self.variables.index(var)
            return self.units[index]
        else:
            return "1"

    # Returns the shape of the given variable array.
    def get_shape(self,var):
        return []

    # Returns the missing value for the givan variable
    def get_missval(self,var):
        return None


# ASCII file reader. Currently, we only support single time step text files
# (initial profiles) with all text commented.
class ascii_reader(dales_reader):

    commentchar = '#'

    def __init__(self,filepath):
        super(ascii_reader,self).__init__(filepath)
        firstrow = ""
        with open(filepath,'r') as f:
            for line in f:
                row = line.strip()
                if(not any(row) or row[0] == ascii_reader.commentchar):
                    firstrow = row
                    continue
                else:
                    break
        cols = firstrow[1:].split()
        for col in cols:
            tokens = [s for s in re.split("\(|\)|\[|\]",col)]
            self.variables.append(tokens[0])
            self.units.append(tokens[1] if len(tokens) > 1 else "1")
        if("height" in self.variables):
            self.dimensions.append("height")
        self.data = numpy.loadtxt(self.filepath,comments = ascii_reader.commentchar).transpose()

    def __getitem__(self,varname,*args):
        if(varname in self.variables):
            index = self.variables.index(varname)
            newargs = (index,) + args
            return self.data.__getitem__(newargs)
        else:
            return []

    def get_heights(self,var):
        return self.__getitem__("height")

    def get_times(self,var):
        return [0]

    def get_shape(self,var):
        if(var in self.variables):
            return self.data.shape[1:]
        else:
            return super(ascii_reader,self).get_shape(var)


# Netcdf file reader implementation.
class netcdf_reader(dales_reader):

    def __init__(self,filepath):
        super(netcdf_reader,self).__init__(filepath)
        self.dataset = netCDF4.Dataset(filepath,'r')
        self.variables = [v for v in self.dataset.variables if v not in self.dataset.dimensions]
        self.units = [getattr(self.dataset.variables[v],"units",None) for v in self.variables]

    def __getitem__(self,varname,*args):
        v = self.dataset.variables.get(varname,None)
        if(not v): return []
        return v.__getitem__(*args)

    def get_heights(self,var):
        v = self.dataset.variables.get(var,None)
        if(not v): return []
        for levtype in ["zt","zm"]:
            if(levtype in v.dimensions):
                dim = self.dataset.variables.get(levtype,None)
                return dim[:] if dim else [0.]

    def get_times(self,var):
        v = self.dataset.variables.get(var,None)
        if(not v): return []
        if("time" in v.dimensions):
            dim = self.dataset.variables.get("time",None)
            return dim[:] if dim else []

    def get_shape(self,var):
        if(var in self.variables):
            return self.dataset.variables[var].shape
        else:
            return super(ascii_reader,self).get_shape(var)

    def get_missval(self,var):
        if(var in self.dataset.variables):
            return getattr(self.dataset.variables[var],"_FillValue",None)
        else:
            return super(netcdf_reader,self).get_missval(var)
