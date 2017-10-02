import os

# Basic data objects in the package: profiles (possibly time dependent) and time
# series. No xy-dependent data is currently supported.

# Class holding profiles for a single variable.
class profile(object):

    def __init__(self,var,filereader):
        self.variable = var
        self.unit = filereader.get_unit(var)
        self.filereader = filereader
        self.heights = filereader.get_heights(var)
        self.times = filereader.get_times(var)
        self.shape = filereader.get_shape(var)
        self.missvals = {}
        for v in [var,"height","time"]:
            self.missvals[v] = filereader.get_missval(v)

    def __getitem__(self,*args):
        return self.filereader.__getitem__(self.variable,*args)

    def __str__(self,*args):
        shapestr = ""
        if(len(self.shape) > 1):
            shapestr = "(t = {numt},z = {numz})".format(numt = self.shape[0],numz = self.shape[1])
        elif(len(self.shape) == 1):
            shapestr = "(z = {numz})".format(numz = self.shape[0])
        return "{varname}{shape} from {fname}".format(varname = self.variable,shape = shapestr,fname = os.path.basename(self.filereader.filepath))


# Class holding time series for a single variable.
class timeseries(object):

    def __init__(self,var,filereader):
        self.variable = var
        self.unit = filereader.get_unit(var)
        self.filereader = filereader
        self.times = filereader.get_times(var)
        self.shape = filereader.get_shape(var)
        self.missvals = {}
        for v in [var,"time"]:
            self.missvals[v] = filereader.get_missval(v)

    def __getitem__(self,*args):
        return self.filereader.__getitem__(self.variable,*args)

    def __str__(self,*args):
        shapestr = ""
        if(len(self.shape) == 1):
            shapestr = "(t = {numt})".format(numt = self.shape[0])
        return "{varname}{shape} from {fname}".format(varname = self.variable,shape = shapestr,fname = os.path.basename(self.filereader.filepath))
