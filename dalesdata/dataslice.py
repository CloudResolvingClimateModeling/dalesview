import os


# Basic data objects in the package: profiles (possibly time dependent) and time
# series. No xy-dependent data is currently supported.

# Class holding profiles for a single variable.
class Profile(object):

    def __init__(self, var, file_reader):
        self.variable = var
        self.unit = file_reader.get_unit(var)
        self.file_reader = file_reader
        self.heights = file_reader.get_heights(var)
        self.times = file_reader.get_times(var)
        self.shape = file_reader.get_shape(var)
        self.miss_vals = {}
        for v in [var, "height", "time"]:
            self.miss_vals[v] = file_reader.get_missval(v)

    def __getitem__(self, *args):
        return self.file_reader.__getitem__(self.variable, *args)

    def __str__(self, *args):
        shape_string = ""
        if len(self.shape) > 1:
            shape_string = "(t = {numt},z = {numz})".format(numt=self.shape[0], numz=self.shape[1])
        elif len(self.shape) == 1:
            shape_string = "(z = {numz})".format(numz=self.shape[0])
        return "{varname}{shape} from {fname}".format(varname=self.variable, shape=shape_string,
                                                      fname=os.path.basename(self.file_reader.filepath))


# Class holding time series for a single variable.
class TimeSeries(object):

    def __init__(self, var, file_reader):
        self.variable = var
        self.unit = file_reader.get_unit(var)
        self.file_reader = file_reader
        self.times = file_reader.get_times(var)
        self.shape = file_reader.get_shape(var)
        self.miss_vals = {}
        for v in [var, "time"]:
            self.miss_vals[v] = file_reader.get_missval(v)

    def __getitem__(self, *args):
        return self.file_reader.__getitem__(self.variable, *args)

    def __str__(self, *args):
        shape_string = ""
        if len(self.shape) == 1:
            shape_string = "(t = {numt})".format(numt=self.shape[0])
        return "{varname}{shape} from {fname}".format(varname=self.variable, shape=shape_string,
                                                      fname=os.path.basename(self.file_reader.filepath))
