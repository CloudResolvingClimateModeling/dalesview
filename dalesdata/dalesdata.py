import os
import glob
import logging
from . import dalesreader
from . import dataslice


# Data objects for the DALES view package. The objects contain dictionaries for
# profiles and/or time series for all input and output of a given experiment.

# Object holding all data from a DALES run.
class DalesData(object):

    # Creates DALES data object for the given run directory and experiment number
    def __init__(self, dalesdir=".", exp=1):
        if not os.path.isdir(dalesdir):
            raise Exception("Dales directory %s does not exist" % dalesdir)
        self.path = dalesdir
        self.input = DalesInput(dalesdir, exp)
        self.output = DalesOutput(dalesdir, exp)

    def __str__(self):
        return """Data directory: {path}

Input:
------
{input}
Output:
------
{output}""".format(path=self.path, input=self.input.__str__(), output=self.output.__str__())


# Object holding all input data profiles a DALES run.
class DalesInput(object):

    # Creates DALES input data object for the given run directory and experiment number
    def __init__(self, dalesdir=".", exp=1):
        self.filereaders = []
        self.profiles = {}
        make_profiles(glob.glob(os.path.join(dalesdir, ".".join(["*", "inp", str(exp).zfill(3)]))), self.filereaders,
                      self.profiles)

    def __str__(self):
        result = ""
        for p in self.profiles.values():
            result += "Profile " + p.__str__() + os.linesep
        return result


# Object holding all output data from a DALES run.
class DalesOutput(object):

    # Creates DALES output data object for the given run directory and experiment number
    def __init__(self, dalesdir=".", exp=1):
        self.filereaders = []
        self.profiles = {}
        self.timeseries = {}
        expstr = str(exp).zfill(3)
        asciifiles = glob.glob(os.path.join(dalesdir, ".".join(["*", expstr])))
        inputfiles = glob.glob(os.path.join(dalesdir, ".".join(["*", "inp", expstr])))
        netcdfiles = glob.glob(os.path.join(dalesdir, ".".join(["*", expstr, "nc"])))
        outputfiles = list(set(asciifiles) - set(inputfiles)) + netcdfiles
        # we only support netcdf output for now
        timseriesfiles = [f for f in outputfiles if os.path.basename(f).startswith("tmser")]
        make_timeseries([f for f in timseriesfiles if f in netcdfiles], self.filereaders, self.timeseries)
        profilefiles = [f for f in outputfiles if os.path.basename(f).startswith("profiles")]
        make_profiles([f for f in profilefiles if f in netcdfiles], self.filereaders, self.profiles)

    def __str__(self):
        result = ""
        for p in self.profiles.values():
            result += "Profile " + p.__str__() + os.linesep
        for t in self.timeseries.values():
            result += "Time series " + t.__str__() + os.linesep
        return result


# Adds profile data objects to the argument dictionary profiles,
# and appends all created file readers to the argument list readers.
def make_profiles(files, readers, profiles):
    for f in files:
        reader = dalesreader.make_file_reader(f)
        readers.append(reader)
        for v in list(set(reader.variables) - set(reader.dimensions)):
            profiles[v] = dataslice.Profile(v, reader)


# Adds time series data objects to the argument dictionary timeseries,
# and appends all created file readers to the argument list readers.
def make_timeseries(files, readers, timeseries):
    for f in files:
        reader = dalesreader.make_file_reader(f)
        readers.append(reader)
        for v in list(set(reader.variables) - set(reader.dimensions)):
            timeseries[v] = dataslice.TimeSeries(v, reader)
