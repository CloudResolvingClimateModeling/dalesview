import dalesdata
import dataslice
import logging

# DALES view object. Allows simple but easy plotting of dales data objects.

log = logging.getLogger(__name__)

# DALES view class. Contains a data object and a plotting backend class.
class dalesview(object):

    mpl = 1

    def __init__(self,data,backend = mpl):
        self.data = data
        self.plotapi = dalesview.create_plotapi(backend)

    # Plotting method, allows arbitrary argument lists of variable strings or
    # data objects (input/output/profiles etc.)
    def plot(self,*args,**kwargs):
        dataproviders = list(self.data.output.profiles.values())\
                        + list(self.data.output.timeseries.values())\
                        + list(self.data.input.profiles.values())
        dataprovs = []
        for arg in args:
            include = False
            for d in dataproviders:
                if(dalesview.contains_data(arg,d)):
                    dataprovs.append(d)
                    include = True
            if(not include):
                log.error("Variable %s could not be found in given DALES data directory %s" % (str(arg),self.data.path))
        if(not any(dataprovs)):
            return self.plotapi.multiplot(dataproviders,kwargs) if not any(args) else None
        elif(len(dataprovs) == 1):
            return self.plotapi.plot(dataprovs[0],kwargs)
        else:
            return self.plotapi.multiplot(list(set(dataprovs)),**kwargs)

    # Factory method for creating the plotting backend. Currently only netcdf
    # backend is supported.
    @staticmethod
    def create_plotapi(backend_type):
        if(backend_type == dalesview.mpl):
            import mpl4dales
            return mpl4dales.mpl4dales()
        else:
            log.error("Plotting backend type %s is not known" % str(backend_type))
        return None

    # Helper method to match the data argument in the first argument.
    @staticmethod
    def contains_data(arg,data):
        if(arg == data):
            return True
        if(isinstance(arg,str)):
            return arg == data.variable
        if(isinstance(arg,dalesdata.dalesdata)):
            return contains_data(arg.input,data) or contains_data(arg.output,data)
        if(isinstance(arg,dalesdata.dalesinput)):
            return data in arg.profiles.values()
        if(isinstance(arg,dalesdata.dalesoutput)):
            return data in arg.profiles.values() + arg.timeseries.values()
        if(isinstance(arg,dict)):
            return data in arg.values()
        if(isinstance(arg,list)):
            return data in arg
        return False
