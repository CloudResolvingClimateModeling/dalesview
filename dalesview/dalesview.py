import logging

from dalesdata import dalesdata

# DALES view object. Allows simple but easy plotting of dales data objects.

log = logging.getLogger(__name__)


# DALES view class. Contains a data object and a plotting backend class.
class DalesView(object):
    mpl = 1

    def __init__(self, data, backend=mpl):
        self.data = data
        self.plot_api = DalesView.create_plot_api(backend)

    # Plotting method, allows arbitrary argument lists of variable strings or
    # data objects (input/output/profiles etc.)
    def plot(self, *args, **kwargs):
        data_providers = list(self.data.output.profiles.values()) \
                         + list(self.data.output.timeseries.values()) \
                         + list(self.data.input.profiles.values())
        dataprovs = []
        for arg in args:
            include = False
            for d in data_providers:
                if DalesView.contains_data(arg, d):
                    dataprovs.append(d)
                    include = True
            if not include:
                log.error(
                    "Variable %s could not be found in given DALES data directory %s" % (str(arg), self.data.path))
        if not any(dataprovs):
            return self.plot_api.multiplot(data_providers, **kwargs) if not any(args) else None
        elif len(dataprovs) == 1:
            return self.plot_api.plot(dataprovs[0], **kwargs)
        else:
            return self.plot_api.multiplot(list(set(dataprovs)), **kwargs)

    # Factory method for creating the plotting backend. Currently only netcdf
    # backend is supported.
    @staticmethod
    def create_plot_api(backend_type):
        if backend_type == DalesView.mpl:
            from . import mpl4dales
            return mpl4dales.Mpl4Dales()
        else:
            log.error("Plotting backend type %s is not known" % str(backend_type))
        return None

    # Helper method to match the data argument in the first argument.
    @staticmethod
    def contains_data(arg, data):
        if arg == data:
            return True
        if isinstance(arg, str):
            return arg == data.variable
        if isinstance(arg, dalesdata.DalesData):
            return DalesView.contains_data(arg.input, data) or DalesView.contains_data(arg.output, data)
        if isinstance(arg, dalesdata.DalesInput):
            return data in arg.profiles.values()
        if isinstance(arg, dalesdata.DalesOutput):
            return data in arg.profiles.values() + arg.timeseries.values()
        if isinstance(arg, dict):
            return data in arg.values()
        if isinstance(arg, list):
            return data in arg
        return False
