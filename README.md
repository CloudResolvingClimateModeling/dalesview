# dalesview
A python package for visualization of DALES input/output data. DALES is the Dutch Atmospheric Large Eddy Simulation program (https://github.com/dalesteam/dales) and its input ASCII files and output NetCDF files can be quickly inspected with this library. 

### Installation
Clone the repository, activate your appropriate python environment, and in the root directory of this repo type
```
pip install .
```
or `python setup.py install`. It should automatically install the dependencies numpy, matplotlib and netCDF4 as well.

### Usage
For testing execute the example.py script on the test data. The package allows you to access the input and output data through python objects and plot them easily. E.g. if your run directory is `<rundir>` and experiment number is `<exp>`, you should be able to run
```
$ python
>>> from dalesdata.dalesdata import DalesData
>>> data = DalesData(<rundir>,<exp>)
>>> print data # lists all available variables
>>> from dalesview.dalesview import DalesView
>>> view = DalesView(data)
>>> view.plot(data.input)
>>> view.plot("thl","tke")
>>> # etc.
```
