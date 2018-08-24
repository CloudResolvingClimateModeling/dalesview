#!/usr/bin/env python

from __future__ import print_function
import os
from dalesdata import dalesdata
from dalesview import dalesview
import logging

# Basic example of what dalesview currently can do.

data_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "testdata", "H002")
logging.basicConfig(level=logging.DEBUG)

if __name__ == "__main__":
    data = dalesdata.DalesData(data_path, 2)
    print(data)
    view = dalesview.DalesView(data)

    view.plot("tke", "wmax", "zi", "wthvr", "v2r", "w2r", "thl2r", "thv", "skew", "obukh", "uws", "thlskin")
