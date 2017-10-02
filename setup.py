#!/usr/bin/env python

from setuptools import setup

setup(name = "dalesview",
      version = "0.1",
      description = "A python package for visualizing DALES input/output data",
      author = "Gijs van den Oord",
      author_email = "g.vandenoord@esciencecenter.nl",
      url = "https://github.com/goord/dalesview",
      packages = ["dalesview","dalesdata"],
      setup_requires = ["numpy"],
      install_requires = ["numpy","matplotlib","netCDF4"])
