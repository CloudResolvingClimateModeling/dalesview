import logging
import os
import unittest

import numpy

from dalesdata.dalesdata import DalesOutput

logging.basicConfig(level=logging.DEBUG)
datadir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "testdata"))


class DalesOutputTests(unittest.TestCase):

    @staticmethod
    def test_output_files():
        output = DalesOutput(dalesdir=os.path.join(datadir, "H001"))
        filenames = set([os.path.basename(f.filepath) for f in output.filereaders])
        assert filenames == {"profiles.001.nc", "tmser.001.nc"}

    @staticmethod
    def test_profile_variables():
        output = DalesOutput(dalesdir=os.path.join(datadir, "H001"))
        profvars = set(output.profiles)
        assert "thv" in profvars and "zt" not in profvars

    @staticmethod
    def test_timeseries_variables():
        output = DalesOutput(dalesdir=os.path.join(datadir, "H001"))
        seriesvars = set(output.timeseries)
        assert "ustar" in seriesvars and "time" not in seriesvars

    @staticmethod
    def test_profile_heights():
        output = DalesOutput(dalesdir=os.path.join(datadir, "H001"))
        thv = output.profiles["thv"]
        assert len(thv.heights) == 128

    @staticmethod
    def test_input_times():
        output = DalesOutput(dalesdir=os.path.join(datadir, "H001"))
        thlskin = output.timeseries["thlskin"]
        assert (thlskin.times == numpy.array([60., 120., 180., 240., 300.])).all()
