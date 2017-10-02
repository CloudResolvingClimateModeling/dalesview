import os
import logging
import numpy
import unittest
import dalesdata
import nose.tools

logging.basicConfig(level = logging.DEBUG)
basedir = os.path.dirname(dalesdata.__file__)

class dalesoutput_tests(unittest.TestCase):

    def test_output_files(self):
        output = dalesdata.dalesoutput(dalesdir = os.path.join(basedir,"testdata","H001"))
        filenames = set([os.path.basename(f.filepath) for f in output.filereaders])
        assert filenames == set(["profiles.001.nc","tmser.001.nc"])

    def test_profile_variables(self):
        output = dalesdata.dalesoutput(dalesdir = os.path.join(basedir,"testdata","H001"))
        profvars = set(output.profiles)
        assert "thv" in profvars and "zt" not in profvars

    def test_timeseries_variables(self):
        output = dalesdata.dalesoutput(dalesdir = os.path.join(basedir,"testdata","H001"))
        seriesvars = set(output.timeseries)
        assert "ustar" in seriesvars and "time" not in seriesvars

    def test_profile_heights(self):
        output = dalesdata.dalesoutput(dalesdir = os.path.join(basedir,"testdata","H001"))
        thv = output.profiles["thv"]
        assert len(thv.heights) == 128

    def test_input_times(self):
        output = dalesdata.dalesoutput(dalesdir = os.path.join(basedir,"testdata","H001"))
        thlskin = output.timeseries["thlskin"]
        print thlskin.times
        assert (thlskin.times == numpy.array([60.,120.,180.,240.,300.])).all()
