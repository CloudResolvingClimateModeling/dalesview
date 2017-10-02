import os
import logging
import unittest
import dalesdata
import nose.tools

logging.basicConfig(level = logging.DEBUG)
basedir = os.path.dirname(dalesdata.__file__)

class dalesinput_tests(unittest.TestCase):

    def test_input_files(self):
        input = dalesdata.dalesinput(dalesdir = os.path.join(basedir,"testdata","H001"))
        filenames = set([os.path.basename(f.filepath) for f in input.filereaders])
        assert filenames == set(["baseprof.inp.001","lscale.inp.001","prof.inp.001"])

    def test_input_variables(self):
        input = dalesdata.dalesinput(dalesdir = os.path.join(basedir,"testdata","H001"))
        profvars = set(input.profiles)
        assert profvars == set(["rhobf","ugeo","vgeo","wfls","dqtdtls","dthldt","thl","qt","u","v","tke","not_used"])

    def test_input_heights(self):
        input = dalesdata.dalesinput(dalesdir = os.path.join(basedir,"testdata","H001"))
        wfls = input.profiles["wfls"]
        assert len(wfls.heights) == 128

    def test_input_times(self):
        input = dalesdata.dalesinput(dalesdir = os.path.join(basedir,"testdata","H001"))
        wfls = input.profiles["wfls"]
        assert wfls.times == [0]

    def test_input_profile(self):
        input = dalesdata.dalesinput(dalesdir = os.path.join(basedir,"testdata","H001"))
        tke = input.profiles["tke"]
        print tke[1]
        assert tke[0] == 0.0008 and tke[-1] == 0.
