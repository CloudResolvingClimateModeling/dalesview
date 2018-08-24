import logging
import os
import unittest

from dalesdata.dalesdata import DalesInput

logging.basicConfig(level=logging.DEBUG)
datadir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "testdata"))


class DalesInputTests(unittest.TestCase):

    @staticmethod
    def test_input_files():
        dales_input = DalesInput(dalesdir=os.path.join(datadir, "H001"))
        filenames = set([os.path.basename(f.filepath) for f in dales_input.filereaders])
        assert filenames == {"baseprof.inp.001", "lscale.inp.001", "prof.inp.001"}

    @staticmethod
    def test_input_variables():
        dales_input = DalesInput(dalesdir=os.path.join(datadir, "H001"))
        profvars = set(dales_input.profiles)
        assert profvars == {"rhobf", "ugeo", "vgeo", "wfls", "dqtdtls", "dthldt", "thl", "qt", "u", "v", "tke",
                            "not_used"}

    @staticmethod
    def test_input_heights():
        dales_input = DalesInput(dalesdir=os.path.join(datadir, "H001"))
        wfls = dales_input.profiles["wfls"]
        assert len(wfls.heights) == 128

    @staticmethod
    def test_input_times():
        dales_input = DalesInput(dalesdir=os.path.join(datadir, "H001"))
        wfls = dales_input.profiles["wfls"]
        assert wfls.times == [0]

    @staticmethod
    def test_input_profile():
        dales_input = DalesInput(dalesdir=os.path.join(datadir, "H001"))
        tke = dales_input.profiles["tke"]
        assert tke[0] == 0.0008 and tke[-1] == 0.
