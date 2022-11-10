import unittest
from isopeptide_scan import IsopetideScanner, PoseOps
import pyrosetta_help as ph
import pyrosetta
import os


logger = ph.configure_logger()
pyrosetta.distributed.maybe_init(extra_options=ph.make_option_string(no_optH=False,
                                                ex1=None,
                                                ex2=None,
                                                #mute='all',
                                                ignore_unrecognized_res=True,
                                                load_PDB_components=False,
                                                ignore_waters=False)
                                 )

class Start:
    def setUp(self):
        if os.path.exists('1ubq.pdb'):
            filename = '1ubq.pdb'
        else:
            filename:str = ph.download_pdb('1ubq')
        self.ref:pyrosetta.Pose = pyrosetta.pose_from_file(filename)  # noqa
        pyrosetta.rosetta.core.pose.remove_nonprotein_residues(self.ref)

class TestBaseOps(Start, unittest.TestCase):

    def test_basics(self):
        self.assertGreater(len(self.ref.residues), 0)

    def test_extend(self):
        test = self.ref.clone()
        test_ops = PoseOps(test)
        pass

    def test_has_iso(self):
        iso = IsopetideScanner(self.ref).create_ideal_isopeptide()



