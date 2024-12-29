import unittest
from isopeptide_scan import IsopetideScanner
import isopeptide_scan.pose_operations as pose_ops
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

class TestBaseOps(unittest.TestCase):

    def setUp(self):
        if os.path.exists('1ubq.pdb'):
            filename = '1ubq.pdb'
        else:
            filename:str = ph.download_pdb('1ubq')
        self.pose:pyrosetta.Pose = pyrosetta.pose_from_file(filename)  # noqa
        pyrosetta.rosetta.core.pose.remove_nonprotein_residues(self.pose)

    def test_basics(self):
        self.assertGreater(len(self.pose.residues), 0)

    def test_extend(self):
        test = self.pose.clone()
        pass

    def test_has_iso(self):
        # this is a isopeptide bonded pair of residues only
        scanner = IsopetideScanner(self.pose)
        ideal: pyrosetta.Pose = scanner.create_ideal_isopeptide(8)
        self.assertEqual(ideal.sequence(), 'DK')
        # does it have a non-polymeric connection?
        self.assertEqual(ideal.residue(1).n_non_polymeric_residue_connections(), 1)
        self.assertEqual(ideal.residue(2).n_non_polymeric_residue_connections(), 1)
        # is it okay?
        self.assertFalse(ideal.residue(1).connection_incomplete(3))
        self.assertFalse(ideal.residue(2).connection_incomplete(3))
        # is it mutual
        self.assertEqual(ideal.residue(1).connected_residue_at_resconn(3), 2)
        self.assertEqual(ideal.residue(2).connected_residue_at_resconn(3), 1)
        # is there a constraint?
        cons = [c.to_string().split('\n')[0] for c in ideal.constraint_set().get_all_constraints()]
        print(cons)  # 2x DihedralConstraint, 1x AtomPairConstraint
        self.assertEqual(len(cons), 3)

    def test_pair(self):
        scanner = IsopetideScanner(self.pose)
        ideal: pyrosetta.Pose = scanner.create_ideal_isopeptide(8)
        score = scanner.score_isopeptide_plausibility(ideal.residue(1), ideal.residue(2))
        self.assertLess(score, 0)




