import itertools
import pyrosetta
import pyrosetta_help as ph
import types
from typing import *
from base import IsopetideType

prc: types.ModuleType = pyrosetta.rosetta.core
Residue: types.ModuleType = prc.conformation.Residue




# ----- test isopeptide ------------------------------









class QuickLayoutMetric:
    """
    This is not a Rosetta score function, but a wrapper to it
    """

    # no... this should not be stored:
    # original_test_pose: Union[None, pyrosetta.Pose] = None

    def __init__(self, **options):
        self.ideal_distance: int = options.get('ideal_distance', 9)
        self.test_template = self.create_ideal_isopeptide(self.ideal_distance)
        self.test_template.remove_constraints()
        constrain_isopeptide(self.test_template)


