import types
from typing import (Optional, Dict)
import pyrosetta
from .common import IsopetideType, create_weighted_scorefxn
from .pose_operations import PoseOps

prc: types.ModuleType = pyrosetta.rosetta.core
Residue: types.ModuleType = prc.conformation.Residue


# ---------------------------------------------------------------------------------------------------------

class BaseScanner:
    residue_names = {'ASP': 'D[ASP:SidechainConjugation]',
                     'LYS': 'K[LYS:SidechainConjugation]',
                     }
    atom_names = {'ASP': 'CG', 'LYS': 'NZ'}

    def __init__(self,
                 pose: pyrosetta.Pose,
                 mode: IsopetideType = IsopetideType.ASP_LYS_NONE,
                 distance: int = 9):
        self.pose = pose
        self.mode: IsopetideType = mode
        self.distance = distance
        self.ideal_template = self.create_ideal_isopeptide(distance)
        self.scorefxn = create_weighted_scorefxn({prc.scoring.ScoreType.dihedral_constraint: 5})

    @property
    def acid(self):
        return self.mode.name.split('_')[0]

    @property
    def base(self):
        return self.mode.name.split('_')[1]

    @property
    def catalyst(self):
        return self.mode.name.split('_')[2]

    # ---------------------------------------------------------------------------------------------------------

    def create_ideal_isopeptide(self, distance: Optional[float] = None) -> pyrosetta.Pose:
        """
        THis is just two residues connected and constrained
        """
        distance = distance if distance else self.distance
        pose = self._create_template()
        pose_ops = PoseOps(pose)
        pose_ops.remove_constraints()
        pose_ops.constrain_isopeptide()
        pose_ops.constrain_CA_distance(distance)
        pose_ops.relax_weighted_cart(self.scorefxn)
        return pose

    def _create_template(self) -> pyrosetta.Pose:
        """
        A test pose for the isopeptide bond containing solely
        an ASP:SidechainConjugation and LYS:SidechainConjugation
        residues that are bonded by an isopeptide and separated by a jump in terms of BB.
        """
        test = pyrosetta.pose_from_sequence(self.residue_names[self.acid])
        new = pyrosetta.pose_from_sequence(self.residue_names[self.base])
        test.append_pose_by_jump(new, 1)
        test.conformation().declare_chemical_bond(seqpos1=1,
                                                  atom_name1=self.atom_names[self.acid],
                                                  seqpos2=2,
                                                  atom_name2=self.atom_names[self.base])
        PoseOps(test).remove_termini((1, 2))
        return test

    # ---------------------------------------------------------------------------------------------------------

    @classmethod
    def copy_backbone(cls, donor_residue: Residue, acceptor_residue: Residue) -> None:
        """
        Copy the coordinates of ``donor_residue`` into ``acceptor_residue``.
        """
        donor_bb_atom_names = [donor_residue.atom_name(ai) for ai in donor_residue.all_bb_atoms()]
        acceptor_bb_atom_names = [acceptor_residue.atom_name(ai) for ai in acceptor_residue.all_bb_atoms()]
        for atomname in set(donor_bb_atom_names).intersection(set(acceptor_bb_atom_names)):
            acceptor_residue.set_xyz(atomname, donor_residue.xyz(atomname))
