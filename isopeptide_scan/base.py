"""
The class ``BaseScanner`` is a base class for the ``IsopetideScanner`` class.
"""

import types
from typing import (Optional, Dict)
import pyrosetta
from .common import IsopetideType, create_weighted_scorefxn
from . import pose_operations as pose_ops

prc: types.ModuleType = pyrosetta.rosetta.core
Residue: types.ModuleType = prc.conformation.Residue


# ---------------------------------------------------------------------------------------------------------

class BaseScanner:
    """
    :cvar residue_names: dictionary name3 to the Rosetta patched name3 for sidechain conjugation
    :cvar atom_names: dictionary name3 to the atom name for the isopeptide bond connection
    :ivar mode:  the wanted isopeptide type, the enum ``IsopetideType``, from __init__ argument
    :ivar distance: the default distance of the isopeptide bond for the method ``create_ideal_isopeptide``
    :ivar ideal_template: a pose containing only the relevant residues at the stated distance. Gets used as playdough
    :ivar scorefxn: a score function with a dihedral constraint term, see ``common.create_weighted_scorefxn``
    """
    residue_names = {'ASP': 'D[ASP:SidechainConjugation]',
                     'LYS': 'K[LYS:SidechainConjugation]',
                     }
    atom_names = {'ASP': 'CG', 'LYS': 'NZ'}

    def __init__(self,
                 pose: pyrosetta.Pose,
                 mode: IsopetideType = IsopetideType.ASP_LYS_NONE,
                 distance: int = 9):
        """
        :param pose: the pose to scan
        :param mode: the wanted isopeptide type, the enum ``IsopetideType``
        :param distance: the default distance of the isopeptide bond for the method ``create_ideal_isopeptide``
        """
        self.pose = pose
        self.mode: IsopetideType = mode
        self.distance = distance
        self.scorefxn_cart = create_weighted_scorefxn({
                                                    prc.scoring.ScoreType.coordinate_constraint: 1,
                                                    prc.scoring.ScoreType.atom_pair_constraint: 5,
                                                    prc.scoring.ScoreType.angle_constraint: 5,
                                                    prc.scoring.ScoreType.dihedral_constraint: 5,
                                                  },
                                                 score_name='ref2015_cart')
        self.scorefxn_cart_ex = create_weighted_scorefxn({
            prc.scoring.ScoreType.coordinate_constraint: 1,
            prc.scoring.ScoreType.atom_pair_constraint: 5,
            prc.scoring.ScoreType.angle_constraint: 5,
            prc.scoring.ScoreType.dihedral_constraint: 5,
            prc.scoring.ScoreType.fa_sol: 0,  # ignore solvent
            prc.scoring.ScoreType.fa_dun: 0,  # ignore rotamers
        },
            score_name='ref2015_cart')
        self.scorefxn_internal = create_weighted_scorefxn({
            prc.scoring.ScoreType.coordinate_constraint: 1,
            prc.scoring.ScoreType.atom_pair_constraint: 5,
            prc.scoring.ScoreType.angle_constraint: 5,
            prc.scoring.ScoreType.dihedral_constraint: 5,
        },
            score_name='ref2015')
        # this gets used as playdough:
        self.ideal_template: pyrosetta.Pose = self.create_ideal_isopeptide(distance)

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
        The returned pose is just two residues connected and constrained.
        See ``_create_template`` for more.

        The distance is enforced by a constraint on the CA atoms.
        """
        distance = distance if distance else self.distance
        pose = self._create_template()
        pose_ops.remove_constraints(pose)
        pose_ops.constrain_CA_distance(pose, distance)
        pose_ops.constrain_isopeptide(pose)
        pose_ops.relax_weighted_cart(pose, self.scorefxn_cart_ex)
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
        pose_ops.remove_termini(test, (1, 2))
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
