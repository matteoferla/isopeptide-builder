import pyrosetta
import types
from typing import (Tuple)
from .base import BaseScanner
from .pose_operations import PoseOps

prc: types.ModuleType = pyrosetta.rosetta.core
Residue: types.ModuleType = prc.conformation.Residue


class IsopetideScanner(BaseScanner):
    """
    :ivar mode:  isopeptide or not, the enum ``IsopetideType``
    :cvar residue_names: shorthand for the sidechain annotated names... mod if new residues needed
    :ivar distance: distance of _ideal_ isopetide.
    :ivar ideal_template: a pose containing only the relevant residues at the stated distance
    """

    def __call__(self, asx: Residue, lys: Residue) -> float:
        """
        Given two residues return the score of an isopeptide at those positions.
        The two residues givem of type ``pyrosetta.rosetta.core.conformation.Residue``,
        have to be L alpha amino acid residues, but not necessarily asparax (asx) and lysine (lys).
        The two are not interchangeable, as the amide bond is closer to the backbone of the
        asparax than the lysine.

        Calling ``.is_suitable`` first is recommended.

        :param asx:
        :param lys:
        :return:
        """
        # if not is_suitable(i, j):
        #     return float('nan')
        test = self.ideal_template.clone()
        test_ops = PoseOps(test)
        self.copy_backbone(donor_residue=asx, acceptor_residue=test.residue(1))
        self.copy_backbone(donor_residue=lys, acceptor_residue=test.residue(2))
        test_ops.relax_sidechains(self.scorefxn, 3)
        return self.scorefxn(test)

    def is_suitable(self,
                    resi_A: int, resi_B: int,
                    distance_interval: Tuple[int, int] = (3.5, 10.5),
                    tolerance: int = 1) -> bool:
        """
        Do not waste time calculating the score of a pair of residues that are not suitable.

        * Is it a peptide?
        * Is distance between alpha atoms between the specified?
        * Is distance between beta atoms is less than the distance of alpha atom (with added tolerance)?
        """
        if not self.pose.residue(resi_A).is_protein() or not self.pose.residue(resi_B).is_protein():
            return False
        ca_d = (self.pose.residue(resi_A).atom('CA').xyz() - self.pose.residue(resi_B).atom('CA').xyz()).norm()
        if distance_interval[0] > ca_d or ca_d > distance_interval[1]:
            return False
        # glycine in Rosetta has 1HA and 2HA, not HA1 & HA2
        beta_A = 'CB' if self.pose.residue(resi_A).name3() != 'GLY' else '1HA'
        beta_B = 'CB' if self.pose.residue(resi_B).name3() != 'GLY' else '1HA'
        cb_d = (self.pose.residue(resi_A).atom(beta_A).xyz() - self.pose.residue(resi_B).atom(beta_B).xyz()).norm()
        return ca_d + tolerance >= cb_d
