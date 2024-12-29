import pandas as pd
import pyrosetta
import types
import itertools
import pandas as pd
from typing import Tuple, Union, Sequence
from .base import BaseScanner
from . import pose_operations as pose_ops

prc: types.ModuleType = pyrosetta.rosetta.core
Residue: types.ModuleType = prc.conformation.Residue


class IsopetideScanner(BaseScanner):
    """
    :ivar mode:  isopeptide or not, the enum ``IsopetideType``
    :cvar residue_names: shorthand for the sidechain annotated names... mod if new residues needed
    :ivar distance: distance of _ideal_ isopetide.
    :ivar ideal_template: a pose containing only the relevant residues at the stated distance
    """

    def score_isopeptide_plausibility(self, asx: Residue, lys: Residue) -> float:
        """
        Given two residues return the score of a potential isopeptide at those positions.
        The two residues givem of type ``pyrosetta.rosetta.core.conformation.Residue``,
        have to be L alpha amino acid residues, but not necessarily asparax (asx) and lysine (lys).
        The two are not interchangeable, as the amide bond is closer to the backbone of the
        asparax than the lysine.

        Calling ``.is_suitable`` first is recommended.
        This costs 200 ms.

        :param asx:
        :param lys:
        :return:
        """
        test = self.adapt_isopeptide(asx, lys)
        return self.scorefxn_cart_ex(test)

    def adapt_isopeptide(self, asx: Residue, lys: Residue) -> float:
        """
        See ``score_isopeptide_plausibility`` for details.
        """
        asx = self.get_residue(asx)
        lys = self.get_residue(lys)
        d: float = asx.xyz('CA').distance(lys.xyz('CA'))
        test = self.create_ideal_isopeptide(d)
        self.copy_backbone(donor_residue=asx, acceptor_residue=test.residue(1))
        self.copy_backbone(donor_residue=lys, acceptor_residue=test.residue(2))
        pose_ops.relax_weighted_cart(test, self.scorefxn_cart_ex, bb=False, jump=False)
        #pose_ops.relax_sidechains(test, self.scorefxn_internal, 3)
        return test

    def get_residue(self, res: Union[Residue, int, Tuple[str, int]]) -> Residue:
        if isinstance(res, Residue):
            return res
        elif isinstance(res, int):
            return self.pose.residue(res)
        elif isinstance(res, tuple):
            return self.pose.pdb_rsd(res)
        else:
            raise TypeError('res must be Residue, int or tuple')

    def is_suitable(self,
                    res_A: Union[Residue, int, Tuple[str, int]],
                    res_B: Union[Residue, int, Tuple[str, int]],
                    distance_interval: Tuple[int, int] = (3.5, 10.5),
                    tolerance: int = 1) -> bool:
        """
        Do not waste time calculating the score of a pair of residues that are not suitable.

        * Is it a peptide?
        * Is distance between alpha atoms between the specified?
        * Is distance between beta atoms is less than the distance of alpha atom (with added tolerance)?
        """
        res_A: Residue = self.get_residue(res_A)
        res_B: Residue = self.get_residue(res_B)
        if not res_A.is_protein() or not res_B.is_protein():
            return False
        ca_d = (res_A.atom('CA').xyz() - res_B.atom('CA').xyz()).norm()
        if distance_interval[0] > ca_d or ca_d > distance_interval[1]:
            return False
        # glycine in Rosetta has 1HA and 2HA, not HA1 & HA2
        beta_A = 'CB' if res_A.name3() != 'GLY' else '1HA'
        beta_B = 'CB' if res_B.name3() != 'GLY' else '1HA'
        cb_d = (res_A.atom(beta_A).xyz() - res_B.atom(beta_B).xyz()).norm()
        return ca_d + tolerance >= cb_d

    def __call__(self, indices: Sequence[int]=()) -> pd.DataFrame:
        """
        If no indices are given, all pairs are scored.
        """
        if not indices:
            indices = range(1, 1 + self.pose.total_residue())
        scores = []
        _pose2pdb = self.pose.pdb_info().pose2pdb
        pose2pdb = lambda i: _pose2pdb(i).strip().replace(' ', ':')
        for i, j in itertools.product(indices, indices):
            score = {'Asx_idx1': i, 'Lys_idx1': j,
                     'Asx_pdb': pose2pdb(i), 'Lys_pdb': pose2pdb(j)}
            if i == j:
                score['verdict'] = 'identity'
            elif not self.is_suitable(i, j):
                score['verdict'] = 'unsuitable'
            else:
                s = self.score_isopeptide_plausibility(i, j)
                score['score_of_ideal'] = s
                if s > 0:
                    score['verdict'] = 'unfavorable'
                else:
                    score['verdict'] = 'favorable'
            scores.append(score)
        return pd.DataFrame(scores).sort_values('score_of_ideal')