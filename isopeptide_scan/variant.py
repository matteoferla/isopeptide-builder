import pyrosetta
import functools
from types import ModuleType
from typing import Callable
from . import pose_operations as pose_ops

prc: ModuleType = pyrosetta.rosetta.core
pr_res: ModuleType = prc.select.residue_selector
MutateResidue:Callable = pyrosetta.rosetta.protocols.simple_moves.MutateResidue


class Variant:

    def __init__(self, pose, asx_i: int, lys_i: int):
        self.pose = pose.clone()
        self.asx_i = asx_i
        self.lys_i = lys_i
        self.history = []

    @functools.cached_property
    def iso_selector(self):
        iso_selector = pr_res.ResidueIndexSelector()
        iso_selector.append_index(self.asx_i)
        iso_selector.append_index(self.lys_i)
        return iso_selector

    def mutate_to_isopeptide(self):
        MutateResidue(target=self.asx_i, new_res='ASP:SidechainConjugation').apply(self.pose)
        MutateResidue(target=self.lys_i, new_res='LYS:SidechainConjugation').apply(self.pose)
        self.pose.conformation().declare_chemical_bond(seqpos1=self.lys_i, atom_name1='NZ', seqpos2=self.asx_i, atom_name2='CG')
        pose_ops.constrain_isopeptide(self.pose, self.asx_i, self.lys_i)
        self.history.append('isopeptide')

    def get_aa_neigh_selector(self, distance=6):
        neigh_sele = pr_res.NeighborhoodResidueSelector(selector=self.iso_selector,
                                                        distance=distance,
                                                        include_focus_in_subset=False)
        canonical_sele = pr_res.ResiduePropertySelector(prc.chemical.ResidueProperty.CANONICAL_AA)
        return pr_res.AndResidueSelector(canonical_sele, neigh_sele)

    def get_adaptable_idxs(self, distance=6):
        idxs = pr_res.ResidueVector(self.get_aa_neigh_selector(distance).apply(self.pose))
        return [idx1 for idx1 in idxs if self.pose.residue(idx1).name3() not in ('PRO', 'GLY')]

    def shave_neighbors(self, distance=6):
        for idx1 in self.get_adaptable_idxs(distance):
            MutateResidue(target=idx1, new_res='ALA').apply(self.pose)
        self.history.append('shaven')

    def make_isopeptide(self, distance=6):
        self.mutate_to_isopeptide()
        self.shave_neighbors(distance)
        scorefxn_cart = pose_ops.create_weighted_scorefxn({
            prc.scoring.ScoreType.coordinate_constraint: 1,
            prc.scoring.ScoreType.atom_pair_constraint: 5,
            prc.scoring.ScoreType.angle_constraint: 5,
            prc.scoring.ScoreType.dihedral_constraint: 5,
            prc.scoring.ScoreType.fa_sol: 0,  # ignore solvent
            prc.scoring.ScoreType.fa_dun: 0,  # ignore rotamers
        },
            score_name='ref2015_cart')
        pose_ops.relax_weighted_cart(self.pose,
                                    scorefxn_cart,
                                    chi=self.iso_selector.apply(self.pose),
                                    bb=False, jump=False)
        return self.pose

    def assess_catalyst(self, idx1, resn='GLU'):
        assert resn.upper() in ('GLU', 'ASP'), 'Only Glu and Asp are allowed'
        MutateResidue(target=idx1, new_res=resn).apply(self.pose)
        # constrain the carboxylate to the O of the isopeptide
        ...
        raise NotImplementedError('This is not implemented yet')

