import types
from typing import (Callable, Tuple, Dict)
import pyrosetta
from .common import IsopetideType, create_weighted_scorefxn

prc: types.ModuleType = pyrosetta.rosetta.core
Residue: types.ModuleType = prc.conformation.Residue


# ------------------------------------------------------------------------------------------------

class PoseOps:
    def __init__(self, pose: pyrosetta.Pose):
        self.pose = pose

    def constrain_isopeptide(self, asx_i: int = 1, lys_i: int = 2):
        """
        Add a dihedral constraint to the pose ``self.pose``
        for the isopeptide bond between the residues with pose indices ``asx_i`` and ``lys_i``

        Does not check for pre-existing so ``pose.remove_constraints()`` may be required...
        """
        CircularHarmonicFunc: Callable = prc.scoring.func.CircularHarmonicFunc  # noqa
        # ## Declare the atoms for constraints
        D_fore_atom = pyrosetta.rosetta.core.id.AtomID(atomno_in=self.pose.residue(asx_i).atom_index('CG'), rsd_in=asx_i)
        D_aft_atom = pyrosetta.rosetta.core.id.AtomID(atomno_in=self.pose.residue(asx_i).atom_index('CB'), rsd_in=asx_i)
        # sidechain conjugated has no OD1 but a OD:
        D_side_atom = pyrosetta.rosetta.core.id.AtomID(atomno_in=self.pose.residue(asx_i).atom_index('OD'), rsd_in=asx_i)
        K_fore_atom = pyrosetta.rosetta.core.id.AtomID(atomno_in=self.pose.residue(lys_i).atom_index('NZ'), rsd_in=lys_i)
        K_aft_atom = pyrosetta.rosetta.core.id.AtomID(atomno_in=self.pose.residue(lys_i).atom_index('CE'), rsd_in=lys_i)
        # sidechain conjugated has no 1HZ but a HZ:
        K_side_atom = pyrosetta.rosetta.core.id.AtomID(atomno_in=self.pose.residue(lys_i).atom_index('HZ'), rsd_in=lys_i)
        fun = CircularHarmonicFunc(3.14159, 0.2)
        # ## Declare constraints and add
        # get the functions
        DihedralConstraint = pyrosetta.rosetta.core.scoring.constraints.DihedralConstraint
        # Dihedral CE 9A NZ 9A CG 121A CB 121A CIRCULARHARMONIC 3.14159 0.2
        con = DihedralConstraint(K_aft_atom, K_fore_atom, D_fore_atom, D_aft_atom, fun)
        self.pose.add_constraint(con)
        # Dihedral 1HZ 9A NZ 9A CG 121A OD 121A CIRCULARHARMONIC 3.14159 0.2
        con = DihedralConstraint(K_side_atom, K_fore_atom, D_fore_atom, D_side_atom, fun)
        self.pose.add_constraint(con)


    def remove_termini(self, indices=(1, 2)):
        rm_upper = prc.conformation.remove_upper_terminus_type_from_conformation_residue
        rm_lower = prc.conformation.remove_lower_terminus_type_from_conformation_residue
        for i in indices:
            rm_upper(self.pose.conformation(), i)
            rm_lower(self.pose.conformation(), i)


    def constrain_CA_distance(self, d: float, indices: Tuple[int, int] = (1, 2), atom_name: str = 'CA'):
        """
        Constrain the CA atoms of the residues with indices ``indices`` by a distance of ``d``.
        This is used for the test pose basically to simulate a protein conformation.
        """
        A_atom = prc.id.AtomID(atomno_in=self.pose.residue(indices[0]).atom_index(atom_name), rsd_in=indices[0])
        B_atom = prc.id.AtomID(atomno_in=self.pose.residue(indices[1]).atom_index(atom_name), rsd_in=indices[1])
        DihedralConstraint = prc.scoring.constraints.AtomPairConstraint
        HarmonicFunc = prc.scoring.func.HarmonicFunc
        fun = HarmonicFunc(d, 0.2)
        con = DihedralConstraint(A_atom, B_atom, fun)
        self.pose.add_constraint(con)

    def remove_constraints(self):
        """
        This is silly, but stops confusion...
        """
        self.pose.remove_constraints()


    def relax_weighted_cart(self, scorefxn):
        movemap = pyrosetta.MoveMap()
        movemap.set_bb(True)
        movemap.set_chi(True)
        movemap.set_jump(True)
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 5)  # noqa
        # relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
        relax.set_movemap(movemap)
        relax.minimize_bond_angles(True)
        relax.minimize_bond_lengths(True)
        relax.cartesian(True)
        relax.apply(self.pose)


    def relax_sidechains(self, scorefxn, cycles=3):
        """
        Relax only sidechains
        """
        if cycles == 0:
            return
        movemap = pyrosetta.MoveMap()
        movemap.set_bb(False)
        movemap.set_chi(True)
        movemap.set_jump(False)
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)  # noqa
        # relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
        relax.set_movemap(movemap)
        relax.apply(self.pose)
