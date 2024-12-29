import types
from typing import (Callable, Tuple, Dict)
import pyrosetta
from pyrosetta.rosetta.core.scoring.constraints import AtomPairConstraint

from .common import IsopetideType, create_weighted_scorefxn

prc: types.ModuleType = pyrosetta.rosetta.core
Residue: types.ModuleType = prc.conformation.Residue

# ------------------------------------------------------------------------------------------------


def constrain_isopeptide(pose, asx_i: int = 1, lys_i: int = 2,
                         dist_sd=0.3,
                         angle_sd_radians=0.2,
                         dihedral_sd_radians=0.2):
    """
    Add a dihedral constraint to the pose ``pose``
    for the isopeptide bond between the residues with pose indices ``asx_i`` and ``lys_i``

    Does not check for pre-existing so ``pose.remove_constraints()`` may be required...
    """
    # ## Declare the atoms for constraints
    D_fore_atom = pyrosetta.rosetta.core.id.AtomID(atomno_in=pose.residue(asx_i).atom_index('CG'), rsd_in=asx_i)
    D_aft_atom = pyrosetta.rosetta.core.id.AtomID(atomno_in=pose.residue(asx_i).atom_index('CB'), rsd_in=asx_i)
    # sidechain conjugated has no OD1 but a OD:
    D_side_atom = pyrosetta.rosetta.core.id.AtomID(atomno_in=pose.residue(asx_i).atom_index('OD'), rsd_in=asx_i)
    K_fore_atom = pyrosetta.rosetta.core.id.AtomID(atomno_in=pose.residue(lys_i).atom_index('NZ'), rsd_in=lys_i)
    K_aft_atom = pyrosetta.rosetta.core.id.AtomID(atomno_in=pose.residue(lys_i).atom_index('CE'), rsd_in=lys_i)
    # sidechain conjugated has no 1HZ but a HZ:
    K_side_atom = pyrosetta.rosetta.core.id.AtomID(atomno_in=pose.residue(lys_i).atom_index('HZ'), rsd_in=lys_i)
    # ## Declare constraints and add
    # ### Distance
    HarmonicFunc: Callable = prc.scoring.func.HarmonicFunc
    AtomPairConstraint: Callable = prc.scoring.constraints.AtomPairConstraint
    # AtomPair NZ 9A CG 121A HARMONIC 1.30 0.3
    fun = HarmonicFunc(1.30, dist_sd)
    con = AtomPairConstraint(K_fore_atom, D_fore_atom, fun)
    pose.add_constraint(con)
    # AtomPair CE 9A CG 121A HARMONIC 2.40 0.3
    fun = HarmonicFunc(2.40, dist_sd)
    con = AtomPairConstraint(K_aft_atom, D_fore_atom, fun)
    pose.add_constraint(con)
    # AtomPair 1HZ 9A OD 121A HARMONIC 3.20 0.3
    fun = HarmonicFunc(3.20, dist_sd)
    con = AtomPairConstraint(K_side_atom, D_side_atom, fun)
    pose.add_constraint(con)
    # ### Angle
    AngleConstraint: Callable = prc.scoring.constraints.AngleConstraint
    CircularHarmonicFunc: Callable = prc.scoring.func.CircularHarmonicFunc  # noqa
    # Angle CE 9A NZ 9A CG 121A CIRCULARHARMONIC 2.08 0.2 #119 deg
    fun = CircularHarmonicFunc(2.08, angle_sd_radians)
    con = AngleConstraint(K_aft_atom, K_fore_atom, D_fore_atom, fun)
    pose.add_constraint(con)
    # Angle NZ 9A CG 121A OD 121A CIRCULARHARMONIC 2.08 0.2
    con = AngleConstraint(K_fore_atom,D_fore_atom,D_side_atom, fun)
    pose.add_constraint(con)
    # ### Dihedral
    DihedralConstraint = prc.scoring.constraints.DihedralConstraint
    fun = CircularHarmonicFunc(3.14159, dihedral_sd_radians)
    # Dihedral CE 9A NZ 9A CG 121A CB 121A CIRCULARHARMONIC 3.14159 0.2
    con = DihedralConstraint(K_aft_atom, K_fore_atom, D_fore_atom, D_aft_atom, fun)
    pose.add_constraint(con)
    # Dihedral 1HZ 9A NZ 9A CG 121A OD 121A CIRCULARHARMONIC 3.14159 0.2
    con = DihedralConstraint(K_side_atom, K_fore_atom, D_fore_atom, D_side_atom, fun)
    pose.add_constraint(con)


def remove_termini(pose, indices=(1, 2)):
    rm_upper = prc.conformation.remove_upper_terminus_type_from_conformation_residue
    rm_lower = prc.conformation.remove_lower_terminus_type_from_conformation_residue
    for i in indices:
        rm_upper(pose.conformation(), i)
        rm_lower(pose.conformation(), i)


def constrain_CA_distance(pose, d: float, indices: Tuple[int, int] = (1, 2), atom_name: str = 'CA'):
    """
    Constrain the CA atoms of the residues with indices ``indices`` by a distance of ``d``.
    This is used for the test pose basically to simulate a protein conformation.
    """
    A_atom = prc.id.AtomID(atomno_in=pose.residue(indices[0]).atom_index(atom_name), rsd_in=indices[0])
    B_atom = prc.id.AtomID(atomno_in=pose.residue(indices[1]).atom_index(atom_name), rsd_in=indices[1])
    AtomPairConstraint = prc.scoring.constraints.AtomPairConstraint
    HarmonicFunc = prc.scoring.func.HarmonicFunc
    fun = HarmonicFunc(d, 0.2)
    con = AtomPairConstraint(A_atom, B_atom, fun)
    pose.add_constraint(con)

def remove_constraints(pose):
    """
    This is silly, but stops confusion...
    """
    pose.remove_constraints()


def relax_weighted_cart(pose, scorefxn, cycles=5, chi=True, bb=True, jump=True):
    movemap = pyrosetta.MoveMap()
    movemap.set_bb(bb)
    movemap.set_chi(chi)
    movemap.set_jump(jump)
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)  # noqa
    # relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
    relax.set_movemap(movemap)
    relax.minimize_bond_angles(True)
    relax.minimize_bond_lengths(True)
    relax.cartesian(True)
    relax.apply(pose)


def relax_sidechains(pose, scorefxn, cycles=3):
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
    relax.apply(pose)
