# isopeptide-builder
Given a protein without an isopeptide predict where one would work best.

## Rationale
Isopeptide bonds are a type of covalent bond that can form between the side chains of lysine and aspartate/glutamate residues.
This irreversible bond, were it engineered into a protein, could be used to stabilise the protein structure,
and could be used to **create a binding protein that once bound and crosslinked to its target could not be removed**.

## proof of concept / abandoned project
> This is a proof of concept

This currently works by molecular mechanics (physics),
but with a neural network operating from the onset in 3D space and at an atom-level of detail it would be much faster.
This would be a denoising diffusion network or a flow-based model.
Unfortunately, the over-arching project was not funded.

## Scan

```python
import itertools
import pyrosetta
...
prc = pyrosetta.rosetta.core
from isopeptide_scan import IsopetideScanner

# get only core and boundary residues
selector = prc.select.residue_selector.LayerSelector()
selector.set_layers(pick_core=True, pick_boundary=True, pick_surface=False)
idxs = prc.select.residue_selector.ResidueVector(selector.apply(pose))

scanner = IsopetideScanner(pose)
df = scanner(idxs)
df.loc[df.verdict == 'favorable']
```

## Isopeptide not transition state

The transition state of an isopeptide bond is a hemiacetal and is chiral,
that means two versions would need to be tested when scanning.
This while the order is

1. distance
2. ideal isopeptide
3. mutant w/ isopeptide
4. mutant w/ isopeptide and catalytic residue
5. mutant w/ transition state(s) and catalytic residue

## Reference isopeptide

The first step given two residues in a pose is to check an isopeptide would work.
`IsopetideScanner(..).is_suitable(res1, res2)` will do a crude distance check,
while `IsopetideScanner(..).score_isopeptide_plausibility` will check
that the backbones of the residues would work by playing with an ideal isopeptide.
This underhood is made with `create_ideal_isopeptide`.

```python
scanner = IsopetideScanner(pyrosetta.Pose())
ideal: pyrosetta.Pose = scanner.create_ideal_isopeptide(7.4)
print(scanner.scorefxn_cart_ex(ideal))  # ``scorefxn_cart_ex`` has constraints but no fa_sol and fa_dun

import pyrosetta_help as ph
display( ph.pose2pandas(ideal, scanner.scorefxn_cart).T )
display( ph.constraints2pandas(ideal) )
```
Regarding the constraints present (second table), see `constrain_isopeptide` function.
At present this adds four bond lengths, 2 angle and 2 dihedral constraints.
The code is written to be able to support non Asx-lys isopeptides,
but the constraints are for Asx-lys for now.

The constraints for the amide are strict.
They could be tailored to reflex values in the CCDC database.

As a test case, a two residues with a known isopeptide bond gets a negative score (good)
for the potential to have an isopeptide!

```python
filename = '2x5p.pdb'  # Streptococcus pyogenes fibronectin binding protein Fbab-B = spyCatcher parent
ref:pyrosetta.Pose = pyrosetta.pose_from_file(filename)  # noqa
pyrosetta.rosetta.core.pose.remove_nonprotein_residues(ref)
print( ref.sequence() )
asp117 = ref.pdb_rsd(('A', 117))
lys31 = ref.pdb_rsd(('A', 31))
assert asp117.name3() == 'ASP'
assert lys31.name3() == 'LYS'
scanner = IsopetideScanner(ref)
print('isopeptide residues are correct... ', scanner.is_suitable(lys31, asp117))
d: float = lys31.xyz('CA').distance( asp117.xyz('CA') )
print('residue are at', d, 'Å')
print('score', scanner.score_isopeptide_plausibility(lys31, asp117))
```

## Challenges

I have written previously [how to encode an isopeptide and isopetide hemiacetal transition state](https://blog.matteoferla.com/2018/09/everything-you-wanted-to-know-about.html),
one omission is that the energy of the residues as calculated in PyRosetta are rather high.

### Validity

Scanning combinatorially is a challenge.
The distance and dihedral angle between the two C&alpha;-C&beta; bonds would be a good start.

> Angle not coded.

### Tangle
First, a crystal structure is the product of a diffraction pattern of a crystalline lattice of protein 
that ought to be in identical conformations. When there is an ensemble what one gets from the density is an average,
which can result in a [tangle](https://bl831.als.lbl.gov/~jamesh/challenge/twoconf/),
wherein the solved structure is an energetically unhealthy average of two or more conformations.

For example see ``fa_rep`` (Lenard-Jones repulsion) of the residues in the isopeptide bond below:

```python
import pyrosetta
import pyrosetta_help as ph

logger = ph.configure_logger()
pyrosetta.distributed.maybe_init(extra_options=ph.make_option_string(no_optH=False,
                                                ex1=None,
                                                ex2=None,
                                                #mute='all',
                                                ignore_unrecognized_res=True,
                                                load_PDB_components=False,
                                                ignore_waters=False)
                               )
filename = '2x5p.pdb'  # Streptococcus pyogenes fibronectin binding protein Fbab-B = spyCatcher parent
ref:pyrosetta.Pose = pyrosetta.pose_from_file(filename)  # noqa
pyrosetta.rosetta.core.pose.remove_nonprotein_residues(ref)
print( ref.sequence() )

scorefxn = pyrosetta.get_fa_scorefxn()
ph.pose2pandas(ref, scorefxn).set_index('residue').T[['D117:A', 'K31:A']]
```

|                       |   D117:A  |   K31:A  |
|:----------------------|----------:|---------:|
| fa_atr                |      -6.4 |     -9.1 |
| fa_rep                |     110.9 |     84.5 |
| fa_sol                |       8.3 |      9.3 |
| fa_intra_rep          |       1.4 |      3.7 |
| fa_intra_sol_xover4   |       0.6 |      0.6 |
| lk_ball_wtd           |       0.4 |     -0.3 |
| fa_elec               |      -0.7 |     -2.6 |
| hbond_sr_bb           |       0   |      0   |
| hbond_lr_bb           |      -0.5 |     -0.6 |
| hbond_bb_sc           |       0   |      0   |
| hbond_sc              |       0   |      0   |
| dslf_fa13             |       0   |      0   |
| atom_pair_constraint  |       0   |      0   |
| coordinate_constraint |       0   |      0   |
| angle_constraint      |       0   |      0   |
| dihedral_constraint   |       0   |      0   |
| omega                 |       0.3 |      0.1 |
| fa_dun                |       4.7 |     14.6 |
| p_aa_pp               |       0.4 |      0.6 |
| yhh_planarity         |       0   |      0   |
| ref                   |      -2.1 |     -0.7 |
| rama_prepro           |      -0.3 |      0.3 |
| cart_bonded           |       0.6 |      2.4 |
| total_score           |      64.4 |     55.1 |

The Dunbrack potential is a statistical term based on how frequently the chi angles go that way:
isopeptides are uncommon so it is not surprising that the Dunbrack potential is high.

## RDKit

> This snippet is here for now...

What does RDKit MMFF94 think of the isopeptide bond?

```python
from rdkit import Chem
from typing import Tuple
from rdkit.Chem import AllChem

#isopeptide = Chem.MolFromSmiles("C([C@@H](C(=O))N)C(=O)NCCCC[C@@H](C(=O))N")
#isopeptide: Chem.Mol = Chem.MolFromMolFile('/Users/user/iso.mol')
isopeptide: Chem.Mol = Chem.MolFromPDBBlock(ph.get_pdbstr(test))
    
bad_amide = Chem.MolFromSmiles('C(-O)-N')
if isopeptide.HasSubstructMatch(bad_amide):
    # PDB bond order is wrong
    c, o, n = isopeptide.GetSubstructMatch(bad_amide)
    isopeptide.GetBondBetweenAtoms(c, o).SetBondType(Chem.BondType.DOUBLE)
    isopeptide.GetAtomWithIdx(o).SetNumExplicitHs(0)
    isopeptide.GetAtomWithIdx(n).SetNumExplicitHs(1)
    isopeptide.GetAtomWithIdx(o).SetFormalCharge(0)
    isopeptide.GetAtomWithIdx(n).SetFormalCharge(0)

AllChem.SanitizeMol(isopeptide)
isopeptide = AllChem.AddHs(isopeptide, addCoords=True)
isopeptide.RemoveAllConformers()
AllChem.EmbedMolecule(isopeptide)
p = AllChem.MMFFGetMoleculeProperties(isopeptide, 'MMFF94')
ff = AllChem.MMFFGetMoleculeForceField(isopeptide, p, ignoreInterfragInteractions=False)
ff.Initialize()
ff.Minimize()
ff.CalcEnergy()
```