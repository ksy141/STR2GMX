import MDAnalysis as mda
import STR2GMX    as str2gmx

u = mda.Universe('step2.gro')
toppar = str2gmx.ReadToppars('/Users/siyoungkim/SMDAnalysis/FF/CHARMM_LJPME_r/toppar.str')

mol_TRIO = str2gmx.Molecule(u.select_atoms('resname TRIO'), toppar)
mol_TIP3 = str2gmx.Molecule(u.select_atoms('resname TIP3'), toppar)

mols = str2gmx.Molecules([mol_TRIO, mol_TIP3], toppar)
mols.write()

