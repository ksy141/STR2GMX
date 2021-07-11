import MDAnalysis as mda
import STR2GMX    as str2gmx

u = mda.Universe('../C/step5_input.gro')
u.atoms.write('step5_input.gro')

toppar = str2gmx.ReadToppars(['../../../FF/C36/toppar/top_all36_lipid.rtf',
                              '../../../FF/C36/toppar/par_all36_lipid.prm',
                              '../../../FF/C36/toppar/toppar_water_ions.str'])

mol_POPC = str2gmx.Molecule(u.select_atoms('resname POPC'), toppar)
mol_POT  = str2gmx.Molecule(u.select_atoms('resname POT'),  toppar)
mol_CLA  = str2gmx.Molecule(u.select_atoms('resname CLA'),  toppar)
mol_TIP3 = str2gmx.Molecule(u.select_atoms('resname TIP3'), toppar)

mols = str2gmx.Molecules([mol_POPC, mol_POT, mol_CLA, mol_TIP3], toppar)
mols.write()

with mda.selections.gromacs.SelectionWriter('index.ndx', mode='w') as ndx:
    ndx.write(u.select_atoms('resname POPC'), name='MEMB')
    ndx.write(u.select_atoms('resname POT CLA TIP3'), name='SOLV')
   
