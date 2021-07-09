import MDAnalysis as mda
import STR2GMX    as str2gmx

u = mda.Universe('../C36.charmm-gui.gromacs/step5_input.gro')

toppar = str2gmx.ReadToppars('../../../FF/C36/toppar.str', verbose=False)
POPC   = str2gmx.Molecule(u.select_atoms('resname POPC'), toppar)
TIP3   = str2gmx.Molecule(u.select_atoms('resname TIP3'), toppar, generate_angles=False, generate_dihedrals=False)

mols = str2gmx.Molecules([POPC, TIP3], toppar)
mols.write()

