import MDAnalysis as mda
import STR2GMX    as str2gmx

u = mda.Universe('../C36.charmm-gui.gromacs/step5_input.gro')

toppar = str2gmx.ReadToppars('../../../FF/C36/toppar.str', verbose=False)
resnames = []
mollist  = []

### MEMB
residues = u.select_atoms('not resname SOD CLA TIP3').residues
for residue in residues:
    if residue.resname in resnames: continue
    print(residue.resname)
    mol = str2gmx.Molecule(residue.atoms, toppar)
    mollist.append(mol)
    resnames.append(residue.resname)

### SOLV
residues = u.select_atoms('resname SOD CLA TIP3').residues
for residue in residues:
    if residue.resname in resnames: continue
    mol = str2gmx.Molecule(residue.atoms, toppar, generate_angles=False, generate_dihedrals=False)
    mollist.append(mol)
    resnames.append(residue.resname)

mols = str2gmx.Molecules(mollist, toppar)
mols.write()

