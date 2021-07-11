import MDAnalysis as mda
import STR2GMX    as str2gmx

u = mda.Universe('../C36/step5_input.gro')

POPC = u.select_atoms('resname POPC')

### The order of atoms has changed from C36 to C36_LJPME.
### C36:       N, C12, H12A, H12B, C13, H13A, H13B, H13C, C14, C14A, C14B, C14C, C15, H15A, H15B, H15C
###               C11, H11A, H11B, P,   O13,  O14,  O12,  O11
### C36_LJPME: N, C13, H13A, H13B, H13C, C14, H14A, H14B, H14C, C15, H15A, H15B, H15C, C12, H12A, H12B
###               C11, H11A, H11B, P,   O13,  O14,  O11,  O12

hg = ['N', 'C13', 'H13A', 'H13B', 'H13C', 'C14', 'H14A', 'H14B', 
      'H14C', 'C15', 'H15A', 'H15B', 'H15C', 'C12', 'H12A', 'H12B',
      'C11', 'H11A', 'H11B', 'P', 'O13', 'O14', 'O11', 'O12']

ags = []
for res in POPC.residues:
    ag  = res.atoms
    ag2 = ag.select_atoms('name N') + \
          ag.select_atoms('name C13 H13A H13B H13C C14 H14A H14B H14C C15 H15A H15B H15C') + \
          ag.select_atoms('name C12 H12A H12B') + \
          ag.select_atoms('name C11 H11A H11B P O13 O14') + \
          ag.select_atoms('name O11') + \
          ag.select_atoms('name O12') + \
          ag.select_atoms('not name ' + ' '.join(hg))
    ags.append(ag2)

ags.append(u.select_atoms('resname POT CLA TIP3'))
newu = mda.Merge(*ags)
newu.dimensions = u.dimensions
newu.atoms.write('step5_input.gro')

toppar = str2gmx.ReadToppars(['../../../FF/C36_LJPME/top_all36_lipid_ljpme.rtf',
                              '../../../FF/C36_LJPME/par_all36_lipid_ljpme.prm',
                              '../../../FF/C36_LJPME/main/toppar_water_ions.str'])

mol_POPC = str2gmx.Molecule(newu.select_atoms('resname POPC'), toppar)
mol_POT  = str2gmx.Molecule(newu.select_atoms('resname POT'),  toppar)
mol_CLA  = str2gmx.Molecule(newu.select_atoms('resname CLA'),  toppar)
mol_TIP3 = str2gmx.Molecule(newu.select_atoms('resname TIP3'), toppar)

mols = str2gmx.Molecules([mol_POPC, mol_POT, mol_CLA, mol_TIP3], toppar)
mols.write()


with mda.selections.gromacs.SelectionWriter('index.ndx', mode='w') as ndx:
    ndx.write(newu.select_atoms('resname POPC'), name='MEMB')
    ndx.write(newu.select_atoms('resname POT CLA TIP3'), name='SOLV')
 
