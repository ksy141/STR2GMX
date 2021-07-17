import MDAnalysis as mda
import numpy as np
import STR2GMX as str2gmx

u = mda.Universe('step1.gro')
pbc = u.dimensions[0:3]
u.atoms.positions += pbc / 2
u.dimensions = \
[pbc[0]*2, pbc[1]*2, pbc[2]*2, 
90, 90, 90]

s = str2gmx.Solvate(u)
newu = s.run(cutoff=5.0)
newu.atoms.write('step2.gro')

ag1 = newu.select_atoms('byres (resname TRIO and prop x > %f)' %(pbc[0]/2))
ag2 = newu.select_atoms('resname TIP3 and prop x > %f' %(pbc[0]/2 + 15))
mda.Merge(ag1, ag2).atoms.write('vis.gro')

