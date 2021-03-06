import MDAnalysis as mda
import numpy as np
import STR2GMX as str2gmx

u = mda.Universe('step1.gro')
pbc = u.dimensions[0:3]
u.atoms.positions += pbc / 2
u.dimensions = [pbc[0]*2, pbc[1]*2, pbc[2]*2, 90, 90, 90]

s = str2gmx.Solvate(u)
newu = s.run()
newu.atoms.write('step2.gro')


