import MDAnalysis as mda
import numpy as np
import STR2GMX as str2gmx

u = mda.Universe('step2.gro')
dz = 20

x, y, z = u.dimensions[0:3]
u.atoms.positions += np.array([0, 0, dz])
u.dimensions = [x, y, z+dz*2, 90, 90, 90]

s = str2gmx.Solvate(u)
newu = s.run(cutoff=3.0, zUP=z+dz, zDW=dz)
newu.atoms.write('step3.gro')


