import MDAnalysis as mda
import STR2GMX as str2gmx

u = mda.Universe('step1.gro')
z = u.select_atoms('name P').center_of_mass()[2]
UP = u.select_atoms('name P and prop z > %f' %z).center_of_mass()[2]
LP = u.select_atoms('name P and prop z < %f' %z).center_of_mass()[2]

s = str2gmx.Solvate(u)
newu = s.run(zUP=UP, zDW=LP, cutoff=3.0)
newu.atoms.write('step2.gro')

