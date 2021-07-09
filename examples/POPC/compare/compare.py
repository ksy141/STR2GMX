import pandas as pd
import numpy  as np

def ff2df(ff):
    readprm = False

    bondtypes     = []
    pairtypes     = []
    angletypes    = []
    dihedraltypes = []

    f = open(ff)
    for line in f:
        if not line: continue
        sl = line.split(';')[0]
        if not sl.rstrip(): continue

        ### WHAT TO READ
        if sl.startswith('[ bondtypes ]'):
            readprm = 'bond'
            continue

        if sl.startswith('[ pairtypes ]'):
            readprm = 'pair'
            continue

        if sl.startswith('[ angletypes ]'):
            readprm = 'angle'
            continue

        if sl.startswith('[ dihedraltypes ]'):
            readprm = 'dihedral'
            continue


        ### READ
        if readprm == 'bond':
            t1, t2, *c = sl.split()
            if t1 > t2: t1, t2 = t2, t1
            bondtypes.append([t1, t2, *c])

        if readprm == 'pair':
            t1, t2, *c = sl.split()
            if t1 > t2: t1, t2 = t2, t1
            pairtypes.append([t1, t2, *c])

        if readprm == 'angle':
            t1, t2, t3, *c = sl.split()
            if t1 > t3: t1, t3 = t3, t1
            angletypes.append([t1, t2, t3, *c])

        if readprm == 'dihedral':
            t1, t2, t3, t4, *c = sl.split()
            if t1 > t4:
                t1, t2, t3, t4 = t4, t3, t2, t1

            if t1 == t4 and t2 > t3:
                t2, t3 = t3, t2


            dihedraltypes.append([t1, t2, t3, t4, *c])


    dfb = pd.DataFrame(bondtypes).sort_values(by=[0,1], ignore_index=True)
    dfp = pd.DataFrame(pairtypes).sort_values(by=[0,1], ignore_index=True)
    dfa = pd.DataFrame(angletypes).sort_values(by=[0,1,2], ignore_index=True)
    dfd = pd.DataFrame(dihedraltypes).sort_values(by=[0,1,2,3], ignore_index=True)

    return dfb, dfp, dfa, dfd


sb, sp, sa, sd = ff2df('../C36.str2gmx.gromacs/toppar/forcefield.itp')
cb, cp, ca, cd = ff2df('../C36.charmm-gui.gromacs/toppar/forcefield.itp')

sb.to_csv('sbonds', '\t')
cb.to_csv('cbonds', '\t')

sa.to_csv('sangles', '\t')
ca.to_csv('cangles', '\t')

sp.to_csv('spairs', '\t')
cp.to_csv('cpairs', '\t')

sd.to_csv('sdihedrals', '\t')
cd.to_csv('cdihedrals', '\t')

