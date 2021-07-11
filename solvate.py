import MDAnalysis as mda
import numpy as np
import os

LL = 18.856
NN = 216

class Solvate:

    def __init__(self, u):
        self.pbc = u.dimensions[0:3]
        if np.sum(self.pbc) == 0:
            raise ValueError('pbc defined?')

        self.u = u

        # Read tip216.crd coordinates
        pos = []
        path = os.path.dirname(os.path.realpath(__file__))
        with open(path + '/FF/C36/toppar/tip216.crd') as f:
            for line in f.readlines():
                if line.startswith('*'): continue
                sl = line.split()
                if len(sl) == 7:
                    x  = float(sl[4])
                    y  = float(sl[5])
                    z  = float(sl[6])
                    pos.append([x, y, z])
        pos = np.array(pos)
        pos += np.array([LL/2, LL/2, LL/2])
        self.pos = pos


    def run(self, zUP = None, zDW = None):
        
        # Existing water
        ew = self.u.select_atoms('resname TIP3')
        if ew.n_atoms == 0:
            maxresid = 0
        else:
            maxresid = max(ew.resids)

        Nx = int(self.pbc[0] / LL) + 1
        Ny = int(self.pbc[1] / LL) + 1
        Nz = int(self.pbc[2] / LL) + 1
        
        # can you make this more efficient?
        # the dimension of 216 water box is pretty small,
        # so you have to do this a lot of times...
        water_pos = []
        for i in range(Nx):
            for j in range(Ny):
                for k in range(Nz):
                    dr = np.array([LL * i, LL * j, LL * k])
                    water_pos.append(self.pos + dr)


        water_pos = np.concatenate(water_pos, axis=0)
        bAx = water_pos[:,0][::3] > self.pbc[0]
        bAy = water_pos[:,1][::3] > self.pbc[1]
        
        wpz = water_pos[:,2][::3]
        if zUP == None and zDW == None:
            bAz = wpz > self.pbc[2]
        elif zDW == None:
            bAz = (wpz > self.pbc[2]) | (wpz < zUP)
        elif zUP == None:
            bAz = (wpz > self.pbc[2]) | (wpz > zDW)
        else:
            bAz = (wpz > self.pbc[2]) | ((zDW < wpz) & (wpz < zUP))

        bA  = bAx | bAy | bAz
        waterbox_pos = water_pos[np.repeat(~bA, 3)]

        n_atoms = len(waterbox_pos)
        n_res   = int(n_atoms / 3)

        sol = mda.Universe.empty(n_atoms = n_atoms,
                                 n_residues = n_res,
                                 atom_resindex = np.repeat(np.arange(n_res), 3),
                                 residue_segindex = [0] * n_res,
                                 trajectory = True)

        sol.add_TopologyAttr('resnames', ['TIP3'] * n_res)
        sol.add_TopologyAttr('resids', np.arange(1, n_res + 1))
        sol.add_TopologyAttr('names', ['OH2', 'H1', 'H2'] * n_res)
        sol.atoms.positions = waterbox_pos

        newu = mda.Merge(self.u.atoms, sol.atoms)
        newu.dimensions = self.u.dimensions

        newu.add_TopologyAttr('segids', ['ORI', 'NEW'])
        sel = '(segid ORI) or (segid NEW and not byres (name OH2 and around 5 (segid ORI)))'
        ag = newu.select_atoms(sel)
        
        newu2 = mda.Merge(ag)
        newu2.dimensions = self.u.dimensions

        newWater = newu2.select_atoms('segid NEW and resname TIP3')
        newWater.residues.resids = np.arange(maxresid + 1, maxresid + newWater.n_residues + 1)
        assert newWater.n_atoms == newWater.n_residues * 3, 'atoms missing?'
        
        allWater = newu2.select_atoms('resname TIP3')
        assert allWater.n_atoms == allWater.n_residues * 3, 'atoms missing?'

        return newu2


