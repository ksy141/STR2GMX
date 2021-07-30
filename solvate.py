import MDAnalysis as mda
import numpy as np
import os
import random
from   MDAnalysis.analysis.distances import distance_array

LL = 18.856
NN = 216

class Solvate:

    def __init__(self):
        pass

    def sol(self, u, zUP = None, zDW = None, cutoff=5):
        pbc = u.dimensions[0:3]
        if np.sum(pbc) == 0:
            raise ValueError('pbc defined?')

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

        
        # Existing water
        ew = u.select_atoms('resname TIP3')
        if ew.n_atoms == 0:
            maxresid = 0
        else:
            maxresid = max(ew.resids)

        Nx = int(pbc[0] / LL) + 1
        Ny = int(pbc[1] / LL) + 1
        Nz = int(pbc[2] / LL) + 1
        
        # can you make this more efficient?
        # the dimension of 216 water box is pretty small,
        # so you have to do this a lot of times...
        water_pos = []
        for i in range(Nx):
            for j in range(Ny):
                for k in range(Nz):
                    dr = np.array([LL * i, LL * j, LL * k])
                    water_pos.append(self.pos + dr)

        dxdydz = 1
        water_pos = np.concatenate(water_pos, axis=0)
        bAx = water_pos[:,0][::3] > pbc[0] - dxdydz
        bAy = water_pos[:,1][::3] > pbc[1] - dxdydz
        
        wpz = water_pos[:,2][::3]
        if zUP == None and zDW == None:
            bAz =  wpz > pbc[2] - dxdydz
        elif zDW == None:
            bAz = (wpz > pbc[2] - dxdydz) | (wpz < zUP)
        elif zUP == None:
            bAz = (wpz > pbc[2] - dxdydz) | (wpz > zDW)
        else:
            bAz = (wpz > pbc[2] - dxdydz) | ((zDW < wpz) & (wpz < zUP))

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

        newu = mda.Merge(u.atoms, sol.atoms)
        newu.dimensions = u.dimensions

        newu.add_TopologyAttr('segids', ['ORI', 'NEW'])
        sel = '(segid ORI) or (segid NEW and not byres (name OH2 and around %d (segid ORI)))' %cutoff
        ag = newu.select_atoms(sel)
        
        newu2 = mda.Merge(ag)
        newu2.dimensions = u.dimensions

        newWater = newu2.select_atoms('segid NEW and resname TIP3')
        newWater.residues.resids = np.arange(maxresid + 1, maxresid + newWater.n_residues + 1)
        assert newWater.n_atoms == newWater.n_residues * 3, 'atoms missing?'
        
        allWater = newu2.select_atoms('resname TIP3')
        assert allWater.n_atoms == allWater.n_residues * 3, 'atoms missing?'

        return newu2



    def ion(self, ag, qtot, conc = 0.15,  pos = 'POT', neg = 'CLA'):
        print("##### ION #####")
        ### ag = u.select_atoms('resname TIP3')
        u = ag.universe
        n_res = ag.n_residues

        if pos in ['MG', 'CAL', 'BAR', 'ZN2', 'CD2']:
            factor = 2
        else:
            factor = 1

        # Existing Positive Ions of the same type
        ep = u.select_atoms('resname %s' %pos)
        if ep.n_atoms == 0:
            max_ep_resid = 0
        else:
            max_ep_resid = max(ep.resids)

        # Existing Negative Ions of the same type
        en = u.select_atoms('resname %s' %neg)
        if en.n_atoms == 0:
            max_en_resid = 0
        else:
            max_en_resid = max(en.resids)
        
        pos_add = int(conc * (n_res * 18 / 1000))
        neg_add = int(pos_add * factor)

        if qtot < 0:
            pos_add += int(- qtot / factor)

        qfinal = qtot + pos_add * factor - neg_add
        
        if qfinal < 0:
            pos_add += 1

        qfinal = qtot + pos_add * factor - neg_add
        
        if qfinal > 0:
            neg_add += qfinal

        qfinal = qtot + pos_add * factor - neg_add
        
        res0 = ag.residues[ random.sample(range(n_res), pos_add + neg_add) ]
        assert res0.n_residues == pos_add + neg_add, 'residue?'
        resp = res0[:pos_add]
        resn = res0[pos_add:]

        assert resp.atoms.n_atoms == 3 * resp.n_residues, 'water?'
        assert resn.atoms.n_atoms == 3 * resn.n_residues, 'water?'

        resp.atoms[::3].names = pos
        resp.resnames = pos

        resn.atoms[::3].names = neg
        resn.resnames = neg
        
        newPOS = u.select_atoms('resname %s and name H1 H2' %pos)
        newPOS.residues.resids = np.arange(max_ep_resid + 1, max_ep_resid + pos_add + 1)

        newNEG = u.select_atoms('resname %s and name H1 H2' %neg)
        newNEG.residues.resids = np.arange(max_en_resid + 1, max_en_resid + neg_add + 1)

        newu = mda.Merge(u.select_atoms('not (resname %s %s and name H1 H2)' %(pos, neg)))
        
        ag1 = newu.select_atoms('not (resname TIP3 %s %s)' %(pos, neg))
        ag2 = newu.select_atoms('resname %s' %pos)
        ag3 = newu.select_atoms('resname %s' %neg)
        ag4 = newu.select_atoms('resname TIP3')
        
        newconc = min([ag2.n_atoms, ag3.n_atoms]) / (ag4.n_residues * 18 / 1000)
        print('CONC: %.3f M' %newconc)

        newu2 = mda.Merge(ag1, ag2, ag3, ag4)
        newu2.dimensions = u.dimensions
         
        newPOS = newu2.select_atoms('resname %s' %pos)
        newNEG = newu2.select_atoms('resname %s' %neg)
        
        if newPOS.n_atoms * newNEG.n_atoms > 0:
            d1 = distance_array(newPOS.positions, newNEG.positions, box=newu2.dimensions)
            d2 = distance_array(newPOS.positions, newPOS.positions, box=newu2.dimensions)
            d3 = distance_array(newNEG.positions, newNEG.positions, box=newu2.dimensions)

            dshort = min([np.min(d1), min(d2[d2 != 0]), min(d3[d3 != 0])])
            print("THE SHORTEST ION PAIR DISTANCE IS %.3f A" %dshort)
        
        return newu2

