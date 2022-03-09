import MDAnalysis as mda
import numpy as np
import os
import random
from   MDAnalysis.analysis.distances import distance_array

LL    = 7
intra = 2.8

class Solvate_Sirah:

    def __init__(self):
        pass

    def sol(self, u, zUP = None, zDW = None, cutoff=5):
        pbc = u.dimensions[0:3]
        if np.sum(pbc) == 0:
            raise ValueError('pbc defined?')

        # Existing water
        ew = u.select_atoms('resname WT4')
        if ew.n_atoms == 0:
            maxresid = 0
        else:
            maxresid = max(ew.resids)

        Nx = int(pbc[0] / LL)
        Ny = int(pbc[1] / LL)
        Nz = int(pbc[2] / LL)
        NN = Nx * Ny * Nz

        water_pos = []
        for i in range(Nx):
            for j in range(Ny):
                for k in range(Nz):
                    water_pos.append([LL*i, LL*j, LL*k])

        waterbox_pos = np.array(water_pos)
        n_atoms = len(waterbox_pos)

        solv = mda.Universe.empty(n_atoms = n_atoms * 4,
                                  n_residues = n_atoms,
                                  atom_resindex = np.repeat(np.arange(n_atoms), 4),
                                  residue_segindex = [0] * n_atoms,
                                  trajectory = True)

        solv.add_TopologyAttr('resnames', ['W'] * n_atoms)
        solv.add_TopologyAttr('resids', np.arange(maxresid + 1, maxresid + n_atoms + 1))
        solv.add_TopologyAttr('names', ['WN1', 'WN2', 'WP1', 'WP2'] * n_atoms)

        water_positions = np.repeat(waterbox_pos, 4, axis=0)
        first  = [0, 0, intra]
        second = [0, 2 * np.sqrt(2) / 3 * intra, -intra / 3]
        third  = [+ np.sqrt(2) / np.sqrt(3) * intra, - np.sqrt(2) / 3 * intra, - intra / 3]
        fourth = [- np.sqrt(2) / np.sqrt(3) * intra, - np.sqrt(2) / 3 * intra, - intra / 3]
        single = [first, second, third, fourth]

        solv.atoms.positions = water_positions + np.array(single * n_atoms) + np.array([LL, LL, LL])/2
        
        if u.atoms.n_atoms != 0:
            newu = mda.Merge(u.atoms, solv.atoms)
            newu.dimensions = u.dimensions

            newu.add_TopologyAttr('segids', ['ORI', 'NEW'])
            sel = '(segid ORI) or (segid NEW and not byres (resname WT4 and around %d (segid ORI)))' %cutoff
            ag = newu.select_atoms(sel)
            
            newu2 = mda.Merge(ag)
            newu2.dimensions = u.dimensions

            newWater = newu2.select_atoms('segid NEW and resname WT4')
            newWater.residues.resids = np.arange(maxresid + 1, maxresid + newWater.n_residues + 1)
            assert newWater.n_atoms == newWater.n_residues, 'n_atoms != n_residues?'
            
            allWater = newu2.select_atoms('resname WT4')
            assert allWater.n_atoms == allWater.n_residues, 'n_atoms != n_residues?'
            return newu2

        else:
            solv.dimensions = u.dimensions
            return solv



    def ion(self, ag, qtot, conc = 0.15,  pos = 'POT', neg = 'CLA'):
        print("##### ION #####")
        ### ag = u.select_atoms('resname W')
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

        assert resp.atoms.n_atoms == resp.n_residues, 'water?'
        assert resn.atoms.n_atoms == resn.n_residues, 'water?'

        resp.atoms.names = pos
        resp.resnames = pos

        resn.atoms.names = neg
        resn.resnames = neg
        
        newPOS = u.select_atoms('resname %s' %pos)
        newPOS.residues.resids = np.arange(max_ep_resid + 1, max_ep_resid + pos_add + 1)

        newNEG = u.select_atoms('resname %s' %neg)
        newNEG.residues.resids = np.arange(max_en_resid + 1, max_en_resid + neg_add + 1)

        newu = mda.Merge(u.atoms)
        
        ag1 = newu.select_atoms('not (resname W %s %s)' %(pos, neg))
        ag2 = newu.select_atoms('resname %s' %pos)
        ag3 = newu.select_atoms('resname %s' %neg)
        ag4 = newu.select_atoms('resname W')
        
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

