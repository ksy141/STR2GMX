import numpy  as np
import pandas as pd
import MDAnalysis as mda
import os

class Molecule:

    def __init__(self, ag, toppar, generate_angles = True, generate_dihedrals = True):
        assert ag.n_atoms != 0, 'AtomGroup contains no atoms.'
        
        self.system = ag.universe
        newag = ag.residues[0].atoms
        self.resname = newag.resnames[0]
        self.u = mda.Merge(newag)
        self.RESIDUE = toppar.RESI[self.resname]
        self.ATOMS   = toppar.ATOMS

        name2index = {}
        index2name = {}
        for i, atom in enumerate(self.u.atoms):
            name2index[atom.name] = i
            index2name[i] = atom.name

        self.name2index = name2index
        self.index2name = index2name
        self.generate_angles = generate_angles
        self.generate_dihedrals = generate_dihedrals

        self.generate()

    def generate(self):
        self.ag = self.u.atoms
        resname = self.resname

        assert np.all(self.ag.names == np.array(self.RESIDUE['names'])), 'are atoms missing in ' + resname + '?'

        self.ag.types = self.RESIDUE['types']
        self.u.add_TopologyAttr('charges', self.RESIDUE['charges'])
        
        for atom in self.ag:
            atom.mass = self.ATOMS[atom.type]['mass']
        
        self.bsorted = []
        self.asorted = []
        self.dsorted = []
        self.isorted = []
        
        ### BONDS
        if len(self.RESIDUE['bonds']) > 0:
            bonds_str = np.array(self.RESIDUE['bonds']).reshape(-1)
            bonds = np.array([self.name2index[b] for b in bonds_str]).reshape(-1, 2)

            bb = []
            for i, j in bonds:
                if i > j:
                    i, j = j, i
                bb.append([i, j])

            df = pd.DataFrame(bb)
            self.bsorted = df.sort_values(by=[0,1]).to_numpy()
            self.u.add_TopologyAttr('bonds', self.bsorted)
        
        else:
            self.u.add_TopologyAttr('bonds', [])

       

        ### IMPROPERS
        if len(self.RESIDUE['imprs']) > 0:
            imprs = self.RESIDUE['imprs']
            for impr in imprs:
                atom1, atom2, atom3, atom4 = impr
                i = self.name2index[atom1]
                j = self.name2index[atom2]
                k = self.name2index[atom3]
                l = self.name2index[atom4]
                self.isorted.append([i, j, k, l])
 
        
        ### ANGLES
        angles = np.zeros((0,3), dtype=np.int64)
        if len(self.RESIDUE['angles']) > 0:
            angles_str = np.array(self.RESIDUE['angles']).reshape(-1)
            predefined = np.array([self.name2index[a] for a in angles_str]).reshape(-1, 3)
            angles     = np.concatenate([angles, predefined], axis=0)
        
        if self.generate_angles:
            guessed = np.array(mda.topology.guessers.guess_angles(self.u.bonds)).reshape(-1, 3)
            angles  = np.concatenate([angles, guessed], axis=0)

        if len(angles) > 0:
            df = pd.DataFrame(angles)
            self.asorted = df.sort_values(by=[1,0,2]).to_numpy()
            self.u.add_TopologyAttr('angles', self.asorted)
        else:
            self.u.add_TopologyAttr('angles', [])
        

        ### DIHEDRALS
        dihedrals = np.zeros((0, 4), dtype=np.int64)
        if len(self.RESIDUE['dihedrals']) > 0:
            dihe_str   = np.array(self.RESIDUE['dihedrals']).reshape(-1)
            predefined = np.array([self.name2index[d] for d in dihe_str]).reshape(-1, 4)
            dihedrals = np.concatenate([dihedrals, predefined], axis=0)
        
        if self.generate_dihedrals:
            guessed   = np.array(mda.topology.guessers.guess_dihedrals(self.u.angles)).reshape(-1, 4)
            dihedrals = np.concatenate([dihedrals, guessed], axis=0)

        if len(dihedrals) > 0:
            df = pd.DataFrame(dihedrals)
            self.dsorted = df.sort_values(by=[1,2,0,3]).to_numpy()
            self.u.add_TopologyAttr('dihedrals', self.dsorted)


        ### PAIRS1-4
        ### should be followed by obtaining self.dsorted
        ### because if len(self.dsorted) == 0, then it will return []
        self.psorted = self._guess_pairs14(self.bsorted)


    def _guess_pairs14(self, bonds):
        
        # build 1-4 pairs
        pairs14 = []
        if len(self.dsorted) == 0:
            return pairs14

        for b in bonds:
            for i in [0, 1]:
                idx1 = b[i]
                bA1l = bonds[:,0] == idx1
                bA1r = bonds[:,1] == idx1

                ids2 = np.concatenate([bonds[:,1][bA1l], bonds[:,0][bA1r]])
                if len(ids2) == 0: continue
                
                bA2l = np.isin(bonds[:,0], ids2)
                bA2r = np.isin(bonds[:,1], ids2)

                ids3 = np.concatenate([bonds[:,1][bA2l], bonds[:,0][bA2r]])
                ids3 = np.delete(ids3, np.isin(ids3, idx1))
                if len(ids3) == 0: continue
                
                bA3l = np.isin(bonds[:,0], ids3)
                bA3r = np.isin(bonds[:,1], ids3)

                ids4 = np.concatenate([bonds[:,1][bA3l], bonds[:,0][bA3r]])
                ids4 = np.delete(ids4, np.isin(ids4, ids2))
                ids4 = np.delete(ids4, np.isin(ids4, ids3)) #cholesterol pentagon
                if len(ids4) == 0: continue

                for idx4 in ids4:
                    if idx1 > idx4:
                        result = [idx4, idx1]

                    else:
                        result = [idx1, idx4]

                    if result not in pairs14:
                        pairs14.append(result)
 
        df = pd.DataFrame(pairs14)
        return df.sort_values(by=[0,1]).to_numpy()

