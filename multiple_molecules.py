import pandas as pd
from   .molecule import Molecule

class MultipleMolecules:

    def __init__(self):
        pass

    def generate(self, AtomicGroup, toppar, 
                 generate_angles = True, 
                 generate_dihedrals = True):
        
        resnames = pd.unique(AtomicGroup.residues.resnames)
        
        mols = []
        for resname in resnames:
            #print(resname)
            ag = AtomicGroup.select_atoms('resname %s' %resname)
            if ag.n_atoms == 0: continue
            mol = Molecule(ag, toppar, generate_angles, generate_dihedrals)
            mols.append(mol)

        return mols

