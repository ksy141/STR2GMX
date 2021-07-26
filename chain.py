import MDAnalysis as mda
import numpy as np

class Chain:

    def __init__(self, ag, toppar):

        assert ag.select_atoms('protein').n_atoms == ag.n_atoms, 'no protein atoms?'
        oneletterseq = ag.residues.sequence().seq

        self.toppar = toppar
        self.ag     = ag

        self.generate()

    def generate(self):

        agnames = self.ag.names
        
        ### N-TERM
        if np.all(np.isin(['HT1', 'HT2', 'HT3'], agnames)):
            if self.ag.residues[0].resname == 'GLY':
                NPATCH = 'GLYP'
            else:
                NPATCH = 'NTER'

        if np.all(np.isin(['HY1', 'HY2', 'HY3'], agnames)):
            if self.ag.residues[0].resname == 'PRO':
                NPATCH = 'ACP'
            else:
                NPATCH = 'ACE'
        
        ### C-TERM
        if np.all(np.isin(['OT1', 'OT2'], agnames)):
            if np.isin(['HT2B'], agnames):
                CPATCH = 'CNEU'
            else:
                CPATCH = 'CTER'

        if np.all(np.isin(['CT', 'HT1', 'HT2', 'HT3'], agnames)):
            CPATCH = 'CT1'

        if np.all(np.isin(['NT', 'HT1', 'HT2'], agnames)):
            CPATCH = 'CT2'

        if np.all(np.isin(['CAT', 'HT1', 'HT2', 'HT3'], agnames)):
            CPATH = 'CT3'

        print('NPATCH: ', NPATCH)
        print('CPATCH: ', CPATCH, '\n')



        
        #names = np.zeros((0))
        #for residue in self.ag.residues:
        #    resn  = residue.resname
        #    names = np.concatenate([names, self.toppar.RESI[resn]['names']])

        #
        #agnames = self.ag.names
        #idx_firstN = np.where(agnames == 'N')[0][0]
        #idx_lastC  = np.where(agnames == 'C')[0][-1]
        #agnamescut = agnames[idx_firstN:idx_lastC+1] #C
        #
        #self.n1 = names
        #self.n2 = agnamescut
        #assert np.all(agnamescut == names[:-1]), 'missing atoms?'

        #
        #Nterm = agnames[:idx_firstN]
        #Cterm = agnames[idx_lastC:]
        

        
