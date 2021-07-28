import MDAnalysis as mda
import numpy as np
import pandas as pd

class Chain:

    def __init__(self, ag, toppar, segname = 'PROA', add_bonds = []):

        assert ag.select_atoms('protein').n_atoms == ag.n_atoms, 'no protein atoms?'
        oneletterseq = ag.residues.sequence().seq
        
        self.resname = segname
        self.toppar = toppar
        self.system = ag.universe
        self.ag     = ag
        self.u      = mda.Merge(ag)
        self.ag     = self.u.atoms
        self.add_bonds = add_bonds
        self.u.add_TopologyAttr('charges', [0] * self.u.atoms.n_atoms)
        self.generate()

    def generate(self):
        ### CURRENTLY SUPPORT
        ### NTER/GLYP, ACP/ACE
        ### CTER, CNEU, CT1, CT2, CT3

        atnames  = self.u.atoms.names
        residues = self.u.residues
        
        ### N-TERM
        if np.all(np.isin(['HT1', 'HT2', 'HT3'], atnames)):
            if self.u.residues[0].resname == 'GLY':
                NPATCH = 'GLYP'
            else:
                NPATCH = 'NTER'


        if ('HT1' in atnames) and ('HT2' in atnames) and ('HT3' not in atnames):
            if self.u.residues[0].resname == 'GLY':
                NPATCH = 'NGNE'
            else:
                NPATCH = 'NNEU'

        if np.all(np.isin(['HY1', 'HY2', 'HY3'], atnames)):
            if self.u.residues[0].resname == 'PRO':
                NPATCH = 'ACP'
            else:
                NPATCH = 'ACE'
        
        ### C-TERM
        if np.all(np.isin(['OT1', 'OT2'], atnames)):
            if np.isin(['HT2B'], atnames):
                CPATCH = 'CNEU'
            else:
                CPATCH = 'CTER'

        if np.all(np.isin(['CT', 'HT1', 'HT2', 'HT3'], atnames)):
            CPATCH = 'CT1'

        if np.all(np.isin(['NT', 'HT1', 'HT2'], atnames)):
            CPATCH = 'CT2'

        if np.all(np.isin(['CAT', 'HT1', 'HT2', 'HT3', 'HNT'], atnames)):
            CPATCH = 'CT3'

        print('NPATCH: ', NPATCH)
        print('CPATCH: ', CPATCH, '\n')
        

        ### BOND / ANGLES / IMPRS
        bb = []; bnames = []; btypes = [];
        cc = []; ctypes = [];
        self.isorted = []

        
        for i in range(residues.n_residues):
            residue  = residues[i]
            resn     = residue.resname
            atoms    = residue.atoms
            names    = atoms.names

            if resn == 'CYS' and 'HG1' not in names:
                resn = 'CYSD'

            top_resi = self.toppar.RESI[resn]
            top_nter = self.toppar.RESI[NPATCH]
            top_cter = self.toppar.RESI[CPATCH]
            
            ### IF ONLY ONE RESIDUE
            if residues.n_residues == 1:
                for atom in atoms:
                    if atom.name in   top_nter['names']:
                        idx         = top_nter['names'] == atom.name
                        atom.type   = top_nter['types'][idx][0]
                        atom.charge = top_nter['charges'][idx][0]
                        atom.mass   = top_nter['masses'][idx][0]

                    elif atom.name in top_cter['names']:
                        idx         = top_cter['names'] == atom.name
                        atom.type   = top_cter['types'][idx][0]
                        atom.charge = top_cter['charges'][idx][0]
                        atom.mass   = top_cter['masses'][idx][0]

                    elif atom.name in top_resi['names']:
                        idx         = top_resi['names'] == atom.name
                        atom.type   = top_resi['types'][idx][0]
                        atom.charge = top_resi['charges'][idx][0]
                        atom.mass   = top_resi['masses'][idx][0]

                    else:
                        print('Missing atoms', atom)

                    assert np.sum(idx) == 1, 'Duplicate atoms?'

                b_tmp = np.concatenate([top_nter['bonds'], top_resi['bonds'], top_cter['bonds']])
                for ai, aj in b_tmp:
                    atomi = atoms[ai == names]
                    atomj = atoms[aj == names]
                    if atomi.n_atoms * atomj.n_atoms == 0: continue
                    assert atomi.n_atoms * atomj.n_atoms == 1, 'Duplicate atoms?'
                    bb.append([atomi[0].index, atomj[0].index])
                    bnames.append([atomi[0].name, atomj[0].name])
                    btypes.append([atomi[0].type, atomj[0].type])
 
                for ai, aj, ak, al in top_nter['imprs']:
                    atomi = atoms[ai == names]
                    atomj = atoms[aj == names]
                    atomk = atoms[ak == names]
                    atoml = atoms[al == names]
                    assert atomi.n_atoms * atomj.n_atoms * atomk.n_atoms * atoml.n_atoms == 1, 'Check N-PATCH improper?'
                    self.isorted.append([atomi[0].index, atomj[0].index, atomk[0].index, atoml[0].index])

                for ai, aj, ak, al in top_cter['imprs']:
                    atomi = atoms[ai == names]
                    atomj = atoms[aj == names]
                    atomk = atoms[ak == names]
                    atoml = atoms[al == names]
                    assert atomi.n_atoms * atomj.n_atoms * atomk.n_atoms * atoml.n_atoms == 1, 'Check C-PATCH improper?'
                    self.isorted.append([atomi[0].index, atomj[0].index, atomk[0].index, atoml[0].index])

                for ai, aj, ak, al in top_resi['imprs']:
                    atomi = atoms[ai == names]
                    atomj = atoms[aj == names]
                    atomk = atoms[ak == names]
                    atoml = atoms[al == names]
                    if atomi.n_atoms * atomj.n_atoms * atomk.n_atoms * atoml.n_atoms == 1:
                        self.isorted.append([atomi[0].index, atomj[0].index, atomk[0].index, atoml[0].index])

 
                self.bnames = bnames
                self.btypes = btypes
                
                df = pd.DataFrame(np.unique(bb, axis=0))
                self.bsorted = df.sort_values(by=[0,1]).to_numpy()
                self.u.add_TopologyAttr('bonds', self.bsorted)
                
                aa = mda.topology.guessers.guess_angles(self.u.bonds)
                df = pd.DataFrame(aa)
                self.asorted = df.sort_values(by=[1,0,2]).to_numpy()
                self.u.add_TopologyAttr('angles', self.asorted)

                dd = mda.topology.guessers.guess_dihedrals(self.u.angles)
                df = pd.DataFrame(dd)
                self.dsorted = df.sort_values(by=[1,2,0,3]).to_numpy()
                self.u.add_TopologyAttr('dihedrals', self.dsorted)

                self.psorted = self._guess_pairs14(self.bsorted)
                
                self.csorted = cc
                self.ctypes  = ctypes

                self.isorted = np.unique(self.isorted, axis=0)
                return

 
           
            ### FIRST RESIDUE
            if i == 0:
                top_nter = self.toppar.RESI[NPATCH]

                nextN    = residues[i+1].atoms.select_atoms('name N')[0]
                currC    = atoms.select_atoms('name C')[0]

                for atom in atoms:
                    if atom.name in   top_nter['names']:
                        idx         = top_nter['names'] == atom.name
                        atom.type   = top_nter['types'][idx][0]
                        atom.charge = top_nter['charges'][idx][0]
                        atom.mass   = top_nter['masses'][idx][0]

                    elif atom.name in top_resi['names']:
                        idx         = top_resi['names'] == atom.name
                        atom.type   = top_resi['types'][idx][0]
                        atom.charge = top_resi['charges'][idx][0]
                        atom.mass   = top_resi['masses'][idx][0]

                    else:
                        print('Missing atoms in N-PATCH?')

                    assert np.sum(idx) == 1, 'Duplicate atoms?'
                
                b_tmp = np.concatenate([top_nter['bonds'], top_resi['bonds']])
                for ai, aj in b_tmp:
                    atomi = atoms[ai == names]
                    atomj = atoms[aj == names]
                    if atomi.n_atoms * atomj.n_atoms == 0: continue
                    assert atomi.n_atoms * atomj.n_atoms == 1, 'Duplicate atoms?'
                    bb.append([atomi[0].index, atomj[0].index])
                    bnames.append([atomi[0].name, atomj[0].name])
                    btypes.append([atomi[0].type, atomj[0].type])
                
                bb.append([currC.index, nextN.index])
                bnames.append([currC.name, nextN.name])
                btypes.append([currC.type, nextN.type])
                

                for ai, aj, ak, al in top_nter['imprs']:
                    atomi = atoms[ai == names]
                    atomj = atoms[aj == names]
                    atomk = atoms[ak == names]
                    atoml = atoms[al == names]
                    assert atomi.n_atoms * atomj.n_atoms * atomk.n_atoms * atoml.n_atoms == 1, 'Check N-PATCH improper?'
                    self.isorted.append([atomi[0].index, atomj[0].index, atomk[0].index, atoml[0].index])

                impr1 = np.concatenate([atoms.select_atoms('name C').indices, \
                                        atoms.select_atoms('name CA').indices, \
                                        [nextN.index], \
                                        atoms.select_atoms('name O').indices])

                #impr2 = np.concatenate([atoms.select_atoms('name N').indices, \
                #                        atoms.select_atoms('name CAY').indices, \
                #                        atoms.select_atoms('name CA').indices, \
                #                        atoms.select_atoms('name HN CD').indices])

                assert len(impr1) == 4, 'check impr1'
                self.isorted.append(impr1)
                #if len(impr2) == 4: self.isorted.append(impr2)

            

                for ai, aj, ak, al in top_resi['imprs']:
                    atomi = atoms[ai == names]
                    atomj = atoms[aj == names]
                    atomk = atoms[ak == names]
                    atoml = atoms[al == names]
                    if atomi.n_atoms * atomj.n_atoms * atomk.n_atoms * atoml.n_atoms == 1:
                        self.isorted.append([atomi[0].index, atomj[0].index, atomk[0].index, atoml[0].index])

             
            ### LAST RESIDUE
            elif i == residues.n_residues - 1:
                top_cter = self.toppar.RESI[CPATCH]
                prevC    = residues[i-1].atoms.select_atoms('name C')[0]
                currN    = atoms.select_atoms('name N')[0]
                
                for atom in atoms:
                    if atom.name in   top_cter['names']:
                        idx         = top_cter['names'] == atom.name
                        atom.type   = top_cter['types'][idx][0]
                        atom.charge = top_cter['charges'][idx][0]
                        atom.mass   = top_cter['masses'][idx][0]

                    elif atom.name in top_resi['names']:
                        idx         = top_resi['names'] == atom.name
                        atom.type   = top_resi['types'][idx][0]
                        atom.charge = top_resi['charges'][idx][0]
                        atom.mass   = top_resi['masses'][idx][0]

                    else:
                        print('Missing atoms in C-PATCH?')

                    assert np.sum(idx) == 1, 'Duplicate atoms?'
                
                b_tmp = np.concatenate([top_cter['bonds'], top_resi['bonds']])
                for ai, aj in b_tmp:
                    atomi = atoms[ai == names]
                    atomj = atoms[aj == names]
                    if atomi.n_atoms * atomj.n_atoms == 0: continue
                    assert atomi.n_atoms * atomj.n_atoms == 1, 'Duplicate atoms?'
                    bb.append([atomi[0].index, atomj[0].index])
                    bnames.append([atomi[0].name, atomj[0].name])
                    btypes.append([atomi[0].type, atomj[0].type])
                
                bb.append([prevC.index, currN.index])
                bnames.append([prevC.name, currN.name])
                btypes.append([prevC.type, currN.type])
                
                for ai, aj, ak, al in top_cter['imprs']:
                    atomi = atoms[ai == names]
                    atomj = atoms[aj == names]
                    atomk = atoms[ak == names]
                    atoml = atoms[al == names]
                    assert atomi.n_atoms * atomj.n_atoms * atomk.n_atoms * atoml.n_atoms == 1, 'Check C-PATCH improper?'
                    self.isorted.append([atomi[0].index, atomj[0].index, atomk[0].index, atoml[0].index])

                impr1 = np.concatenate([atoms.select_atoms('name N').indices, \
                                        [prevC.index], \
                                        atoms.select_atoms('name CA').indices, \
                                        atoms.select_atoms('name HN CD').indices])

                #impr2 = np.concatenate([atoms.select_atoms('name C').indices, \
                #                        atoms.select_atoms('name CA').indices, \
                #                        atoms.select_atoms('name NT').indices, \
                #                        atoms.select_atoms('name O').indices])

                assert len(impr1) == 4, 'check impr1'
                self.isorted.append(impr1)
                #if len(impr2) == 4: self.isorted.append(impr2)

                for ai, aj, ak, al in top_resi['imprs']:
                    atomi = atoms[ai == names]
                    atomj = atoms[aj == names]
                    atomk = atoms[ak == names]
                    atoml = atoms[al == names]
                    if atomi.n_atoms * atomj.n_atoms * atomk.n_atoms * atoml.n_atoms == 1:
                        self.isorted.append([atomi[0].index, atomj[0].index, atomk[0].index, atoml[0].index])

 

            ### THE REST
            else:
                assert np.all(names == top_resi['names']), 'missing atoms?'
                residue.atoms.masses  = top_resi['masses']
                residue.atoms.types   = top_resi['types']
                residue.atoms.charges = top_resi['charges']
                
                ### BONDS
                for ai, aj in top_resi['bonds']:
                    atomp = atoms[ai == names][0].index
                    if aj == '+N':
                        atomq = residues[i+1].atoms.select_atoms('name N')[0].index
                    else:
                        atomq = atoms[aj == names][0].index
                    bb.append([atomp, atomq])
                

                ### IMPRS
                for ai, aj, ak, al in top_resi['imprs']:
                    atomi = atoms[ai == names][0].index
                    atoml = atoms[al == names][0].index

                    if aj == '-C':
                        atomj = residues[i-1].atoms.select_atoms('name C')[0].index
                    else:
                        atomj = atoms[aj == names][0].index

                    if ak == '+N':
                        atomk = residues[i+1].atoms.select_atoms('name N')[0].index
                    else:
                        atomk = atoms[ak == names][0].index

                    self.isorted.append([atomi, atomj, atomk, atoml])
                
                for ai, aj, ak, al in top_resi['imprs']:
                    atomi = atoms[ai == names]
                    atomj = atoms[aj == names]
                    atomk = atoms[ak == names]
                    atoml = atoms[al == names]
                    if atomi.n_atoms * atomj.n_atoms * atomk.n_atoms * atoml.n_atoms == 1:
                        self.isorted.append([atomi[0].index, atomj[0].index, atomk[0].index, atoml[0].index])

 
        
        ### CMAPS
        if residues.n_residues == 1:
            pass

        elif residues.n_residues == 2:
            cm = ['CY', 'N', 'CA', 'C', '+N']
            tt, nn, ii = self.name2atom(cm, [residues[0], residues[0], residues[1]])
            if len(tt) > 0:
                cc.append(ii)
                ctypes.append(tt)
            
            if CPATCH != 'CT2':
                # CT2 has NT but CT2 doesn't have CMAP
                cm = ['-C', 'N', 'CA', 'C', 'NT']
                tt, nn, ii = self.name2atom(cm, [residues[0], residues[1], residues[1]])
                if len(tt) > 0:
                    cc.append(ii)
                    ctypes.append(tt)
 

        else:
            cm = ['CY', 'N', 'CA', 'C', '+N']
            tt, nn, ii = self.name2atom(cm, [residues[0], residues[0], residues[1]])
            if len(tt) > 0:
                cc.append(ii)
                ctypes.append(tt)
            
            if CPATCH != 'CT2':
                # CT2 has NT but CT2 doesn't have CMAP
                cm = ['-C', 'N', 'CA', 'C', 'NT']
                tt, nn, ii = self.name2atom(cm, [residues[-2], residues[-1], residues[-1]])
                if len(tt) > 0:
                    cc.append(ii)
                    ctypes.append(tt)
       
            for i in range(residues.n_residues):
                if (0 <= i - 1) and (i + 1 < residues.n_residues):
                    for cm in top_resi['cmaps']:
                        tt, nn, ii = self.name2atom(cm, [residues[i-1], residues[i], residues[i+1]])
                        if len(tt) > 0:
                            cc.append(ii)
                            ctypes.append(tt)


                
        self.bnames = bnames
        self.btypes = btypes
        
        if len(self.add_bonds) > 0:
            print('Adding additional bonds (DISU)')
            try:
                ### [[ag1, ag2], [ag3, ag4]] or [[atom1, atom2], [atom3, atom4]]
                for atomi, atomj in self.add_bonds:
                    try:
                        assert atomi.n_atoms * atomj.n_atoms == 1, 'more than two atoms?'
                        bb.append([atomi.indices[0], atomj.indices[0]])
                        self.bnames.append([atomi.names[0], atomj.names[0]])
                        self.btypes.append([atomi.types[0], atomj.types[0]])
                    except:
                        bb.append([atomi.index, atomj.index])
                        self.bnames.append([atomi.name, atomj.name])
                        self.btypes.append([atomi.type, atomj.type])
 
            except:
                ### [ag1, ag2] or [atom1, atom2]
                atomi, atomj = self.add_bonds
                try:
                    assert atomi.n_atoms * atomj.n_atoms == 1, 'more than two atoms?'
                    bb.append([atomi.indices[0], atomj.indices[0]])
                    self.bnames.append([atomi.names[0], atomj.names[0]])
                    self.btypes.append([atomi.types[0], atomj.types[0]])
                except:
                    bb.append([atomi.index, atomj.index])
                    self.bnames.append([atomi.name, atomj.name])
                    self.btypes.append([atomi.type, atomj.type])
 
        df = pd.DataFrame(np.unique(bb, axis=0))
        self.bsorted = df.sort_values(by=[0,1]).to_numpy()
        self.u.add_TopologyAttr('bonds', self.bsorted)
        
        aa = mda.topology.guessers.guess_angles(self.u.bonds)
        df = pd.DataFrame(aa)
        self.asorted = df.sort_values(by=[1,0,2]).to_numpy()
        self.u.add_TopologyAttr('angles', self.asorted)

        dd = mda.topology.guessers.guess_dihedrals(self.u.angles)
        df = pd.DataFrame(dd)
        self.dsorted = df.sort_values(by=[1,2,0,3]).to_numpy()
        self.u.add_TopologyAttr('dihedrals', self.dsorted)

        self.psorted = self._guess_pairs14(self.bsorted)
        
        self.csorted = cc
        self.ctypes  = ctypes

        self.isorted = np.unique(self.isorted, axis=0)

        
    def name2atom(self, names = [], resgroup = []):
        # names = ['-C', 'N', 'CA', 'C', '+N']
        # resgroup = [previous residue, current residue, next residue]
        ag0 = resgroup[0].atoms
        ag1 = resgroup[1].atoms
        ag2 = resgroup[2].atoms
        
        atoms = []
        types = []
        indes = []
        
        for name in names:
            ## +N [C-term]
            if name.startswith('+'):
                name = name[1:]
                atom = ag2.select_atoms('name %s' %name)

            ## -C [N-term]
            elif name.startswith('-'):
                name = name[1:]
                atom = ag0.select_atoms('name %s' %name)

            else:
                atom = ag1.select_atoms('name %s' %name)
            
            if atom.n_atoms == 0:
                return types, [], indes
            
            atoms += [atom[0]]

       
        for atom in atoms:
            types.append(atom.type)
            indes.append(atom.index)

        return types, [], indes
            


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

    