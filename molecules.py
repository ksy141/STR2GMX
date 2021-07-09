import MDAnalysis as mda
import numpy      as np
import pandas     as pd
import os
from   .molecule  import Molecule

funcForBonds     = 1
funcForAngles    = 5    # Urey-Bradley angle type
funcForDihedrals = 9    # special type for treating multiple entries (modification in source code)
funcForImpropers = 2
funcFor14        = 1    # 1-4 interaction pair type
funcForCmap      = 1
funcForLonepairs = 3
funcForExclusions = 1
ptype = 'A'


class Molecules:
    
    def __init__(self, mols, toppar, prefix='./toppar/'):
        
        if not os.path.exists(prefix):
            os.makedirs(prefix)

        for mol in mols:
            assert isinstance(mol, Molecule), 'Molecule instance?'

        self.mols   = mols
        self.toppar = toppar
        self.prefix = prefix

        self.collect_used()
        self.write_forcefield()
        self.write_molecule()


    def collect_used(self):
        
        self.dbbonds     = []
        self.dbpairs     = []
        self.dbnbfix     = []
        self.dbangles    = []
        self.dbdihedrals = []
        self.dbdihewilds = []
        self.dbimpropers = []
        self.dbimprwilds = []
        self.atomtypes   = []
        

        ### CONNECTIVITY
        for mol in self.mols:
            ag = mol.ag
            self.atomtypes.extend(ag.types)
            
            # BONDS
            for i, j in mol.bsorted:
                t1, t2 = ag[i].type, ag[j].type
                if t1 > t2: t1, t2 = t2, t1
                app = [t1, t2, *self.toppar.BONDS[(t1, t2)]]
                if app not in self.dbbonds:
                    self.dbbonds.append(app)
            
            # ANGLES
            for i, j, k in mol.asorted:
                t1, t2, t3 = ag[i].type, ag[j].type, ag[k].type
                if t1 > t3: t1, t3 = t3, t1

                if (t1,t2,t3) not in self.toppar.ANGLES: continue
                app = [t1, t2, t3, *self.toppar.ANGLES[(t1,t2,t3)]]
                if app not in self.dbangles:
                    self.dbangles.append(app)


            # DIHEDRALS
            for i, j, k, l in mol.dsorted:
                t1, t2, t3, t4 = ag[i].type, ag[j].type, ag[k].type, ag[l].type
                
                if t1 > t4:
                    t1, t2, t3, t4 = t4, t3, t2, t1
                
                i4 = (t1, t2, t3, t4)
                if i4 in self.toppar.DIHEDRALS:
                    for m in self.toppar.DIHEDRALS[i4]:
                        app = [*i4, *m]
                        if app not in self.dbdihedrals:
                            self.dbdihedrals.append([*i4, *m])
                

                if t2 > t3: t2, t3 = t3, t2
                i4 = ('X', t2, t3, 'X')
                if i4 in self.toppar.DIHEDRALS:
                    for m in self.toppar.DIHEDRALS[i4]:
                        app = [*i4, *m]
                        if app not in self.dbdihewilds:
                            self.dbdihewilds.append([*i4,*m])

               
            
            # IMPROPERS
            for i, j, k, l in mol.isorted:
                t1, t2, t3, t4 = ag[i].type, ag[j].type, ag[k].type, ag[l].type

                if t1 > t4:
                    t1, t2, t3, t4 = t4, t3, t2, t1

                i4 = (t1, t2, t3, t4)
                if i4 in self.toppar.IMPROPER:
                    q0, cq = self.toppar.IMPROPER[i4]
                    app = [*i4, q0, cq]
                    if app not in self.dbimpropers:
                        self.dbimpropers.append(app)
                
                i4 = (t1, 'X', 'X', t4)
                if i4 in self.toppar.IMPROPER:
                    q0, cq = self.toppar.IMPROPER[i4]
                    app = [*i4, q0, cq]
                    if app not in self.dbimprwilds:
                        self.dbimprwilds.append([*i4, q0, cq])
            



        ### NONBONDED14 [pairtypes]
        self.atomtypes = set(self.atomtypes)
        fix = []; unfix = [];
        
        for attype in self.atomtypes:
            if attype in self.toppar.NB14:
                fix.append(attype)
            else:
                unfix.append(attype)
        
        
        for type1 in fix:
            isigma14, ieps14 = self.toppar.NB14[type1]
        
            for type2 in unfix:
                jsigma14, jeps14 = self.toppar.NONBONDED[type2]
        
                if ((type1, type2) in self.toppar.NBFIX) or ((type2, type1) in self.toppar.NBFIX):
                    continue
        
                sigma14 = (isigma14 + jsigma14)/2.0
                eps14   = (ieps14 * jeps14)**0.5
                app     = [type1, type2, sigma14, eps14]
                self.dbpairs.append(app)
        
        
        for i, type1 in enumerate(fix):
            isigma14, ieps14 = self.toppar.NB14[type1]
        
            for j, type2 in enumerate(fix):
                jsigma14, jeps14 = self.toppar.NB14[type2]
                
                if i > j: continue
                if ((type1, type2) in self.toppar.NBFIX) or ((type2, type1) in self.toppar.NBFIX):
                    continue
        
                sigma14 = (isigma14 + jsigma14)/2.0
                eps14   = (ieps14 * jeps14)**0.5
                app     = [type1, type2, sigma14, eps14]
                self.dbpairs.append(app)


        ### NBFIX
        for nbfix in self.toppar.NBFIX:
            type1, type2 = nbfix
            sigma, eps   = self.toppar.NBFIX[type1,type2]

            if type1 in self.atomtypes and type2 in self.atomtypes:
                if not [type1, type2, sigma, eps] in self.dbnbfix:
                    self.dbnbfix.append([type1,type2,sigma,eps])
    

        #### dihedraltypes
        #for dihedral in self.toppar.DIHEDRALS:
        #    type1, type2, type3, type4 = dihedral
        #    if type1 == type4 == 'X' and all([attype in psf.types for attype in [type2, type3]]):
        #        for m in toppar.params['dihedrals'][type1, type2, type3, type4]:
        #            phi, cp, mult = m
        #            dbdihewilds.append([type1,type2,type3,type4,phi,cp,mult])
        #    elif all([attype in psf.types for attype in [type1, type2, type3, type4]]):
        #        for m in toppar.params['dihedrals'][type1, type2, type3, type4]:
        #            phi, cp, mult = m
        #            dbdihedrals.append([type1,type2,type3,type4,phi,cp,mult])
    
        ## impropertypes
        #for improper in toppar.params['impropers']:
        #    type1, type2, type3, type4 = improper
        #    if type2 == type3 == 'X' and all([attype in psf.types for attype in [type1, type4]]):
        #        q0, cq = toppar.params['impropers'][type1,type2,type3,type4]
        #        dbimprwilds.append([type1,type2,type3,type4,q0,cq])
        #    elif all([attype in psf.types for attype in [type1, type2, type3, type4]]):
        #        q0, cq = toppar.params['impropers'][type1,type2,type3,type4]
        #        dbimpropers.append([type1,type2,type3,type4,q0,cq])





    def write_forcefield(self):
        itpFile = open(self.prefix + 'forcefield.itp', 'w')
        itpFile.write(';;\n')
        itpFile.write(';; Generated by str2gmx\n')
        itpFile.write(';; A previous CHARMM-GUI script, psf2itp.py, has been used\n')
        itpFile.write(';; Correspondance:\n')
        itpFile.write(';; siyoungkim@uchicago.edu\n')
        itpFile.write(';;\n')
        itpFile.write(';; CHARMM FF in GROMACS format\n')
        itpFile.write(';;\n\n')
     
        # defaults
        itpFile.write('\n[ defaults ]\n')
        itpFile.write('; nbfunc\tcomb-rule\tgen-pairs\tfudgeLJ\tfudgeQQ\n')
        itpFile.write('1\t2\tyes\t1.000000\t1.000000\n')


    
        # atomtypes
        fmt = ' %7s %5d %10.4f %10.3f %5s %20.11e %15.6e ; %19.11e %15.6e \n' 

        itpFile.write('\n[ atomtypes ]\n')
        itpFile.write('; name\tat.num\tmass\tcharge\tptype\tsigma\tepsilon\t;\tsigma_14\tepsilon_14\n')
        
        series = []
        for attype in sorted(self.atomtypes):
            mass   = self.toppar.ATOMS[attype]['mass']
            elem   = self.toppar.ATOMS[attype]['elem']
            numb   = self.toppar.ATOMS[attype]['numb']
            charge = 0
            sigma, eps = self.toppar.NONBONDED[attype]

            if attype in self.toppar.NB14:
                sigma14, eps14 = self.toppar.NB14[attype]
            else:
                sigma14, eps14 = self.toppar.NONBONDED[attype]

            series.append([attype, numb, mass, charge, ptype, sigma, eps, sigma14, eps14])

        df = pd.DataFrame(series)
        ff = df.sort_values(by=[1,0], ignore_index=True)
        for i in range(len(ff)):
            itpFile.write(fmt %tuple(ff.loc[i]))


        # nbfix
        fmt = '%7s %7s %5d %18.11e %18.11e \n'
        if len(self.dbnbfix) > 0:
            itpFile.write('\n[ nonbond_params ]\n')
            itpFile.write('; i\tj\tfunc\tsigma\tepsilon\n')

            df = pd.DataFrame(self.dbnbfix)
            ff = df.sort_values(by=[0,1], ignore_index=True)

            for i in range(len(ff)):
                itpFile.write(fmt %tuple(ff.loc[i]))
    

        # bondtypes
        fmt = '%7s %7s %5d %13.6e %13.6e\n'
        if len(self.dbbonds) > 0:
            df = pd.DataFrame(self.dbbonds)
            ff = df.sort_values(by=[0,1]).to_numpy()

            itpFile.write('\n[ bondtypes ]\n')
            itpFile.write('; i\tj\tfunc\tb0\tKb\n')
            for bond in ff:
                type1, type2, b0, Kb = bond
                itpFile.write(fmt %(type1, type2, funcForBonds, b0, Kb))
    
        # pairtypes
        fmt = '%7s %7s %5d %18.11e %18.11e \n'
        if len(self.dbpairs) > 0:
            itpFile.write('\n[ pairtypes ]\n')
            itpFile.write('; i\tj\tfunc\tsigma1-4\tepsilon1-4\n')
            for pair in self.dbpairs:
                type1, type2, sigma14, eps14 = pair
                itpFile.write(fmt %(type1, type2, funcFor14, sigma14, eps14))
    
        # angletypes
        fmt = '%7s %7s %7s %5d %14.7e %14.7e %14.7e %14.7e\n'
        if len(self.dbangles) > 0:
            itpFile.write('\n[ angletypes ]\n')
            itpFile.write('; i\tj\tk\tfunc\tth0\tKth\ts0\tKub\n')
            for angle in self.dbangles:
                type1, type2, type3, th0, cth, S0, Kub = angle
                itpFile.write(fmt %(type1, type2, type3, funcForAngles, th0, cth, S0, Kub))
    
        # dihedraltypes
        fmt = '%7s %7s %7s %7s %5d %13.6e %13.6e %6d\n'
        if len(self.dbdihedrals) > 0 or len(self.dbdihewilds) > 0:
            itpFile.write('\n[ dihedraltypes ]\n')
            itpFile.write('; i\tj\tk\tl\tfunc\tphi0\tKphi\tmult\n')
            for dihedral in self.dbdihedrals:
                type1, type2, type3, type4, phi, cp, mult = dihedral
                itpFile.write(fmt %(type1, type2, type3, type4, funcForDihedrals, phi, cp, mult))
            for dihedral in self.dbdihewilds:
                type1, type2, type3, type4, phi, cp, mult = dihedral
                itpFile.write(fmt %(type1, type2, type3, type4, funcForDihedrals, phi, cp, mult))
    
        # impropertypes
        fmt = '%7s %7s %7s %7s %5d %13.6e %13.6e\n'
        if len(self.dbimpropers) > 0 or len(self.dbimprwilds) > 0:
            itpFile.write('\n[ dihedraltypes ]\n')
            itpFile.write('; i\tj\tk\tl\tfunc\tq0\tKq\n')
            for improper in self.dbimpropers:
                type1, type2, type3, type4, q0, cq = improper
                itpFile.write(fmt %(type1, type2, type3, type4, funcForImpropers, q0, cq))
            for improper in self.dbimprwilds:
                type1, type2, type3, type4, q0, cq = improper
                itpFile.write(fmt %(type1, type2, type3, type4, funcForImpropers, q0, cq))

        itpFile.close()

                
    def write_molecule(self):
        for mol in self.mols:
            nrexcl = 1
            if len(mol.bsorted) > 1:
                nrexcl = 2
            if len(mol.dsorted) > 0:
                nrexcl = 3
            
            itpFile = open(self.prefix + mol.resname.upper() + '.itp', 'w')
            itpFile.write(';;\n')
            itpFile.write(';; Generated by str2gmx\n')
            itpFile.write(';; A previous CHARMM-GUI script, psf2itp.py, has been used\n')
            itpFile.write(';; Correspondance:\n')
            itpFile.write(';; siyoungkim@uchicago.edu\n')
            itpFile.write(';;\n')
            itpFile.write(';; GROMACS topology file for %s\n' % mol.resname.upper())
            itpFile.write(';;\n\n')
            
            itpFile.write('\n[ moleculetype ]\n')
            itpFile.write('; name\tnrexcl\n')
            itpFile.write('%s\t %5d\n' % (mol.resname.upper(), nrexcl))
            
    
            # ATOMS
            itpFile.write('\n[ atoms ]\n')
            itpFile.write('; nr	type	resnr	residu	atom	cgnr	charge	mass\n')
    
            qtot = 0
            fmt = ' %5d %10s %6s %8s %6s %6d %12.6f %10.4f   ; qtot %6.3f\n'
            
            for i, atom in enumerate(mol.ag):
                qtot += atom.charge
                itpFile.write(fmt %(i+1, atom.type, 1, atom.resname, atom.name, i+1, atom.charge, atom.mass, qtot))
    
    
            # BONDS
            if len(mol.bsorted) > 0:
                itpFile.write('\n[ bonds ]\n')
                itpFile.write('; ai\taj\tfunct\tb0\tKb\n')
                for i, j in mol.bsorted:
                    itpFile.write('%5d %5d %5d\n' % (i + 1, j + 1, funcForBonds))
    
    
            # PAIRS
            if len(mol.psorted) > 0:
                itpFile.write('\n[ pairs ]\n')
                itpFile.write('; ai\taj\tfunct\tc6\tc12\n')
                for i, j in mol.psorted:
                    itpFile.write('%5d %5d %5d\n' % (i + 1, j + 1, funcFor14))
    
    
            # ANGLES
            if len(mol.asorted) > 0:
                itpFile.write('\n[ angles ] \n')
                itpFile.write('; ai\taj\tak\tfunct\tth0\tcth\tS0\tKub\n')
                for angle in mol.asorted:
                    i, j, k = angle
                    itpFile.write('%5d %5d %5d %5d\n' % (i+1, j+1, k+1, funcForAngles))
    
    
            # DIHEDRALS
            if len(mol.dsorted) > 0:
                itpFile.write('\n[ dihedrals ]\n')
                itpFile.write('; ai\taj\tak\tal\tfunct\tphi0\tcp\tmult\n')
                for dihedral in mol.dsorted:
                    i, j, k, l = dihedral
                    itpFile.write('%5d %5d %5d %5d %5d\n' % (i+1, j+1, k+1, l+1, funcForDihedrals))
    
    
            # IMPROPER
            if len(mol.isorted) > 0:
                itpFile.write('\n[ dihedrals ]\n')
                itpFile.write('; ai\taj\tak\tal\tfunct\tq0\tcq\n')
                for improper in mol.isorted:
                    i, j, k, l = improper
                    itpFile.write('%5d %5d %5d %5d %5d\n' % (i+1, j+1, k+1, l+1, funcForImpropers))
            
            itpFile.close()

