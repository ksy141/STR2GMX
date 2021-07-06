import numpy  as np
import pandas as pd
import glob


comment = '!'
kcal2kJ = 4.184
items   = ('RESI', 'ATOMS', 'BONDS', 'ANGLES', 'DIHEDRALS', 'IMPROPER', 'NONBONDED', 'NBFIX', 'CMAP') 


class ReadToppars:

    def __init__(self, toppars=None):
        
        # a list of files to read
        if isinstance(toppars, list):
            pass
        
        # indicate a path
        elif isinstance(toppars, str):
            toppars = glob.glob(toppars + '/*rtf') + \
                      glob.glob(toppars + '/*prm') + \
                      glob.glob(toppars + '/*str') 

        if not toppars:
            raise ValueError('please provide a list of top/prm/str or a path to them')
        
        print(pd.DataFrame([toppar.split('/')[-1] for toppar in toppars]))
        self.toppars = toppars
        self.parse()
        self.collect()



    def parse(self):
        '''Remove comments and parse useful lines'''

        save  = ''
        saves = []

        read_type = False
        for toppar in self.toppars:
            for line in open(toppar, 'r'):

                rline = line.rstrip()
                if rline: # Pass empty line

                    sline = rline.split(comment)[0].rstrip()
                    if sline: # Pass a line starting with !

                        if sline.startswith('HBOND'):
                            continue

                        if sline.startswith(items):
                            read_type = sline.split()[0]
                            saves.append(save)
                            save = ''

                        if sline.startswith(('END', 'end')):
                            read_type = False
                            saves.append(save)
                            save = ''
                        
                        # READ EACH COMPONENT
                        if read_type:
                            if sline.endswith('-'):
                                save += sline
                            else:
                                save += sline + '\n'

        self.saves = saves



    def collect(self):
        self.RESI      = {}
        self.ATOMS     = {}
        self.BONDS     = {}
        self.ANGLES    = {}
        self.DIHEDRALS = {}
        self.IMPROPER  = {}
        self.NONBONDED = {}
        self.NB14      = {}
        self.NBFIX     = {}


        for save in self.saves:
            if not save: continue
            ssave = save.rstrip().split('\n')

            if save.startswith('RESI'):
                self._read_RESI(ssave)

            elif save.startswith('ATOMS'):
                self._read_ATOMS(ssave)

            elif save.startswith('BONDS'):
                self._read_BONDS(ssave)

            elif save.startswith('ANGLES'):
                self._read_ANGLES(ssave)

            elif save.startswith('DIHEDRALS'):
                self._read_DIHEDRALS(ssave)

            elif save.startswith('IMPROPER'):
                self._read_IMPROPER(ssave)
            
            elif save.startswith('NONBONDED'):
                self._read_NONBONDED(ssave)

            elif save.startswith('NBFIX'):
                self._read_NBFIX(ssave)



    def _read_RESI(self, ssave):
        resname = ssave[0].split()[1]

        names   = []
        types   = []
        charges = []
        bonds   = []
        imprs   = []

        for line in ssave:
            segments = line.split()

            if line.startswith('ATOM'):
                names.append(segments[1])
                types.append(segments[2])
                charges.append(float(segments[3]))

            elif line.startswith(('BOND', 'DOUBLE')):
                for i in range(0, len(segments[1:]), 2):
                    bonds.append([segments[i+1], segments[i+2]])

            elif line.startswith('IMPR'):
                for i in range(0, len(segments[1:]), 4):
                    imprs.append([segments[i+1], segments[i+2], segments[i+3], segments[i+4]])


        self.RESI[resname] = {'names':   names, 
                              'types':   types, 
                              'bonds':   bonds, 
                              'imprs':   imprs,
                              'charges': charges}



    def _read_ATOMS(self, ssave):
        for line in ssave:
            if not line.startswith('MASS'): continue

            segments = line.split()
            if segments[2] not in self.ATOMS.keys():
                self.ATOMS[segments[2]] = float(segments[3])


    def _read_BONDS(self, ssave):
        for line in ssave[1:]:
            segments = line.split()
            type1 = segments[0]
            type2 = segments[1]

            if type1 > type2:
                type1, type2 = type2, type1

            # kcal/mol/A**2 -> kJ/mol/nm**2 incl factor 2
            Kb    = float(segments[2]) * 2 * kcal2kJ * 1000 / 10
            b0    = float(segments[3]) / 10

            self.BONDS[(type1, type2)] = [b0, Kb]

        
        
    def _read_ANGLES(self, ssave):
        for line in ssave[1:]:
            segments = line.split()

            type1 = segments[0]
            type2 = segments[1]
            type3 = segments[2]

            if type1 > type3:
                type1, type3 = type3, type1

            cth   = float(segments[3]) * 2 * kcal2kJ
            th0   = float(segments[4])

            Kub   = 0.0
            S0    = 0.0
            
            if len(segments) > 6: # check for Urey-Bradley parameters
                Kub = float(segments[5]) * 2 * kcal2kJ * 1000 / 10
                S0  = float(segments[6]) / 10


            self.ANGLES[(type1, type2, type3)] = [th0, cth, S0, Kub]
 



    def _read_DIHEDRALS(self, ssave):
        for line in ssave[1:]:
            segments = line.split()

            type1 = segments[0]
            type2 = segments[1]
            type3 = segments[2]
            type4 = segments[3]
            
            cp   = float(segments[4]) * kcal2kJ
            mult = int(segments[5])
            phi0 = float(segments[6])
            
            if (type1, type2, type3, type4) in self.DIHEDRALS:
                for imult, dihedral in enumerate(self.DIHEDRALS[type1, type2, type3, type4]):
                    if dihedral[2] == mult: del self.DIHEDRALS[type1, type2, type3, type4][imult]
                self.DIHEDRALS[type1, type2, type3, type4].append([phi0, cp, mult])


            elif (type4, type3, type2, type1) in self.DIHEDRALS:
                for imult, dihedral in enumerate(self.DIHEDRALS[type4, type3, type2, type1]):
                    if dihedral[2] == mult: del self.DIHEDRALS[type4, type3, type2, type1][imult]
                self.DIHEDRALS[type4, type3, type2, type1].append([phi0, cp, mult])

            else:
                self.DIHEDRALS[type1, type2, type3, type4] = [[phi0, cp, mult]]
 


    def _read_IMPROPER(self, ssave):
        for line in ssave[1:]:
            segments = line.split()

            type1 = segments[0]
            type2 = segments[1]
            type3 = segments[2]
            type4 = segments[3]

            cq = float(segments[4]) * 2 * kcal2kJ
            q0 = float(segments[6])

            self.IMPROPER[(type1, type2, type3, type4)] = [q0, cq]



    def _read_NONBONDED(self, ssave):
        ### LENNARD JONES
        ### [CHARMM]   V_LJ = 1 * eps * ( (r_0/ r)^12  - 2 * (r_0 / r)^6)
        ### [GROMACS]  V_LJ = 4 * eps * ( (sigma/r)^12 - (sigma/r)^6 )
        ### r_0 = RminHalf * 2 = 2^(1/6) * sigma

        for line in ssave[1:]:
            segments = line.split()
        
            attype   = segments[0]
            epsilon  = float(segments[2])
            RminHalf = float(segments[3])
            eps      = abs(epsilon*kcal2kJ)                # conversion to kJ and positive
            sigma    = 2*RminHalf/(10.0*2.0**(1.0/6.0))    # -> nm, double distance and rmin2sigma factor

            self.NONBONDED[attype] = [sigma, eps]

            if len(segments)> 6:        # test length to avoid IndexError
                try:                    # if possible, convert element 5 to float
                    segments[5] = float(segments[5])
                except:
                    None

                # is segment 5 and 6 floats => there's 1-4 defined
                if not isinstance(segments[5], str):                  # not string?
                    epsilon14  = float(segments[5])                   # read charmm epsilon
                    eps14      = abs(epsilon14*kcal2kJ)               # conversion to gromacs units
                    Rmin14Half = float(segments[6])                   # read charmm Rmin*1/2
                    sigma14    = 2*Rmin14Half/(10.0*2.0**(1.0/6.0))   # conversion to gromacs units
                    self.NB14[attype] = [sigma14, eps14]              # add to list
    


    def _read_NBFIX(self, ssave):
        for line in ssave[1:]:
            if not line: continue
            if line.startswith('HBOND'): continue

            segments = line.split()
            type1 = segments[0]
            type2 = segments[1]

            if type1 > type2:
                type1, type2 = type2, type1

            epsilon = float(segments[2])
            Rmin    = float(segments[3])
            eps     = abs(epsilon*kcal2kJ)          # conversion to kJ and positive
            sigma   = Rmin/(10.0*2.0**(1.0/6.0))    # -> nm, double distance and rmin2sigma factor
            
            self.NBFIX[(type1, type2)] = [sigma, eps]




