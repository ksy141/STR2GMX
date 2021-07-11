import MDAnalysis as mda
import numpy as np

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
        with open('FF/C36/toppar/tip216.crd') as f:
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


    def run(self, ag = None, z = None, point = None):
        if point == None:
            raise ValueError('define a point')

        if ag == None and z == None:
            raise ValueError('either define AtomGroup or z value')


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
        bAz = water_pos[:,2][::3] > self.pbc[2]
        bA  = bAx | bAy | bAz
        waterbox_pos = water_pos[np.repeat(~bA, 3)]








