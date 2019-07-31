# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 15:15:39 2019

@author: Ali
"""

import numpy as np
from scipy.special import erfc
from numpy.linalg import norm, inv, det
from numpy import dot
        
def birch(V, V0, K0, K1):
    """Birch equation of state, fit to DFT total energy vs V and extract equillibrium V"""
    term1 = 1 / (K1 * (K1 - 1)) * (V / V0 ) ** (1 - K1)
    term2 = V / (K1 * V0) - 1 / (K1 - 1)
    return E0 + K0 * V0 * (term1 + term2)

def uniformK(a1, a2, a3, n):
    """Creates a uniform n by n /by n k-space grid give Direct lattice vectors as input
 
        a1, a2, a3: the Direct lattice vectors
        n: grid spacing
        
        
    """
    a1 = np.array(a1)
    a2 = np.array(a2)
    a3 = np.array(a3)
    
    #unit cell volume
    V = np.dot(a1 ,np.cross(a2, a3))
    
    #reciprocal lattice vectors
    b1 = np.cross(a2, a3) / V
    b2 = np.cross(a3, a1) / V
    b3 = np.cross(a1, a2) / V
    
    #k grid
    grid = np.zeros((n, n, n, 3))
    
    for i in range(n):
        for j in range(n):
            for k in range(n):
                grid[i, j, k, :] = (i / n) * b1 + (j / n) * b2 + (k / n) * b3
                
    return grid

def readPotential(file):
    """"Reads potential, formatted in the form of a gaussian file, into an array"""

    
    with open(file, 'r') as file:
        file.readline()
        file.readline()
    
        natoms = int(file.readline().lstrip().split()[0])

        #grid density
        line = file.readline().lstrip().split()
        nx = int(line[0])
        vx = np.array([line[1], line[2], line[3]])
        
        line = file.readline().lstrip().split()
        ny = int(line[0])
        vy = np.array([line[1], line[2], line[3]])

        line = file.readline().lstrip().split()
        nz = int(line[0])
        vz = np.array([line[1], line[2], line[3]])

        for i in range(natoms):
            file.readline()
            
        data = np.zeros( (nx, ny, nz) )

        for i in range(nx):
        
            for j in range(ny):
            
                k = 0
                while (k < nz):
                
                    line = file.readline().split()

                   
                    for value in line:
                        
                        data[i, j, k] = float(value)
                        k += 1
        
    return (data, vx, vy, vz)


def bandfill(file, VBM, CBM):
    """Computes the correction term due moss burnstein type band filling which occures due to unphysically large defect concentrations.
       Computes the correction for shallow donors and acceptros
       
       Host VBM and CBM values after potential alignment
       file: Name of file containing as columns the occupations, energy eigenvalues and wights"""

    #load energies (eV) occupations and weights
    data = np.loadtxt(file, skiprows=1)
    energies = data[:, 3]
    occupations = data[:, 4]
    weights = data[:, 5]


    #shallow donor correction
    sd = -1 * weights * occupations * (energies - CBM)
   
    #shallow acceptor correction
    sd = np.sum(sd[energies > CBM])

    sa = weights * (1 - occupations) * (energies - VBM)
    sa = np.sum(sa[energies < VBM])



    print("Shallow donor: {}\nShallow acceptor: {}".format(sd, sa))
    
        

class cell(object):
    """unit cell object. Attributes: 
                         1. lattice: 3x3 matrix storing lattice vectros as its columns
                         2. atoms: list of atom object, specifying the atomic basis of the unit cell
                                   each atom object stores the name and position of the atom

                         Contains methods for unit cell manipulation and calculating corrections to the formation energy in the supercell aproach
    """
    
    def __init__(self, file='', atoms=[], lattice=[]):
        """2 ways to initiallize a lattice object: 1. specify the lattice matrix and list of atom objects
                                                   2. specify a file with format: alat
                                                                                   v1(1) v1(2) v1(3)
                                                                                   v2(1) v2(2) v2(3)
                                                                                   v3(1) v3(2) v3(3)
                                                                                   A x1(1) x2(1) x3(1)
                                                                                   B x2(1) x2(2) x2(3)
                                                                                      ...........
                                                                                      ...........
                                                      First line containes the lattice parameter, the next 3 the lattice vectors, 
                                                      the rest of the file should specify the atom name and position
      """
        
        if not len(atoms) or not len(lattice):
            
            with open(file, 'r') as file:

                #lattice parameter is the first line of file
                self.alat = float( file.readline() )
                self.lattice = np.zeros( (3,3) )
                self.atoms = []
                
                #next 3 lines specify the lattice vectors
                for i in range(3):
                    line = file.readline()
                    line = line.split()
                    self.lattice[:, i] = np.array([float(line[0]), float(line[1]), float(line[2])])


                #rest of the file specifies atom names and their positions in the lattice
                line = file.readline()
                    
                while(line):
                    fields = line.split()
                    pos = [float(x) for x in fields[1:]]
                    self.atoms.append( atom( fields[0], np.array(pos) ) )

                    line = file.readline()
                        
        else:

            self.alat=0
            self.atoms = atoms
            self.lattice = lattice
                        


    def NN(self, pos):
        """Returns coordeinates of atoms closest to the point pos in space"""
        pass

    def prim(cell):
        """returns position of atoms in the primitive cell"""

        newAtoms = []
        newPos = np.array(3)
        
        for at in cell.atoms:

            #position in crystal coordinates
            newPos[:] = dot(inv(cell.lattice), at.pos)
            
            
            for i in range(3):
                while True:
                    
                    if at.pos[i] >=1:
                        newPos[i] -=1
                        
                    elif at.pos[i] < 0:
                        
                        newPos[i] += 1

                    else:
                        break
                    
            newAtoms.append(atom(at.name, newPos))


            return newAtoms
    
    def RtoO(cel):
        """Converts a rhombohedral primitive lattice to an orthorhombic cell
          !!!!Assumes atomic positions are in crystal coordinates!!!!"""
    
        #matrix of the rhombohedral lattice vectors
        lattice = cel.lattice

        #stores hexagonal lattice vectors
        lattice_H = np.empty( (3,3) )

        #volume of rhombohedral cell
        V_R = np.dot(lattice[:, 0], np.cross(lattice[:, 1], lattice[:, 2]))

        #length of c axis in hexagonal system
        c0 = norm( lattice[:, 0] + lattice[:, 1] + lattice[:, 2] )

        #length of a axis in hexagonal
        a0 = np.sqrt( 6 * V_R / (np.sqrt(3) * c0) )

        #The hexagonal vectors
        lattice_H[:, 0] = [a0, 0, 0]
        lattice_H[:, 1] = [a0 / 2, a0 * np.sqrt(3) / 2, 0]
        lattice_H[:, 2] = [0, 0, c0]

        #Orthorombic lattice vectors
        lattice_O = np.zeros((3,3))

        lattice_O[:, 0] = lattice_H[:, 0]
        lattice_O[:, 1] = [0, lattice_H[0, 0] * np.sqrt(3), 0]
        lattice_O[:, 2] = lattice_H[:, 2]

        superAtoms = []
        #create rhobohedral supercell

        n = 4
        for i in range(-n, n):
            for j in range(-n, n):
                for k in range(-n, n):
                    for at in cel.atoms:
                        superAtoms.append( atom(at.name, at.pos + np.array( [i, j, k] )) )
                    
                    
        atomsO = []
        #change atomic coordinates to orthorombic basis

        m = np.dot(inv(lattice_O), lattice) #change of basis matrix
        atomsO = [atom(at.name, np.dot(m, at.pos)) for at in superAtoms]

        #include only those atoms with coordinates between 0 and 1
        atomsO = [at for at in atomsO if np.sum( at.pos >= 0 ) == 3 and np.sum( at.pos < 0.9999 ) == 3 ]
    
        return cell(0, atomsO, lattice_O)
                
    def super(self, n, l, m):
        """Creates a super cell by scalling each primitive lattice vector by n, l, m"""
        lattice = self.lattice

        lattice[:, 0] *= n
        lattice[:, 1] *= l
        lattice[:, 2] *= m
        
        for at in self.atoms:
            at.pos /= np.array([n, l, m])
            
        #assumes atoms are in crystal coordinates
        atoms = []
        for i in range(n):
            for j in range(l):
                for k in range(m):
                    for at in self.atoms:
                        
                        pos = np.array([(i / n),(j / l), (k / m)]) + at.pos
                        atoms.append(atom(at.name, pos))
        self.atoms = atoms

    def __str__(self):
        """Prints to stdout in the same format as the input file for instantiating object"""
        
        string = "{}\n".format(self.alat)

        for i in range(3):
            string += str(self.lattice[:, i])[1:-1] + "\n"

        for atom in self.atoms:
            string += "{} {} {} {}\n".format(atom.name, atom.pos[0], atom.pos[1], atom.pos[2])
                        
        return string
    
"""This section contains a series of methods for computing the image charge correction"""

    def realSpace(self,convP, eps=np.identity(3)):
        """Computes the real space part of the Ewald sum
           Parameters: convP = Ewald convergence parameter
                        eps   = Dielectric Tensor (if material is isotropic provide the full tensor i.e. eps * np.identiy(3))

        """
        
        ieps = inv(eps)
        deps = det(eps)

        lattice = self.lattice
        sumDir = 0
        sConvP  = np.sqrt(convP)

        #Limits the terms of the sum to less than  erfc(15) (beyond which the terms of the sum are negligible)
        iRange  = int( 15 /  (sConvP * norm(lattice[:, 0])) ) + 1
        jRange =  int( 15 /  (sConvP * norm(lattice[:, 1])) ) + 1
        kRange =  int( 15 /  (sConvP * norm(lattice[:, 2])) ) + 1


        for i in range(-iRange, iRange + 1):
                   
            for j in range(-jRange, jRange + 1):
                   
                for k in range(-kRange, kRange + 1):
                   
                   if i == j == k == 0:
                       continue
                   
                   R = i * lattice[:, 0] + j * lattice[:, 1] + k * lattice[:, 2]
                   arg = np.sqrt( dot(R, dot(ieps, R)) )

                   sumDir += erfc( sConvP * arg) / arg / np.sqrt(deps)  
    
        return  sumDir * norm(lattice[:, 0])

    def reciprocal(self, thr, convP, eps=np.identity(3)):
        """Computes the reciprocal space part of the ewald sum
           Parameters: thr   = Energy cut off in eV (Tune this wrt the conergence parameter to find a reasonable value)
                       convP = Ewald convergence Parameter
                       eps   = Dielectric Tensor
        """

        
        lattice = self.lattice
        V = np.dot(lattice[:, 0], np.cross(lattice[:, 1], lattice[:, 2]))

        #reciprocal lattice vectors
        rec = np.empty((3,3))
        rec[:, 0] = 2 * np.pi * np.cross(lattice[:, 1], lattice[:, 2]) / V
        rec[:, 1] = 2 * np.pi * np.cross(lattice[:, 2], lattice[:, 0]) / V
        rec[:, 2] = 2 * np.pi * np.cross(lattice[:, 0], lattice[:, 1]) / V

        sumRec = 0
        thr = thr / ( 13.605698066 * 2 ) #energy in Hartree atomic units

        #Limits the terms of the sum to the energy cutoff provided by thr
        iRange = int( np.sqrt(2 * thr) / norm(lattice[:, 0]) ) + 1
        jRange = int( np.sqrt(2 * thr) / norm(lattice[:, 1]) ) + 1
        kRange = int( np.sqrt(2 * thr) / norm(lattice[:, 2]) ) + 1
        
        for i in range(-iRange, iRange + 1):

            for j in range(-jRange, jRange + 1):
        
                for k in range(-kRange, kRange + 1):
                    
                    if i == j == k == 0:
                        continue
                    
                    G = i * rec[:, 0] + j * rec[:, 1] + k * rec[:, 2]

                    arg = dot(G, dot(eps, G))
                    normG = norm(G)
                    sumRec += 4 * np.pi / V  * np.exp( -arg / (4 * convP)) / arg


        return sumRec

    def ewaldPot(self, thr, convP, charge, eps=np.identity(3)):
        """Ewald Potential obtained by adding the real and reciprical space contributions"""
        
        deps = det(eps)
        lattice = self.lattice
        V = np.dot(lattice[:, 0], np.cross(lattice[:, 1], lattice[:, 2]))
        rsum = self.realSpace(convP, eps)
        gsum = self.reciprocal(thr, convP, eps)
        sconvP = np.sqrt(convP)

        
        pot = rsum + gsum - np.pi / (V * convP) - (2 * sconvP) / np.sqrt(np.pi * deps)
        return pot * charge

    def firstO(self, thr, convP, charge,eps=np.identity(3)):
        """Screened coloumb energy of a periodic array of point charges in a uniform neutralizing background, obtained by integrating the ewaldPot"""
        return self.ewaldPot(thr, convP, charge, eps) * -1 * charge / 2 


    def integrate(self, size=200):
        """Integrate over a weigner seitz cell of the lattice using simple Reinman sum
        size = grid density for integral (size x size x size)
        Sampeling center of voxels
        
        Note: This is very slow; Better to vectorize
        """
        
        lattice = self.lattice
    

        #cartesian Tensor
        g = np.dot(lattice.transpose(), lattice)

        #volume of cell
        V = np.sqrt(det(g))
        
        #begin integration (sum)

        sum = 0
        for l in range(size):
            for m in range(size):
                for n in range(size):

                    candidate_integrand = []
                    #finds the position between all surrounding cells that is closest
                    #to the origin i.e. insde the weigner seitz cell
                    for i in [-1, 0, 1]:
                        for j in [-1, 0, 1]:
                            for k in [-1, 0, 1]:
                                latV = np.array([l, m, n]) / size #lattice vector
                                pos = np.array([i, j, k]) #postion inside cell
                            
                                vector = pos + latV - 0.5
                            
                                length = dot(vector, dot(g, vector))
                            
                                candidate_integrand.append(length)

                    #the smallest coordinate is in the weigner seitz cell
                    sum += np.min(candidate_integrand)
                
        #integration is done in crystal coordinate (1/V = det(Jacobian))
        return  sum / (V * size ** 3)

    def correction1(self, charge, eps, thr, convP):
        
        thirdOrder = self.integrate(size=100)* (4 * np.pi / 3 * charge ** 2) * (1 - 1 / eps) / eps
        firstOrder = self.firstO(thr, convP, charge, eps * np.identity(3))
    
        return (firstOrder - thirdOrder) * (27.211396)

    def shapeFactor(self, thr, convP, size=200):
        """Shape factor for calcualtion of the image charge correction"""
        
        firstO = self.firstO(thr, convP, 1, eps=np.identity(3))
        
        thirdO = self.integrate(size) * -1 * (2 * np.pi / 3)

        
        return thirdO / firstO

    def LZ_img(self, thr, convP, charge, eps):
        """Expects all quantities in Hartree atomic units (hbar = e = me = 1)"""
        """Lany and Zunger correction to the spurious electrostatic intraction between defects
           Output: Image charge correction in eV"""
        return self.firstO(thr, convP, charge, eps * np.identity(3)) * (1 + self.shapeFactor(thr, convP, 100) * (1 - 1 / eps)) * 27.211396
        
"""Methods for performing atomic sphere averaging and potential alignment"""

    def atomSphereAvg(self, potFile, radii):
        """Atom sphere averaged potential difference between host and defect potentials (PP + Hartree)
           Parameters: potFile: PP + vHartree as a cube file from pp.x output
                        radii: Radius of the sphere about each atom to the averaging at"""
    
        Pot, vx, vy, vz  = readPotential(potFile)

        #Sphere average host potential at atomic positions for host
        spherePots = []
    
        lenTheta = 50
        lenPhi = 50
        vRad = np.zeros(3)
        thetas = np.linspace(0, np.pi, lenTheta)
        phis = np.linspace(0, 2 * np.pi, lenPhi)
        
        for at, normR in zip(self.atoms, radii):

            #generate points on the surface of a sphere
            radAvg = 0
            for theta in thetas:
                for phi in phis:

                    vRad[0] = normR * np.sin(theta) * np.cos(phi)
                    vRad[1] = normR * np.sin(theta) * np.sin(phi)
                    vRad[2] = normR * np.cos(theta)

                    #point on the sphere about the atom in cartesian coordinates
                    pos = np.dot(self.lattice, at.pos) + vRad
                    pos = np.array(pos, dtype='float')

                    #find approximate postion of sphere point interms of array index
                    ipos = np.dot(inv(np.array([vx, vy, vz], dtype='float')), pos)
                    ipos = np.array(list(map(int, ipos)))

                    radAvg += Pot[tuple(ipos)]

            spherePots.append( radAvg  / (lenTheta * lenPhi) )
            

        return spherePots


    def potAllign(host, defect, cellH, cellD, radH, radD, thr, defectPos=np.array([0,0,0])):
        """ Assumes defect is the first element of the atom list"""
        
        """ Currently planar averages
        1. Average potential over same atom type
        2. take difference between average and actual
        3. Exclude atom if difference is beyond threshold (want potential at atomic sites far away from defect)
        4. Take the difference between potentail of defect and host at accepetible atomic sites then average for alignment
        """
        hostPot = cellH.atomSphereAvg(host, radH)
        defectPot = cellD.atomSphereAvg(defect, radD)

        
        #Not considering potential at defect site in average
        if len(cellH.atoms) > len(cellD.atoms):
            print("ppoooop")
            atomsD = cellD.atoms
            atomsH = cellH.atoms[1:]
            hostPot = hostPot[1:]

        elif len(cellH.atoms) < len(cellD.atoms):

            atomsD = cellD.atoms[1:]
            defectPot = defectPot[1:]
            atomsH = cellH.atoms

        else:

            atomsD = cellD.atoms[1:]
            defectPot = defectPot[1:]
            
            atomsH = cellH.atoms[1:]
            hostPot = hostPot[1:]

        lattice = cellD.lattice
        atomPos = []
        print(len(defectPot), len(hostPot))
        #compute distance to defect
        for at in atomsD:
            positions = []
            for i in [-1, 0, 1]:
                for j in [-1, 0, 1]:
                    for k in [-1, 0, 1]:
                        R = i * lattice[:, 0] + j * lattice[:, 1] + k * lattice[:, 2]
                        positions.append( norm(at.pos - (defectPos + R)) )
                        
            atomPos.append(min(positions))
                        
            

        
        avgPot = np.average(defectPot)
        std = np.std(defectPot)
    
        include = list( range(len(defectPot)) )
        
        #include only those atoms in defect cell whose potential is within 1 std of the average defect potential
        copyDefectPot = list(defectPot)

        
        #average the potential for defects at equivalent sites
        #find positions with the same distance from the defect
        sames = []
        for i in range(len(atomPos)):
            same = []
            for j in range(len(atomPos)):

                if abs(atomPos[i] - atomPos[j]) < 0.0001:
                    same.append(j)
            sames.append(same)

        sames = [tuple( sorted(i)  ) for i in sames]
        sames = set(sames)
        sames = list(map(list, sames))
        print(sames)
        
        difference = np.array(defectPot) - np.array(hostPot)
        newPot = []
        newDist = []
        for same in sames:
            newPot.append(np.average(difference[same]))
            newDist.append(atomPos[same[0]])
            
        newPot = np.array(newPot)
        newDist = np.array(newDist)
            
        
        
        import matplotlib.pyplot as plt
        
        #plot against distance to defect
        plt.plot(newDist, newPot, 'o')
        plt.show()
       



class atom(cell):
    """Atom object. Attributes: 1.name: Element name
                                2.pos: position of atom in unit cell"""
    
    def __init__(self, name, pos):
        self.name = name
        self.pos = pos

    def __str__(self):
        return  "{}\t{} {} {}".format(self.name, *self.pos)



