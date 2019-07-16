# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 15:15:39 2019

@author: Ali
"""

import numpy as np
from scipy.special import erfc 

        
def birch(V, V0, K0, K1):
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
        
    return data
        
                
def potAlign(Host, Defect, tol):
    """input: 
        Defect, Host: Local potential from defect and host
        tol: Energy tolerance
        
        1. Average potential over same atom type
        2. take difference between average and actual
        3. Exclude atom if difference is beyond threshold (want potential at atomic sites far away from defect)
        4. Take the difference between potentail of defect and host at accepetible atomic sites
        
    """
    hostPot = readPotential(Host)
    defectPot = readPotential(Defect)
    
    hostAvg = np.average(hostPot, axis=(1,2))
    defectAvg = np.average(defectPot, axis=(1,2))

    return defectAvg - hostAvg
    
    
def imgChage(SC, charge, medalung, L, epsil):
    """Computes the correction term due to unphysical electrostatic interaction between periodic images of charged defects in the super cell method
        see  M Leslie and N J Gillan 1985 J. Phys. C: Solid State Phys. 18 973
        and  Stephan Lany and Alex Zunger 2009 Modelling Simul. Mater. Sci. Eng. 17 084002
        for more details
       The form of the correction used is the one specified in the second paper
        
        input: SC: Shape factor
               epsil: dielectric constant of prestine cell
               L: Linear supercell dimenstion (V ** 1/3)
               Madelung: Madelung constant of lattice"""
    
    bohr_A = 1 / 1.88973
    charge = 1.602e-19 #C
    e0 = 8.8541878e-12 #F/m

    #convert from bohr to meter
    L *= bhor_A * 1e-10

    return  ( 1 + SC * (1 - 1 / epsil) ) * (charge) * medalung / (2 * epsil * L * e0)

def bandfill(file, fermi):
    """Computes the correction term due moss burnstein type band filling which occures due to unphysically large defect concentrations.
       Computes the correction for shallow donors and acceptros
       
       fermi : any energy in the band gap will do
       file: Name of file containing as columns the occupations, energy eigenvalues and wights"""
       
    data = np.loadtxt(file, skiprows=1)
       
    energies = data[:, 3]
    occupations = data[:, 4]
    weights = data[:, 5]
      
       
    #Find CBM, VBM
    #fermi needs to be inside gap for this to work
    CBM = min([x for x in energies if x > fermi])
    VBM = max([x for x in energies if x < fermi])

    sd = weights * occupations * (energies - CBM)
    sd = np.sum(sd[energies > CBM]) / np.sum(weights)

    sa = weights * (1 - occupations) * (energies - VBM)
    sa = np.sum(sa[energies < VBM]) / np.sum(weights)

    print("Shallow donor: {}\nShallow acceptor: {}".format(sd, sa))
    
        
    
    


    
class cell(object):


    
    def __init__(self, file, atoms=[], lattice=[]):
        """Reads a structure file and extracts out the atoms, lattice vectors and lattice parameter"""
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
                        


    def NN():
        pass

    def prim(cell):
        """returns position of atoms in the primitive cell"""

        #convert to crystal coordinates
        #if negative add 1
        #if grater than or equal to 1 subtract 1

        #change of basis matrix from cartesian to crystal
        newAtoms = []
        newPos = np.array(3)
        
        for at in cell.atoms:

            #position in crystal coordinates
            newPos[:] = np.dot(np.inv(cell.lattice), at.pos)
            
            
            for i in range(3):
                while True:
                    
                    if at.pos[i] >= 1:
                        newPos[i] -=1
                        
                    elif at.pos[i] < 0:
                        
                        newPos[i] += 1

                    else:
                        break
                    
            newAtoms.append(atom(at.name, newPos))


            return newAtoms
    
    def RtoO(cel):
        """Converts a rhombohedral primitive lattice to an orthorhombic lattice"""
        ###Atoms should be in cartesian coordinates, multiply lattice vectors by alat
        #bring into primitive cell
    
        #matrix of the hexagonal lattice vectors

        lattice = cel.lattice

        #stores hexagonal lattice vectors
        lattice_H = np.empty( (3,3) )

        #volume of rhombohedral cell
        V_R = np.dot(lattice[:, 0], np.cross(lattice[:, 1], lattice[:, 2]))

        #length of c axis in hexagonal system
        c0 = np.linalg.norm( lattice[:, 0] + lattice[:, 1] + lattice[:, 2] )

        #length of a axis in hexagonal
        a0 = np.sqrt( 6 * V_R / (np.sqrt(3) * c0) )

        #The hexagonal vectors
        lattice_H[:, 0] = [a0, 0, 0]
        lattice_H[:, 1] = [a0 / 2, a0 * np.sqrt(3) / 2, 0]
        lattice_H[:, 2] = [0, 0, c0]

        #Basis for the hexagonal system (contains 3 lattice points instead of 1)
        atoms_H = []
        n = 5
        for at in cel.atoms:
            atoms_H.append(atom(at.name, at.pos + [0, 0, 0]))
            atoms_H.append(atom(at.name, at.pos + lattice[:, 0]))
            atoms_H.append(atom(at.name, at.pos + lattice[:, 1]))
        
        #Orthorombic lattice vectors
        lattice_O = np.zeros((3,3))

        lattice_O[:, 0] = lattice_H[:, 0]
        lattice_O[:, 1] = [0, lattice_H[0, 0] * np.sqrt(3), 0]
        lattice_O[:, 2] = lattice_H[:, 2]


        atoms_O = []

        
        for at in atoms_H:
            atoms_O.append(atom(at.name, at.pos + [0, 0, 0]))
            atoms_O.append(atom(at.name, at.pos + lattice_H[:, 1]))               

        print("h1: {}", lattice_H[:, 1])
        return cell(0, atoms_O, lattice_O)
                
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
        string = "{}\n".format(self.alat)

        for i in range(3):
            string += str(self.lattice[:, i])[1:-1] + "\n"

        for atom in self.atoms:
            string += "{} {} {} {}\n".format(atom.name, atom.pos[0], atom.pos[1], atom.pos[2])
                        
        return string
            
    def makeCube():
      
        #cartesian basis
        cube = np.array([ [1, 0, 0],
                          [0, 1, 0],
                          [0, 0, 1] ])
        
        #change of coordinates matrix from lattice to cartesian
        change = np.dot(np.inv(cube), self.lattice) 
    
                

                

class atom(cell):
    
    def __init__(self, name, pos):
        self.name = name
        self.pos = pos

    def __str__(self):
        return  "{}\t{} {} {}".format(self.name, *self.pos)

def uniqSites(cell):
    """Identifies wychoff sites in the super cell, used as condidate sites for vacancies and interstitials
    
       1.Find space group
       2. Perform symmetry operations on atomic basis
       3. Identify sites equivalent through symmetry within tolerance
       4. Return a set of symmetry inequivalent atoms"""
    
def Ewald(cell, n, convP):
    """Ewald sum for computing the Madelung energy of a periodic array of point charges neutralized        by a uniform background
       if EPS is a tensor, compute anisotropic form of sum
    """

    lattice = cel.lattice()
    #Direct lattice volume
    V = np.dot(lattice[:, 0], np.cross(lattice[:, 1], lattice[:, 2]))

    #reciprocal lattice vectors
    rec = np.empty((3,3))
    rec[:, 0] = np.cross(lattice[:, 1], lattice[:, 2]) / V
    rec[:, 1] = np.cross(lattice[:, 2], lattice[:, 0]) / V
    rec[:, 2] = np.cross(lattice[:, 0], lattice[:, 1]) / V

   
    #reciprocal lattice sum
    sumRec = 0

    for i in range(n):
        for j in range(n):
            for k in range(n):

                if i != 0 and j != 0 and k != 0:
                    
                    G = i * rec[:, 0] + j * rec[:, 1] + k * rec[:, 2]
                    normG = np.linalg.norm(G)
                    sumRec += 4 * np.pi * np.exp( -(normG ** 2) / (4 * convP ** 2) ) / (V * normG ** 2)

    #direct lattice sum
    sumDir = 0

    for i in range(n):
        for j in range(n):
            for k in range(n):

                if i != 0 and j != 0 and k != 0:
                    
                    R = i * lattice[:, 0] + j * lattice[:, 1] + k * lattice[:, 2]
                    normR = np.linalg.norm(R)
                    sumDir += erfc( convP * normR) / normR  

    pot = (sumRec + sumDir - np.pi / (V  * convP ** 2) - 2 * convP / np.sqrt(np.pi)) * q / eps
    energy = -pot * q / 2 
    return pot, energy
    
"""Look into 2d projections"""
"""Python program for visualizing lattice and selecting defect coordinate"""
"""Python program for identifying wyckoff sites"""
