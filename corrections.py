# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 15:15:39 2019

@author: Ali
"""

import numpy as np


        
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

    data = np.zeros( (nx, ny, nz) )

    for i in range(nx):
        
        for j in range(ny):
            
            k = 0
            while (k < nz):
                
                line = file.readline().split()

                for value in line:
                    
                    data[i, j, k] = float(value)
                    k += 1
        
    return potential
        
                
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

    #3 x 3 matrix of the lattice vectors
    lattice = np.zeros((3, 3))

    alat = 0
    
    #Stores the atom objects of the cell
    atoms = []
    
    
    def __init__(file):
        """Reads a structure file and extracts out the atoms, lattice vectors and lattice parameter"""

        with open(file) as file:

            #lattice parameter is the first line of file
            alat = float( file.readline() )

            #next 3 lines specify the lattice vectors
            for i in range(3):
                line = file.readline()
                line = line.split(" ")
                lattice[:, i] = np.array([line[0], line[1], line[2]])

            #rest of the file specifies atom names and their positions in the lattice
            line = readline()
            while(line):
                fields = line.split(" ")
                pos = [float(x) for x in fields[1:]]
                
                atoms.append( atom( fields[0], np.array(fileds) ) )

                line = readline()


    def NN():
        pass
        
                
    def super(self, n, l, m):
        """Creates a super cell by scalling each primitive lattice vector by n, l, m"""

        lattice[:, 0] *= n
        lattice[:, 1] *= l
        lattice[:, 2] *= m
        
        for i in range(n):
            for j in range(l):
                for k in range(m):
                    for atom in range(len(atoms)):
                        pos = (i / n) * lattice[:, 0] + (j / l) * lattice[:, 1] + (k / m) * lattice[:, 2] + atom.pos

    def __str__(self):
        string = "{}\n".format(alat)

        for i in range(3):
            string += str(lattice[:, i])[1:-1] + "\n"

        for atom in atoms:
            string += "{} {} {} {}\n".format(atom.name, atom.pos[0], atom.pos[1], atom.pos[2])
                        
        return string
            
    def makeCube():
        """NVM not so easy :(, need to know crystal symmetries"""
        """Changes atomic basis to cubic for to simplify super cell calculations"""

        #cartesian basis
        cube = np.array([ [1, 0, 0],
                          [0, 1, 0],
                          [0, 0, 1] ])
        
        #change of coordinates matrix from lattice to cartesian
        change = np.dot(np.inv(cube), self.lattice) 
    
                

                

class atom(cell):
        
    name = "atom"
    pos = np.zeros(3)

    def __init__(self, name, pos):
        self.name = name
        self.pos = pos
        
"""Look into 2d projections"""
"""Python program for visualizing lattice and selecting defect coordinate"""
"""Python program for identifying wyckoff sites"""
