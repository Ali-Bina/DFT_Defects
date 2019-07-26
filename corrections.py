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
        
                
def potAlign(Host, Defect, tol):
    """input: 
        Defect, Host: Local potential from defect and host
        tol: Energy tolerance
        Currently planar averages
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

def bandfill(file, VBM, CBM):
    """Computes the correction term due moss burnstein type band filling which occures due to unphysically large defect concentrations.
       Computes the correction for shallow donors and acceptros
       
       Host VBM and CBM values after potential alignment
       file: Name of file containing as columns the occupations, energy eigenvalues and wights"""
       
    data = np.loadtxt(file, skiprows=1)
       
    energies = data[:, 3]
    occupations = data[:, 4]
    weights = data[:, 5]
    
    sd = -1 * weights * occupations * (energies - CBM)
    sd = np.sum(sd[energies > CBM])

    sa = weights * (1 - occupations) * (energies - VBM)
    sa = np.sum(sa[energies < VBM])



    print("Shallow donor: {}\nShallow acceptor: {}".format(sd, sa))
    
        

    
class cell(object):

    def __init__(self, file='', atoms=[], lattice=[]):
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
        string = "{}\n".format(self.alat)

        for i in range(3):
            string += str(self.lattice[:, i])[1:-1] + "\n"

        for atom in self.atoms:
            string += "{} {} {} {}\n".format(atom.name, atom.pos[0], atom.pos[1], atom.pos[2])
                        
        return string
###############################EWALD############################
    def realSpace(self,convP, eps=np.identity(3)):

        #eps = dielectric Tensor
        ieps = inv(eps)
        deps = det(eps)

        lattice = self.lattice
        sumDir = 0
        sConvP  = np.sqrt(convP)

        #ensures that sum does not exceed erfc(15)
        iRange  = int( 15 /  (np.sqrt(convP) * norm(lattice[:, 0])) ) + 1
        jRange =  int( 15 /  (np.sqrt(convP) * norm(lattice[:, 1])) ) + 1
        kRange =  int( 15 /  (np.sqrt(convP) * norm(lattice[:, 2])) ) + 1


        for i in range(-iRange, iRange + 1):
                   
            for j in range(-jRange, jRange + 1):
                   
                for k in range(-kRange, kRange + 1):
                   
                   if i == j == k == 0:
                       continue
                   
                   R = i * lattice[:, 0] + j * lattice[:, 1] + k * lattice[:, 2]
                   arg = np.sqrt( dot(R, dot(ieps, R)) )

                   sumDir += erfc( sConvP * arg) / arg / np.sqrt(deps)  
    
        return  sumDir * norm(lattice[:, 0])

    def reciprocal(self, thr, convP, eps=np.identity(3), r=np.zeros(3)):
        lattice = self.lattice
        V = np.dot(lattice[:, 0], np.cross(lattice[:, 1], lattice[:, 2]))
    
        rec = np.empty((3,3))
        rec[:, 0] = 2 * np.pi * np.cross(lattice[:, 1], lattice[:, 2]) / V
        rec[:, 1] = 2 * np.pi * np.cross(lattice[:, 2], lattice[:, 0]) / V
        rec[:, 2] = 2 * np.pi * np.cross(lattice[:, 0], lattice[:, 1]) / V

        sumRec = 0

        thr = thr / (13.605698066 * 2 ) #energy in Hartree
        iRange = int( np.sqrt(thr) / norm(lattice[:, 0]) ) + 1
        jRange = int( np.sqrt(thr) / norm(lattice[:, 1]) ) + 1
        kRange = int( np.sqrt(thr) / norm(lattice[:, 2]) ) + 1
        
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
        deps = det(eps)
        lattice = self.lattice
        V = np.dot(lattice[:, 0], np.cross(lattice[:, 1], lattice[:, 2]))
        rsum = self.realSpace(convP, eps)
        gsum = self.reciprocal(thr, convP, eps)
        sconvP = np.sqrt(convP)
                   
        pot = rsum + gsum - np.pi / (V * convP) - (2 * sconvP) / np.sqrt(np.pi * deps)
        return pot * charge

    def firstO(self, thr, convP, charge,eps=np.identity(3)):
        """Screened coloumb energy of a periodic array of point charges in a uniform neutralizing background """
        return self.ewaldPot(thr, convP, charge, eps) * -1 * charge / 2 


    def integrate(self, size=100):
        """Integrate over a weiner seitz cell of the lattice
        size = grid size
        Sampeling center of voxels
        dx = """
        
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
                
        #integration is done in crystall coordinate (1/V = det(Jacobian))
        return  sum / (V * size ** 3)

    def correction1(charge, eps, thr, convP):
    
        thirdOrder = integrate(size=100)* (4 * np.pi / 3 * charge ** 2) * (1 - 1 / eps) / eps
        firstOrder = firstO(thr, convP, charge, eps * np.identity(3))
    
        return firstOrder - thirdOrder

    def shapeFactor(self, thr, convP):
        E1 = self.firstO(thr, convP, 1, eps=np.identity(3))
        
        E3 = self.integrate(size=50) * -1 * (2 * np.pi / 3)

        
        return E3 / E1

    def LZ_img(self, thr, convP, charge, eps):
        """Expects all quantities in Hartree atomic units (hbar = e = me = 1)"""
        """Lany and Zunger correction to the spurious electrostatic intraction between defects"""
        return self.firstO(thr, convP, charge, eps * np.identity(3)) * (1 + self.shapeFactor(thr, convP) * (1 - 1 / eps)) * 27.211396
        
    ###################Potential Allignment################

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


    def potAllign(host, defect, cellH, cellD, radii, thr):
        """ Assumes defect is the first element of the atom list"""
        
        hostPot = cellH.atomSphereAvg(host, radii)
        defectPot = cellD.atomSphereAvg(defect, radii)

        
        #Not considering potential at defect site in average
        if len(cellH.atoms) > len(cellD.atoms):

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
            
        avgPot = np.average(defectPot)
        std = np.std(defectPot)
    
        include = list( range(len(defectPot)) )
        
        #include only those atoms in defect cell whose potential is within 1 std of the average defect potential
        copyDefectPot = list(defectPot)

        defectPot2 = []
        hostPot2 = []
        for potH, potD in zip(defectPot, hostPot):

            if abs( (potD - avgPot) / std ) < 0.5:

                defectPot2.append(potD)
                hostPot2.append(potH)
        import matplotlib.pyplot as plt
        #plot against distance to defect
        plt.plot(np.array(defectPot2) - np.array(hostPot2))
        plt.show()
        return np.average( np.array(defectPot2) - np.array(hostPot2) )



class atom(cell):

    
    
    def __init__(self, name, pos, charge=0):
        self.name = name
        self.pos = pos
        self.charge=charge

    def __str__(self):
        return  "{}\t{} {} {}".format(self.name, *self.pos)



