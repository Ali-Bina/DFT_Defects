# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 15:15:39 2019

@author: Ali
"""
import math
from scipy.interpolate import RegularGridInterpolator as interpolate
import numpy as np
from scipy.special import erfc
from numpy.linalg import norm, inv, det
from numpy import dot
import matplotlib.pyplot as plt

def Murn(V, V0, K0, K1):
    """Murnaghan equation of state, fit to DFT total energy vs V and extract equillibrium lattice parameter and Bulk modulus"""
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
    #if occupations dont add to one maybe must devide by their sum, magnitude of PA is used
    """Computes the correction term due moss burnstein type band filling which occures due to unphysically large defect concentrations.
       Computes the correction for shallow donors and acceptros
       
       Host VBM and CBM values after potential alignment
       file: Name of file containing as columns the occupations, energy eigenvalues and wights"""

    #load energies (eV) occupations and weights
    data = np.loadtxt(file, skiprows=1)
    energies = data[:, 3]
    occupations = data[:, 4]
    weights = data[:, 5]

    print(np.sum(weights))
    #shallow donor correction
    sd = -1 * weights * occupations * (energies - CBM)
   
    #shallow acceptor correction
    sd = np.sum(sd[energies > CBM])

    sa = weights * (2 - occupations) * (energies - VBM)
    sa = np.sum(sa[energies < VBM])



    print("Shallow donor: {}\nShallow acceptor: {}".format(sd, sa))
    
        

class cell(object):
    """unit cell object. Attributes: 
                         1. lattice: 3x3 matrix storing lattice vectros as its columns
                         2. atoms: list of atom object, specifying the atomic basis of the unit cell
                                   each atom object stores the name and position of the atom

                         Contains methods for unit cell manipulation and calculating corrections to the formation energy in the supercell aproach
    """
    
    def __init__(self, file='', atoms=[], lattice=[], alat=1):
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
                lines = file.readlines()

                self.alat = 1
                
                #lattice parameter is the first line of file
                if len(lines[0].split()) == 1:
                    self.alat = float( lines[0][:-2] )
                    lines = lines[1:]
                    
                self.lattice = np.zeros( (3,3) )
                self.atoms = []
                
                #next 3 lines specify the lattice vectors
                for i in range(3):
                    line = lines[i][:-2]
                    line = line.split()
                    self.lattice[:, i] = np.array([float(line[0]), float(line[1]), float(line[2])])


                #rest of the file specifies atom names and their positions in the lattice
                lines = lines[3:]
                    
                for line in lines:
                    line = line[:-2]
                    fields = line.split()
                    pos = [float(x) for x in fields[1:]]
                    self.atoms.append( atom( fields[0], np.array(pos) ) )

                        
        else:

            self.alat=alat
            self.atoms = atoms
            self.lattice = lattice
                        

    def prim(position):
        """returns position of atoms in the primitive cell"""

        
        newPos = np.array( [round(pos, 6) for pos in position]  )
        
        for i in range(3):
            while True:
                
                if newPos[i] >=1:
                        
                    newPos[i] -=1
                        
                elif newPos[i] < 0:
                        
                    newPos[i] += 1

                else:
                    break
                
        return newPos
    
    def RtoO(cel):
        """Converts a rhombohedral primitive lattice to an orthorhombic cell
           input: cel object, of rhombohedral cell with atomic positions given in CRYSTAL COORDINATES

           output: New cell object for the orthorhombic super cell
  
"""
    
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
        """Creates a super cell by scalling each primitive lattice vector by n, l, m
           Assumes atomic positions are in crystal coordinates"""
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
        string=''
        if self.alat != 0:
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

        #Limits the terms of the sum to less than  erfc(10) (beyond which the terms of the sum are negligible)
        iRange  = int( 10 /  (sConvP * norm(lattice[:, 0])) ) + 1
        jRange =  int( 10 /  (sConvP * norm(lattice[:, 1])) ) + 1
        kRange =  int( 10 /  (sConvP * norm(lattice[:, 2])) ) + 1


        for i in range(-iRange, iRange + 1):
                   
            for j in range(-jRange, jRange + 1):
                   
                for k in range(-kRange, kRange + 1):
                   
                   if i == j == k == 0:
                       continue
                   
                   R = i * lattice[:, 0] + j * lattice[:, 1] + k * lattice[:, 2]
                   arg = np.sqrt( dot(R, dot(ieps, R)) )
                   #print(erfc( sConvP * arg) / arg / np.sqrt(deps))
                   sumDir += erfc( sConvP * arg) / arg / np.sqrt(deps)  
    
        return  sumDir * norm(lattice[:, 0])

    def reciprocal(self, thr, convP, eps=np.identity(3)):
        """Computes the reciprocal space part of the ewald sum
           Parameters: thr   = Energy cut off in eV (Tune this wrt the conergence parameter to find a reasonable value) for the reciprocal space sum
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
                    #print("Kspace: ", 4 * np.pi / V  * np.exp( -arg / (4 * convP)) / arg)
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

    def shapeFactor(self, thr, convP, size=200):
        """Shape factor for calcualtion of the image charge correction"""
        
        firstOrder = self.firstO(thr, convP, 1, eps=np.identity(3))
        
        thirdOrder = self.integrate(size) * -1 * (2 * np.pi / 3)

        
        return thirdOrder / firstOrder

    def LZ_img(self, thr, convP, charge, eps):
        """Expects all quantities in Hartree atomic units (hbar = e = me = 1), i.e. length in bohr"""
        """Lany and Zunger correction to the spurious electrostatic intraction between defects
           Output: Image charge correction in eV"""
        scalarEps = eps[0, 0]
        return self.firstO(thr, convP, charge, eps) * (1 + self.shapeFactor(thr, convP, 100) * (1 - 1 / scalarEps)) * 27.211396
        
    """Methods for performing atomic sphere averaging and potential alignment"""

    def atomSphereAvg(self, potFile, radii, n):
        """Atom sphere averaged potential difference between host and defect potentials (PP + Hartree)
           Parameters: potFile: PP + vHartree as a cube file from pp.x output
                        radii: Radius of the sphere about each atom to the averaging at
                            n: Grid density for integration """
    
        Pot, vx, vy, vz  = readPotential(potFile)
        print(vx, vy, vz, sep='\n')
        #Sphere average host potential at atomic positions for host
        spherePots = []
    
        lenTheta = n
        lenPhi = n
        lenRad = n
        vRad = np.zeros(3)
        thetas = np.linspace(0, np.pi, lenTheta)
        dtheta = thetas[1] - thetas[0]
        phis = np.linspace(0, 2 * np.pi, lenPhi)
        dphi = phis[1] - phis[0]

        x = np.linspace(0, 1, Pot.shape[0])
        y = np.linspace(0, 1, Pot.shape[1])
        z = np.linspace(0, 1, Pot.shape[2])
        #interpolate the potential grid
        pot_int = interpolate((x, y, z), Pot)
        
        for at, normR in zip(self.atoms, radii):
            rads = np.linspace(0, normR, lenRad)[1:]
            dR = rads[1] - rads[0]
            
            #Integration
            radAvg = 0
            for rad in rads:
                for theta in thetas:
                    for phi in phis:

                        vRad[0] = rad * np.sin(theta) * np.cos(phi)
                        vRad[1] = rad * np.sin(theta) * np.sin(phi)
                        vRad[2] = rad * np.cos(theta)

                        #point on the sphere about the atom in cartesian coordinates
                        pos = np.dot(self.lattice, at.pos) + vRad
                    
                        pos = np.dot( inv(self.lattice), pos  ) #back to crystal
                        pos = cell.prim(pos) #find position in unit cell
                        
                        #Add to riemann sum
                        radAvg += pot_int(pos) * np.sin(theta) * (rad ** 2) * dtheta * dphi * dR

            spherePots.append( radAvg  / (4 / 3 * np.pi * (normR ** 3)) )
            

        return spherePots

    def findDefect(self, host, defect, threshold = 0.1):
        """Finds and returns  the defect atom in the defect cell by comparing to host
           host: cell object for the host
           defect: cell object for the defect
           threshold: threshold for determining equivalent atoms in bohr (since atoms 
                      move after relaxation positions must be compared using a threshold)"""

        
        atomsD = defect.atoms
        atomsH = host.atoms

        atomCountH = dict()
        for at in atomsH:

            if at.name not in atomCountH.keys():

                atomCountH[at.name] = 0

            atomCountH[at.name] += 1
            

        atomCountD = dict()
        for at in atomsD:

            if at.name not in atomCountD.keys():

                atomCountD[at.name] = 0

            atomCountD[at.name] += 1

        numAtomsH = sum(atomCountH.values())
        numAtomsD = sum(atomCountD.values())
    
        defectType = ''
        Defect = ()
        if numAtomsH == numAtomsD:

            if atomCountH.keys() == atomCountD.keys():
                
                defectType = 'anti Site'
                    
            else:

                defectType = 'substitutional'

            for atH, atD in zip(atomsH, atomsD):

                if atH.name != atD.name:

                    #Tuple of defect type and defect position
                    Defect = ("{}_{}".format(atD.name, atH.name), atD.pos)
                
        elif numAtomsH > numAtomsD:
            defectType = 'Vacancy'

            for atH in atomsH:
                #print("--------------------------")
                for atD in atomsD:
                    #print(atD.pos - atH.pos)
                    #same atom found in host and defect (not vacancy)
                    if atH.name == atD.name and np.sum(abs(atD.pos - atH.pos) < threshold) == 3:
                    
                        break

                else:

                    #if the atom in host is not found in defect then it is the vacancy
                    Defect = ("V_{}".format(atH.name), atH.pos)
                    break

        elif numAtomsD > numAtomsH:
            defectType = 'interstitial'

            if atomCountH.keys() == atomCountD.keys():

                defectType = 'Self-interstitial'

                for atD in atomsD:
                    for atH in atomsH:

                        if atH.name == atD.name and np.sum(abs(atD.pos - atH.pos) < threshold) == 3:
                            
                            break

                        else:

                            Defect = (atD.name + "_i", atD.pos)
                            break
                

            else:

                setH = set(atomCountH.key())
                setD = set(atomCountD.keys())

                #extracts the extra defect atom not present in  host
                interstitial = str(setD - setH)[2:-2]
                
                Defect = (interstitial + "_i", atomCountD[interstitial])

        return Defect
        
    def potAlign(self, hostP, defectP, cellH, cellD, radH, radD):
        
        """ Currently planar averages

        host, defect: host and defect potential cube files
        cellH, cellD: host and defect cell objects
        radH, radD:  lists specifying the radius to use in the spherical averaging.
                     The lists must have the same number of elements as atoms in the defect cell

        Plots the sphere averaged potential of the host and defect cells as a function of distance to the defect and the difference.
        A slider on the plot can be used to specify how far from the defect to start averaging the difference.

        returns: a tuple with 4 elements containing the ploted information for further analyis. The first 2 elements are arrays containing atomic sphere averaged potential for host and corresponding distance from the defect. The last 2 elements contain the same information but for the defect cell
        """

        #identify position of defect in unit cell
        defect = self.findDefect(cellH, cellD)
        defectPos = defect[1]
        print("The defect is ", defect)

        #use defect type to modify host and defect atoms list

        thr = 0.001 #threshold for determining equivalent positions
        
        #if the defect is a vacancy
        if defect[0].split("_")[0] == 'V':

            index = 0
            for i, at in enumerate(cellH.atoms):
                
                if np.sum(abs(at.pos - defectPos) < thr) == 3:
                    
                    index = i

            cellH.atoms.pop(index)
            print("The index was ", index)
            
        #if the defect is an interstitial
        elif defect[0].split("_")[1] == 'i':

            index = 0
            for i, at in enumerate(cellD.atoms):
                
                if np.sum(abs(at.pos - defectPos) < thr):

                    index = i

            cellH.atoms.pop(index)
            cellD.atoms.pop(index)

        #substitutional defect
        elif defect[0].split("_")[0].isalpha():

            index = 0
            for i, at in enumerate(cellD.atoms):
                
                if np.sum(abs(at.pos - defectPos) < thr) == 3:

                    index = i

            cellD.atoms.pop(index)
            
        
        else:
            print("error no defects found")
            return



            
        #find sphere averaged potential of host and defect cells
        hostPot = cellH.atomSphereAvg(hostP, radH, 10)
        defectPot = cellD.atomSphereAvg(defectP, radD, 10)
        
        lattice = cellD.lattice * cellD.alat

        #Distance of atoms from defect site
        atomPosD = []
        atomPosH = []
        
        #compute distance to defect for Defect cell
        for at in cellD.atoms:
            positions = []
            for i in [-1, 0, 1]:
                for j in [-1, 0, 1]:
                    for k in [-1, 0, 1]:
                        
                        R = i * lattice[:, 0] + j * lattice[:, 1] + k * lattice[:, 2]
                        positions.append( norm(dot(lattice, at.pos + R - defectPos)) )
                        
            atomPosD.append(min(positions))

            
        #distance of host atoms to defect cell
        for at in cellH.atoms:
            positions = []
            for i in [-1, 0, 1]:
                for j in [-1, 0, 1]:
                    for k in [-1, 0, 1]:
                        
                        R = i * lattice[:, 0] + j * lattice[:, 1] + k * lattice[:, 2]
                        positions.append( norm(dot(lattice, at.pos + R - defectPos)) )
                        
            atomPosH.append(min(positions))
            
        from matplotlib.widgets import Slider

        atomPosD = np.array(atomPosD)
        atomPosH = np.array(atomPosH)
        defectPot = np.array(defectPot)
        hostPot = np.array(hostPot)

        diff = defectPot - hostPot
       
        #plots the defect and host potentials on one axis and the difference on another
        plt.subplot(121)
        p, = plt.plot(atomPosH, hostPot, 'o')
        plt.plot(atomPosD, defectPot, 'o')
        plt.title("Defect and host potentials", size=15)
        plt.xlabel("Distance From Defect (bohr)", size=15)
        plt.ylabel("Sphere averaged $V_{loc} + V_{H}$ (Ry)", size=15)
        plt.subplot(122)
        plt.title("Difference between host and defect potentials", size=15)
        plt.xlabel("Distance From Defect (bohr)", size=15)
        plt.plot(atomPosD, diff, 'o')
        plt.subplots_adjust(bottom=0.25)

        #Slider widget for specifying how far away from the defect to start averaging the difference potential
        AxSliderPos = plt.axes([0.1, 0.05, 0.8, 0.025])
        sliderPos = Slider(ax=AxSliderPos, label='Average From:', valmin=0, valmax=max(atomPosD), valinit=0)
        
        t = plt.gcf().text(0.02, 0.155, "Average: " + str(np.average(diff[atomPosD > 0])) + " Ry", fontsize=8)
        def updateText(val):
            avg = np.average(diff[atomPosD > sliderPos.val])

            #if slider goes beyond furthest atom, there is nothing to average and np.average returns nan
            if math.isnan(avg):
                avg = 0
                
            t.set_text("Average: " + str(avg) + " Ry")
            p.set_data(atomPosH, hostPot + avg)
            plt.draw()

        sliderPos.on_changed(updateText)
        plt.show()
       
        return (hostPot, atomPosH, defectPot, atomPosD)
        
       



class atom(cell):
    """Atom object. Attributes: 1.name: Element name
                                2.pos: position of atom in unit cell"""
    
    def __init__(self, name, pos):
        self.name = name
        self.pos = pos

    def __str__(self):
        return  "{}\t{} {} {}".format(self.name, *self.pos)



