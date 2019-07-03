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


def super(a1, a2, a3, n):
    """Creates a super cell by scalling each primitive lattice vector by n"""
    

                
                
def potAlign(Host, Defect, tol):
    """input: 
        Defect, Host: Local potential from defect and host
        tol: Energy tolerance
        
        1. Average potential over same atom type
        2. take difference between average and actual
        3. Exclude atom if difference is beyond threshold (want potential at atomic sites far away from defect)
        4. Take the difference between potentail of defect and host at accepetible atomic sites
        
    """
    
    
def imgChage(SC, charge, medalung, L, epsil):
    """Computes the correction term due to unphysical electrostatic interaction between periodic images of charged defects in the super cell method
        see  M Leslie and N J Gillan 1985 J. Phys. C: Solid State Phys. 18 973
        and  Stephan Lany and Alex Zunger 2009 Modelling Simul. Mater. Sci. Eng. 17 084002
        for more details
        
        input: SC: Shape factor
               epsil: dielectric constant of prestine cell
               L: Linear supercell dimenstion (V ** 1/3)
               Madelung: Madelung constant of lattice"""
    return SC * (1 - 1 / epsil) * (charge ** 2) * medalung / (2 * epsil * L)

def bandfill(file, Fermi):
    """Computes the correction term due moss burnstein type band filling which occures due to unphysically large defect concentrations.
       Computes the correction for shallow donors and acceptros
       
       file: Name of file containing as columns the occupations, energy eigenvalues and wights"""
       
    data = np.loadtxt(file, skiprows=1)
       
    energies = data[:, 3]
    occupations = data[:, 4]
    weights = data[:, 5]
      
       
       #Find CBM, VBM
    CBM = min([x for x in energies if x > fermi])
    VBM = max([x for x in energies if x < fermi])

    sd = weights * occupations (energies - CBM)
    sd = np.sum(sd[energies > CBM]) / np.sum(weights)

    sa = weights * (1 - occupations) * (energies - VBM)
    sa = np.sum(sa[energies < VBM]) / np.sum(weights)

    print("Shallow donor: {}\nShallow acceptor: {}".format(sd, sa))
        
    
