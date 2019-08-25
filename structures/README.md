## File Format:
All structure parameters are specified in bohr.
If the first line is a non-0 number then the lattice vectors are expressed in unit of that number.
If there is no number or the number is 0 then the lattice vectors are expressed in bohr.
The files have the format:

Alat<br/>
v1x v1y v1z<br/>
v2x v2y v2z<br/>
v3x v3y v3z<br/>
At x y z <br/>
At x y z<br/> 
At x y z<br/>
.....<br/>
.....<br/>
.....<br/>
  

## File Name Convention:
SC = simple cube<br/>
r  = rhombohedral<br/>
O  = Orthorhombic<br/>

The number indicates the facter by which each lattice vector has been scaled to create the super cell, e.g. SC333 indicates a 3x3x3 simple cubic super cell.
Vac denotes a vacancy. A_B denotes an substitutional/antisite defect where atom A replaces B and A_i denotes an interstitial.
Forexample r222_GeVacN1 denotes a 2x2x2 rhombohedral super cell with a Ge vacancy and a -1 charge.
