from corrections import *
from sys import argv

file = argv[1]
with open(file, 'r') as struct:
    
    struct.readline()
    
    a1 = list(map(float, struct.readline().split())) 
    a2 = list(map(float, struct.readline().split()))
    a3 = list(map(float, struct.readline().split()))
             
n = int(argv[2])

grid =  uniformK(a1, a2, a3, n)

for i in range(n):
    for j in range(n):
        for k in range(n):

            print(str(grid[i, j, k, :])[1:-1])
