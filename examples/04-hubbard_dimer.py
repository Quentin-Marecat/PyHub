import numpy as np 
from pyhub.tools.operators import c_dagger_c,n
from pyhub.core.basis import Basis
import os
from pyhub.tools.tools import fidelity

nb_sites = 2
t_matrix = np.zeros((2,2))
t_matrix[0,1] = t_matrix[1,0] = -1.
U = 4.

print('Set the operator')
H = t_matrix[1,0]*c_dagger_c((1,'up'),(0,'up')) + t_matrix[1,0]*c_dagger_c((1,'down'),(0,'down')) +\
    t_matrix[0,1]*c_dagger_c((0,'up'),(1,'up')) + t_matrix[0,1]*c_dagger_c((0,'down'),(1,'down')) +\
    U*n((0,('up')))*n((0,('down'))) + U*n((1,('up')))*n((1,('down')))


print(H)
print('\nSet many-body basis in Fock space')
mbbasis_fock = Basis(nb_sites)
print('To use the operator, set the Basis set using Basis class\n')
H.set_basis(mbbasis_fock)

print('Solve the Operator in Basis space')
ek,Vk = np.linalg.eigh(H.to_matrix)
print(f'Ground-state energy {ek[0]}')
N = n((0,('up'))) + n((0,('down'))) + n((1,('up'))) + n((1,('down')))
N.set_basis(mbbasis_fock)
print(f'Number of particles in the ground state : {np.around(N.avg(Vk[:,0]),3)}')

print('\nWe want to solve the Operator in the 2-particle subspace')
index = np.where(N==2)[0]
ek,Vk = np.linalg.eigh(H[index].to_matrix)
print(f'Ground-state energy in the 2-particle subspace {ek[0]}')
print(f'Number of particles in the ground state : {np.around(N[index].avg(Vk[:,0]),3)}')
print('\nUse the correct chemical potential')

mu = U/2
H_mu = H - mu*N
H_mu.set_basis(mbbasis_fock)
ek,Vk = np.linalg.eigh(H_mu.to_matrix)
print(f'Ground-state energy {ek[0] + 2*mu}')
print(f'Number of particles in the ground state : {np.around(N.avg(Vk[:,0]),3)}')

os.remove('*.h5')
