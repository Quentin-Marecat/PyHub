from pyhub.solver.heisenberg import Heisenberg
import numpy as np
from pyhub.tools.models import heisenberg
from pyhub.tools.operators import n,opesum
from pyhub.core.basis import Basis
from time import perf_counter as pc
from pyhub.tools.tools import fidelity
import os

np.set_printoptions(precision=4)
nb_sites = 6
nup = ndown = nb_sites//2
hilbert = (nup,ndown)


J_matrix = np.diag(np.full(nb_sites-1,1.),k=1) + np.diag(np.full(nb_sites-1,1.),k=-1)
J_matrix[0,-1],J_matrix[-1,0] = 1.,1.  

J_matrix *= 1. ##anti-ferro

t0 = pc()
print('Call Heisenberg class\nStart Lanczos diagonalisation')
Hs = Heisenberg(nb_sites,J_matrix,T=0.)
Hs.kernel(max_lcz=100,acc_lcz = 1.e-8,nb_comp_states=1,\
    compute_rq=True,verbose=True)
print(f'End\ntime {np.around(pc()-t0,4)}\nenergy : {Hs.e0}\n')

print('Set Heisenberg Hamiltonian as an Operator')
H = heisenberg(J_matrix)
print('Set Basis in Hilbert space')
mbbasis = Basis(nb_sites,hilbert=hilbert)
H.set_basis(mbbasis)

if H.nstates<200:
    print('Set as matrix')
    t0 = pc()
    Hmatrix = H.to_matrix
    print('Start matrix diagonalisation')
#    ek,Vk = eigsh(Hmatrix,k=1,which='SA')
    ek,Vk = np.linalg.eigh(Hmatrix)
    print(f'End\ntime {np.around(pc()-t0,4)}\nenergy : {ek[0]}\n')

print('Start Lanczos diagonalisation')
t0 = pc()
ek,Vk = H.lanczos(maxstep=200,acc_lcz=1.e-8)
print(f'End\ntime {np.around(pc()-t0,4)}\nenergy : {ek[0]}')

N2 = opesum([n((i,'u'))*n((i,'d')) for i in range(nb_sites)])
N2.set_basis(mbbasis)
print('\nFind elements in the Hilbert space that belong to a target subspace')
print('Simple occupated subspace')
spin_index = np.where(N2 == 0)[0]
print(f'Size of the Spin basis {len(spin_index)}')
print('\nSet new operator in restricted basis')
print('Set as matrix')
t0 = pc()
print('Start matrix diagonalisation')
ek,Vk = np.linalg.eigh(H[spin_index].to_matrix)
print(f'End\ntime {np.around(pc()-t0,4)}\nenergy : {ek[0]}')
print(f'fidelity with exact gs : {fidelity(Hs.psi0,Vk[:,0])}')


nb_sites=16
print(f'\nSolve Spin Hamiltonian for large number of sites {nb_sites}')
print(f'Possible only for Heisenberg class')
nup = ndown = nb_sites//2
hilbert = (nup,ndown)


J_matrix = np.diag(np.full(nb_sites-1,1.),k=1) + np.diag(np.full(nb_sites-1,1.),k=-1)
J_matrix[0,-1],J_matrix[-1,0] = 1.,1.  

J_matrix *= 1. ##anti-ferro

t0 = pc()
print('Call Heisenberg class\nStart Lanczos diagonalisation')
Hs = Heisenberg(nb_sites,J_matrix,T=0.)
Hs.kernel(max_lcz=300,acc_lcz = 1.e-8,nb_comp_states=1,\
    compute_rq=False,verbose=True)
print(f'End\ntime {np.around(pc()-t0,4)}\nenergy : {Hs.e0}\n')

os.remove('basis.h5')
os.remove('operators.h5')
os.remove('solver.h5')