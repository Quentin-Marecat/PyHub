import numpy as np 
from pyhub.tools.models import heisenberg
from pyhub.core.basis import Basis
from pyhub.tools.operators import n,opesum,sisj
from itertools import product
from time import perf_counter as pc

nb_sites = nb_elec = 2
J = 1.
J_matrix = np.diag(np.full(nb_sites-1,J),k=1) + np.diag(np.full(nb_sites-1,J),k=-1)
J_matrix[0,-1] = J_matrix[-1,0] = J

print('Set Heisenberg Hamiltonian as an Operator')
H = heisenberg(J_matrix)
print(H) 
print('Set Basis in Hilbert space')
mbbasis = Basis(nb_sites,hilbert=(nb_elec,0.))
H.set_basis(mbbasis)

if H.nstates<200:
    print('Set as matrix')
    t0 = pc()
    Hmatrix = H.to_matrix
    print(Hmatrix)
    print('Start matrix diagonalisation in hilbert space')
    ek,Vk = np.linalg.eigh(Hmatrix)
    print(f'End\ntime {np.around(pc()-t0,4)}\nenergy : {ek[0]}\n')

if H.nstates>200:
    print('Start Lanczos diagonalisation in hilbert space')
    t0 = pc()
    ek,Vk = H.lanczos(maxstep=200,acc_lcz=1.e-8)
    print(f'End\ntime {np.around(pc()-t0,4)}\nenergy : {ek[0]}\n')

print(Vk[:,0])
S2 = opesum([sisj((i,),(j,)) for i,j in product(range(nb_sites),range(nb_sites))])
S2.set_basis(mbbasis)
print(f'S2 {S2.avg(Vk[:,0])}')

print('Find heisenberg subspace')
t0 = pc()
N2 = opesum([n((i,'u'))*n((i,'d')) for i in range(nb_sites)])
N2.set_basis(mbbasis)
heisen_index = np.where(N2 == 0)[0]
print(f'End\ntime {np.around(pc()-t0,4)}')
print(f'Size of the Hilbert Vs Heisenberg basis : {mbbasis.nstates}/{len(heisen_index)}')
print('Set Basis in Heisenberg space')
H.set_basis(mbbasis,selected=heisen_index)

if len(heisen_index)<200:
    print('Set as matrix')
    t0 = pc()
    H.set_basis(mbbasis,selected=heisen_index)
    Hmatrix = H.to_matrix
    print('Start matrix diagonalisation in heisenberg space')
    ek,Vk = np.linalg.eigh(Hmatrix)
    print(f'End\ntime {np.around(pc()-t0,4)}\nenergy : {ek[0]}\n')

if len(heisen_index)>200:
    print('Start Lanczos diagonalisation in hilbert space')
    t0 = pc()
    ek,Vk = H.lanczos(maxstep=200,acc_lcz=1.e-8)
    print(f'End\ntime {np.around(pc()-t0,4)}\nenergy : {ek[0]}\n')
