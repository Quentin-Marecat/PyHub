from itertools import product
import numpy as np 
from pyhub.tools.operators import n,sisj,opesum
from pyhub.core.basis import Basis
import os

nb_sites = 2
nup = ndown = nb_sites//2
hilbert = (nup,ndown)
mbbasis_hilbert = Basis(nb_sites,hilbert=hilbert)
print('Define double counting operator over all sites')
N2 = opesum([n((i,'u'))*n((i,'d')) for i in range(nb_sites)])
N2.set_basis(mbbasis_hilbert)
print(N2.bit_repr())
print('\nFind elements in the Hilbert space that belong to a target subspace')
print('\nSimple occupated subspace')
spin_index = np.where(N2 == 0)[0]
print('\nSet new operator in restricted basis')
N2_restricted = N2[spin_index]
## equivalent to : 
# N2_restricted = opesum([n((i,'u'))*n((i,'d')) for i in range(nb_sites)])
# N2_restricted.set_basis(mbbasis_hilbert,selected=spin_index)
print('\nNumber of doubly occupated states')
print(N2[spin_index].bit_repr())
print('\nNow we want to find stable S2 space')
S2 = opesum([sisj((i,),(j,)) for i,j in product(range(nb_sites),range(nb_sites))])
S2.set_basis(mbbasis_hilbert)

# H = np.zeros((S2.nb_states,S2.nb_states))
# for i in range(S2.nb_states):
#     psi = np.zeros(S2.nb_states)
#     psi[i]=1.
#     H[i,:] = S2@psi 
# ek,Vk=np.linalg.eigh(H)
# print(ek)
index_ = S2.find_stable_space()
for i,index in enumerate(index_):
    print(f'\nStable space n°{i}')
    print(S2[index].bit_repr())
    print(f'Eigenvalues')
    ek,Vk = np.linalg.eigh(S2[index].to_matrix) 
    print(ek)

os.remove('basis.h5')
os.remove('operators.h5')