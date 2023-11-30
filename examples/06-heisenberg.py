import numpy as np
from pyhub.tools.operators import n,_n,opesum, splus,sminus,sz, sisj, empty_operator
from pyhub.tools.operators import c_dagger_c, c_dagger,_c
from itertools import product
from pyhub.core.basis import Basis
from time import perf_counter as pc
from pyhub.tools.tools import fidelity
from pyhub.solver.heisenberg import Heisenberg

np.set_printoptions(precision=4)
nb_sites = 2
nup = ndown = nb_sites//2
hilbert = (nup,ndown)

np.random.seed(10)
J_matrix = np.random.random((nb_sites,nb_sites))-0.5
J_matrix = abs(J_matrix + J_matrix.T - 2*np.diag(np.diag(J_matrix)))
J_matrix += np.diag(np.full(nb_sites,-1.e0))

# Hs = Heisenberg(nb_sites,J_matrix)
# Hs.kernel(store_H=True)
# print(Hs.e0)

def spin_ham1(J_matrix):
    return opesum([J_matrix[i,j]*sisj((i,),(j,)) for i,j in product(range(nb_sites),range(nb_sites))])
    
def spin_ham2(J_matrix):
    H = empty_operator()
    for p in range(nb_sites):
        for q in range(nb_sites):
            if p==q :
                H += J_matrix[p,q]* ( n((p,'up'))*_n((p,'down')) + n((p,'down'))*_n((p,'up')) )/2
            else:
                H += J_matrix[p, q] * (\
                    -0.5 * ( c_dagger_c((p,'u'),(q,'u')) * c_dagger_c((q,'d'),(p,'d')) + c_dagger_c((p,'d'),(q,'d')) * c_dagger_c((q,'u'),(p,'u')) ) \
                )
            H += J_matrix[p, q] * ( ((n((p,'u')) - n((p,'d'))) / 2.) * ((n((q,'u')) - n((q,'d'))) / 2.) )
    return H

def spin_ham3(J_matrix):
    return opesum([J_matrix[i,j]*(splus((i,))*sminus((j,))/2 + sminus((j,))*splus((i,))/2 + sz((i,))*sz((j,)) ) for i,j in product(range(nb_sites),range(nb_sites))])

def spin_ham4(J_matrix):
    H = empty_operator()
    for p in range(nb_sites):
        for q in range(nb_sites):
            if p==q :
                H += J_matrix[p,q]* ( c_dagger((p,'u'))*_c((p,'u'))*(1- c_dagger((p,'d'))*_c((p,'d'))) + c_dagger((p,'d'))*_c((p,'d'))*(1-c_dagger((p,'u'))*_c((p,'u'))) )/2
            else:
                H += J_matrix[p, q] * (\
                    -0.5 * ( c_dagger((p,'u'))*_c((q,'u')) * c_dagger((q,'d'))*_c((p,'d')) + c_dagger((p,'d'))*_c((q,'d')) * c_dagger((q,'u'))*_c((p,'u')) ) \
                )
            H += J_matrix[p, q] * ( ((c_dagger((p,'u'))*_c((p,'u')) - c_dagger((p,'d'))*_c((p,'d'))) / 2.) * ((c_dagger((q,'u'))*_c((q,'u')) - c_dagger((q,'d'))*_c((q,'d'))) / 2.) )
    return H

print('Def many-body basis in hilbert space')
mbbasis = Basis(nb_sites,hilbert)

print('Definie spin-hamiltonian in the 1st way:')
H1 = spin_ham1(J_matrix)
print(H1)
print('Solve the Hamiltonian')
H1.set_basis(mbbasis)
ek,Vk = np.linalg.eigh(H1.to_matrix)
#print(H1.to_matrix)
print(f'Ground-state energy {ek[0]}\n')

print('Definie spin-hamiltonian in the 2nd way:')
H2 = spin_ham2(J_matrix)
print(H2)
print('Solve the Hamiltonian')
H2.set_basis(mbbasis)
ek,Vk = np.linalg.eigh(H2.to_matrix)
#print(H2.to_matrix)
print(f'Ground-state energy {ek[0]}\n')

print('Definie spin-hamiltonian in the 3rd way:')
H3 = spin_ham3(J_matrix)
print(H3)
print('Solve the Hamiltonian')
H3.set_basis(mbbasis)
ek,Vk = np.linalg.eigh(H3.to_matrix)
#print(H3.to_matrix)
print(f'Ground-state energy {ek[0]}\n')
print(f'Energy from eigenvectors : {H3.avg(Vk[:,0])}')
print('/!\ Warning : S+ and S- not stable in hilbert space\nSet H3.stable as False -> Much slower calculation')
H3.stable = False
print(f'Energy from eigenvectors : {H3.avg(Vk[:,0])}')

print('Definie spin-hamiltonian in the 4th way:')
H4 = spin_ham4(J_matrix)
print(H4)
print('Solve the Hamiltonian')
H4.set_basis(mbbasis)
ek,Vk = np.linalg.eigh(H4.to_matrix)
#print(H4.to_matrix)
print(f'Ground-state energy {ek[0]}\n')
print(f'Energy from eigenvectors : {H4.avg(Vk[:,0])}')
print('/!\ Warning : c_dagger and c not stable in hilbert space\nSet H4.stable as False -> Much slower calculation')
H4.stable = False
print(f'Energy from eigenvectors : {H4.avg(Vk[:,0])}')