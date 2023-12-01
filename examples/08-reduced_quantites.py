import numpy as np
from pyhub.tools.operators import c_dagger_c, sisj
from pyhub.tools.models import fermi_hubbard
from pyhub.solver.fermi_hubbard import FermiHubbard
from pyhub.core.basis import Basis
import os

nb_sites = 8
nup = ndown = nb_sites//2
hilbert = (nup,ndown)
t_matrix = np.diag(np.full(nb_sites-1,-1.),k=1) + np.diag(np.full(nb_sites-1,-1.),k=-1)
t_matrix[0,-1] = t_matrix[-1,0] = -1.
U = 4.

print('Call FermiHubbard class\nStart Lanczos diagonalisation')
FH = FermiHubbard(nb_sites,*hilbert,t_matrix,U)
FH.kernel(nb_comp_states=2,verbose=False)
print(f'Ground-state density matrix')
FH.reduced_quantities()
print(f'One-RDM up')
print(np.around(FH.one_rdm['up'],3))
print(f'Spin correlation function')
print(FH.spin_cor_matrix)
print(f'S2 : {np.sum(FH.spin_cor_matrix)}')

print('\nCall Operator class\nStart Lanczos diagonalisation')
mbbasis = Basis(nb_sites,hilbert)
H = fermi_hubbard(t_matrix,U)
H.set_basis(mbbasis)
ek,Vk=H.lanczos()

def cicj_up(i,j,psi):
    op = c_dagger_c((i,'up'),(j,'up'))
    op.set_basis(mbbasis)
    return np.around(op.avg(psi),3)

def spin_cor(i,j,psi):
    op = sisj((i,),(j,))
    op.set_basis(mbbasis)
    return np.around(op.avg(psi),3)

print(f'One-RDM up')
print(np.array([[cicj_up(i,j,Vk[:,0]) for i in range(nb_sites)] for j in range(nb_sites)]))
spin_cor_matrix = np.array([[spin_cor(i,j,Vk[:,0]) for i in range(nb_sites)] for j in range(nb_sites)])
print(f'Spin correlation function')
print(spin_cor_matrix)
print(f'S2 : {np.sum(spin_cor_matrix)}')

print('\nFirst excited state using FermiHubbard class')
FH.set_state(1)
FH.reduced_quantities()
print(f'One-RDM up')
print(np.around(FH.one_rdm['up'],3))
print(f'Spin correlation function')
print(FH.spin_cor_matrix)
print(f'S2 : {np.sum(FH.spin_cor_matrix)}')

os.remove('basis.h5')
os.remove('operators.h5')
os.remove('solver.h5')
