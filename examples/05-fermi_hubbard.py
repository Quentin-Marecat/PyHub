import numpy as np 
from pyhub.solver.fermi_hubbard import FermiHubbard
from pyhub.tools.models import fermi_hubbard
from pyhub.core.basis import Basis
from time import perf_counter as pc
from pyhub.tools.tools import fidelity

nb_sites = 12
nup = ndown = nb_sites//2
hilbert = (nup,ndown)
t_matrix = np.diag(np.full(nb_sites-1,-1.),k=1) + np.diag(np.full(nb_sites-1,-1.),k=-1)
t_matrix[0,-1] = t_matrix[-1,0] = -1.
U = 4.

FH = FermiHubbard(nb_sites,*hilbert,t_matrix,U)
print('Call FermiHubbard class\nStart Lanczos diagonalisation')
t0 = pc()
FH.kernel(max_lcz=100,acc_lcz=1.e-8,store_H = False,nb_comp_states=2)
print(f'time {np.around(pc()-t0,4)}\nenergy : {FH.e0}\n')
print(f'Set the first excited state')
FH.set_state(1)
print(f'energy : {FH.e0}\n')

print('Set Fermi Hubbard Hamiltonian as an Operator')
H = fermi_hubbard(t_matrix,U)
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
print(f'time {np.around(pc()-t0,4)}\nenergy : {ek[0]}\nfidelity : {fidelity(FH.psi0,Vk[:,0])}')