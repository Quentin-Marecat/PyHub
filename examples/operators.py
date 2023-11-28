import numpy as np 
from pyhub.tools.models import fermi_hubbard
from pyhub.core.basis import Basis
from pyhub.solver.fermi_hubbard import FermiHubbard
from time import perf_counter as pc
from scipy.sparse.linalg import eigsh
def fidelity(psi, phi):
    return np.absolute(np.dot(np.conjugate(psi), phi))**2

nb_sites = nb_elec = 10
t_matrix = np.diag(np.full(nb_sites-1,-1.),k=1) + np.diag(np.full(nb_sites-1,-1.),k=-1)
t_matrix[0,-1] = t_matrix[-1,0] = -1.
U = 4.

FH = FermiHubbard(nb_sites,nb_elec,0.,t_matrix,U)
#print(FH.basis)
print('Call FermiHubbard class\nStart Lanczos diagonalisation')
t0 = pc()
FH.kernel(max_lcz=200,acc_lcz=1.e-8,store_H = True)
print(f'End\ntime {np.around(pc()-t0,4)}\nenergy : {FH.e0}\n')

print('Set Fermi Hubbard Hamiltonian as an Operator')
H = fermi_hubbard(t_matrix,U)
print('Set Basis in Hilbert space')
mbbasis = Basis(nb_sites,hilbert=(nb_elec,0.))
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
print(f'End\ntime {np.around(pc()-t0,4)}\nenergy : {ek[0]}\nfidelity : {fidelity(FH.psi0,Vk[:,0])}')