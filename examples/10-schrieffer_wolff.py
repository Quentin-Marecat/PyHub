import numpy as np 
from pyhub.tools.models import heisenberg
from pyhub.tools.models import fermi_hubbard
from pyhub.solver.fermi_hubbard import FermiHubbard
from pyhub.tools.tools import fidelity
from pyhub.core.basis import Basis
from pyhub.tools.operators import empty_operator, n,c_dagger_c,_n, opeexp
from time import perf_counter as pc

nb_sites = nb_elec = 2
nup = ndown = nb_sites//2
hilbert = (nup,ndown)
t_matrix = np.diag(np.full(nb_sites-1,-1.),k=1) + np.diag(np.full(nb_sites-1,-1.),k=-1)
t_matrix[0,-1] = t_matrix[-1,0] = -1.
U = 10.
U_list = np.full(nb_sites,U)
#theta = 1
theta = (U/4)*np.arctan(4/U)

FH = FermiHubbard(nb_sites,*hilbert,t_matrix,U)
print('Call FermiHubbard class\nStart Lanczos diagonalisation')
t0 = pc()
FH.kernel(max_lcz=200,acc_lcz=1.e-8,store_H = True)
print(f'End\ntime {np.around(pc()-t0,4)}\nenergy : {FH.e0}\n')


print('Set Heisenberg Hamiltonian as an Operator')
t0 = pc()
Hspin = heisenberg(-t_matrix)
print('Set Basis in Hilbert space')
mbbasis = Basis(nb_sites,hilbert=hilbert)
Hspin.set_basis(mbbasis)

if Hspin.nstates<200:
    print('Set as matrix')
    Hmatrix = Hspin.to_matrix
    print('Start matrix diagonalisation in hilbert space')
    ek,Vk = np.linalg.eigh(Hmatrix)
    psi_spin = Vk[:,0]

if Hspin.nstates>200:
    print('Start Lanczos diagonalisation in hilbert space')
    t0 = pc()
    ek,Vk = Hspin.lanczos(maxstep=200,acc_lcz=1.e-8)
    psi_spin = Vk[:,0]


print(f'Set SW Operator')
S = empty_operator()
for p in range(nb_sites):
    for q in range(p):
        lbd_sw = np.array([t_matrix[p,q]/(t_matrix[q,q]-t_matrix[p,p]) if np.abs(t_matrix[p,p]-t_matrix[q,q]) else 0., \
            t_matrix[p,q]/U_list[q], -t_matrix[p,q]/U_list[p], \
                t_matrix[p,q]/((t_matrix[q,q]-t_matrix[p,p]) + (U_list[p]-U_list[q])) if np.abs(U_list[p]-U_list[q]) else 0.])
        for spin in ['up','down']:
            spin_bar ='down' if spin=='up' else 'up'
            if abs(lbd_sw[1])>1.e-10 or abs(lbd_sw[2])>1.e-10:
                S += ( lbd_sw[1] * n((p,spin_bar))* _n((q,spin_bar)) + \
                        lbd_sw[2] * n((q,spin_bar))* _n((p,spin_bar)) \
                        )\
                    * (c_dagger_c((p,spin),(q,spin)) - c_dagger_c((q,spin),(p,spin)))
            if np.abs(lbd_sw[0]) > 1.e-10:
                S += ( lbd_sw[0] * _n((p,spin_bar)) * _n((q,spin_bar)) \
                        )\
                    * (c_dagger_c((p,spin),(q,spin)) - c_dagger_c((q,spin),(p,spin)))
            if np.abs(lbd_sw[3]) > 1.e-10:
                S += ( lbd_sw[3] * n((p,spin_bar)) * n((q,spin_bar)) \
                        )\
                    * (c_dagger_c((p,spin),(q,spin)) - c_dagger_c((q,spin),(p,spin)))

S*=-theta
print('Print S')
for op in S:
    print(op)
S.set_basis(mbbasis)
print(f'Set Fermi-Hubbard Hamiltonian as Operator')
H = fermi_hubbard(t_matrix,U)
H.set_basis(mbbasis)
print(f'Compute Psi fermi-hubbard')
# print(psi_spin )
psi_fh = opeexp(S,psi_spin,unitary = True)
# print(S.to_matrix)
# print(psi_fh)
# print(FH.psi0)
print(f'End\ntime {np.around(pc()-t0,4)}\nenergy error : {100-100*H.avg(psi_fh)/FH.e0} %\nfidelity : {fidelity(FH.psi0,psi_fh)}')