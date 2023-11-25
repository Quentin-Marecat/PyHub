from pyhub.solver.fermi_hubbard import FermiHubbard
import numpy as np

np.set_printoptions(precision=4)
nb_sites=8
nb_elec=nb_sites
sz=-0.
T = 0.

### generate physical correlation potential
np.random.seed(10)
v = np.random.random(nb_sites)
v/=np.linalg.norm(v)
P = np.identity(nb_sites) - 2*np.einsum('i,j->ij',v,v)
ek,vk=np.linalg.eigh(P)
U_ = np.random.random(nb_sites)
U = np.einsum('m,im,jm,km,lm->ijkl',U_,P,P,P,P)

t_matrix = np.diag(np.full(nb_sites-1,-1.),k=1) + np.diag(np.full(nb_sites-1,-1.),k=-1)
#    t_matrix += np.diag(np.full(nb_sites,-0.1),k=0)
#    t_matrix += np.diag(np.full(nb_sites-2,-0.1),k=2) + np.diag(np.full(nb_sites-2,-.1),k=-2)
t_matrix[0,-1],t_matrix[-1,0] = -1.,-1.  

U = 4.

FH = FermiHubbard(nb_sites,nb_elec,sz,t_matrix,U,T=T,fock=False)
FH.kernel(max_lcz=20000,acc_lcz = 1.e-8,nb_comp_states=1,\
    compute_rq_hubbard=True,compute_two_body=False,compute_spgf=False,verbose=True)
        
print('Ground states')
print(f'e0(beta) {FH.e0} / e_0 {FH.e[0]} ')
#    print(f'boltzmann {FH.boltzmann} / {FH.Z}')
print(f'one body energy {FH.e_one_body}\ntwo body energy {FH.e_two_body}\ne {FH.e_one_body+FH.e_two_body}')
print(f'density matrix up\n{FH.density_matrix["up"]}')
print(f'density matrix down\n{FH.density_matrix["down"]}')
etak,ek = np.linalg.eigh(FH.density_matrix["up"])
#    print(f'density matrix calc from gf\n{FH.density_matrix_from_gf["up"]}')
print(f'ni {FH.ni}')
print(f'SiSj\n{FH.spin_cor_matrix}')
print(f'S^2 {FH.s2}')


if (FH.nb_comp_states>1):
    print('First excited states')
    print(f'e1 {FH.e[1]}')
    FH.reduced_quantities(1)
    print(f'one body energy {FH.e_one_body}\ntwo body energy {FH.e_two_body}\ne {FH.e_one_body+FH.e_two_body}')
    print(f'density matrix\n{FH.density_matrix["up"]+FH.density_matrix["down"]}')
    print(f'ni {FH.ni}')
    print(f'SiSj\n{FH.spin_cor_matrix}')
    print(f'S^2 {FH.s2}')

FH.__basisremove__()
FH.__hubbardremove__()
