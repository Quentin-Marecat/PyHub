from pyhub.solver import hubbard
import numpy as np

np.set_printoptions(precision=4)
nb_sites=10
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

U = 0.
J=np.zeros((nb_sites,nb_sites))

#    J = -np.copy(t_matrix)
#    U,t_matrix = 0.,np.zeros((nb_sites,nb_sites))
S = hubbard.Solve(nb_sites,nb_elec,sz,t_matrix,U,J=J,T=T)
S.run(max_lcz=1000,acc_lcz = 1.e-8,nb_comp_states=2,\
    compute_reduced_quantities=True,compute_two_body=False,compute_spgf=False,verbose=False)
        
    print('Ground states')
    print(f'e0(beta) {S.e0} / e_0 {S.e[0]} ')
#    print(f'boltzmann {S.boltzmann} / {S.Z}')
    print(f'one body energy {S.e_one_body}\ntwo body energy {S.e_two_body}\ne {S.e_one_body+S.e_two_body}')
    print(f'density matrix up\n{S.density_matrix["up"]}')
    print(f'density matrix down\n{S.density_matrix["down"]}')
    etak,ek = np.linalg.eigh(S.density_matrix["up"])
    print(etak)
#    print(f'density matrix calc from gf\n{S.density_matrix_from_gf["up"]}')
    print(f'ni {S.ni}')
    print(f'SiSj\n{S.spin_cor_matrix}')
    print(f'S^2 {S.s2}')
    print(f'quasi particle weight {S.qpw[0,0]}')


    if (S.nb_comp_states>1):
        print('First excited states')
        print(f'e1 {S.e[1]}')
        S.reduced_quantities(1)
        print(f'one body energy {S.e_one_body}\ntwo body energy {S.e_two_body}\ne {S.e_one_body+S.e_two_body}')
        print(f'density matrix\n{S.density_matrix["up"]}')
        print(f'ni {S.ni}')
        print(f'SiSj\n{S.spin_cor_matrix}')
        print(f'S^2 {S.s2}')
