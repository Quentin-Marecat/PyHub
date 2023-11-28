from pyhub.solver.heisenberg import Heisenberg
import numpy as np

np.set_printoptions(precision=4)
nb_sites=10
T = 0.


J_matrix = np.diag(np.full(nb_sites-1,1.),k=1) + np.diag(np.full(nb_sites-1,1.),k=-1)
#    t_matrix += np.diag(np.full(nb_sites,-0.1),k=0)
#    t_matrix += np.diag(np.full(nb_sites-2,-0.1),k=2) + np.diag(np.full(nb_sites-2,-.1),k=-2)
J_matrix[0,-1],J_matrix[-1,0] = 1.,1.  

J_matrix *= 1. ##anti-ferro

H = Heisenberg(nb_sites,J_matrix,T=T)
H.kernel(max_lcz=390,acc_lcz = 1.e-8,nb_comp_states=1,\
    compute_rq=True,verbose=True)
print('Ground states')
print(f'e0(T)   {H.e0} \ne_0(T=0) {H.e[0]} ')
#print(f'boltzmann {H.boltzmann} / {H.Z}')
print(f'SiSj\n{H.spin_cor_matrix}')
print(f'S^2 {H.s2}')
print(f'number of state elements : {H.nb_states}')
if nb_sites<=6: 
    print(H.bit_repr())


if (H.nb_comp_states>1):
    H.set_state(1)
    print('First excited states')
    print(f'e1 {H.e[1]}')
    H.reduced_quantities()
    print(f'SiSj\n{H.spin_cor_matrix}')
    print(f'S^2 {H.s2}')

H.__fileremove__()
