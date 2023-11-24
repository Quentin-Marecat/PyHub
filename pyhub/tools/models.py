import numpy as np 
from pyhub.tools.operators import c_dagger_c, n, sisj, empty_operator

def fermi_hubbard(t_matrix,U):
    nb_sites=t_matrix.shape[0]
    H = empty_operator()
    for i,_t in enumerate(t_matrix):
        for j,t in enumerate(_t):
            if np.abs(t)>1.e-14:
                for spin in ['up','down']:
                    H += t*c_dagger_c((i,spin),(j,spin))
    for i in range(nb_sites):
        H += U*n((i,'up'))*n((i,'down'))
    return H

def heisenberg(J_matrix):
    H = empty_operator()
    for i,_t in enumerate(J_matrix):
        for j,t in enumerate(_t):
            if np.abs(t)>1.e-14:
                H += t*sisj((i,),(j,))
    # for p in range(nb_sites):
    #     for q in range(nb_sites):
    #         if np.abs(J_matrix[p, q])>1.e-10:
    #             H += J_matrix[p, q] * (\
    #                 -0.5 * ( c_dagger_c((p,'u'),(q,'u')) * c_dagger_c((q,'d'),(p,'d')) + c_dagger_c((p,'d'),(q,'d')) * c_dagger_c((q,'u'),(p,'u')) )  + \
    #                 ((n((p,'u')) - n((p,'d'))) / 2.) * ((n((q,'u')) - n((q,'d'))) / 2.) \
    #             )
    return H