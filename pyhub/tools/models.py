import numpy as np 
from pyhub.tools.operators import c_dagger_c, n, sisj,  opesum
from itertools import product

def fermi_hubbard(t_matrix,U):
    nb_sites=t_matrix.shape[0]
    return opesum([t_matrix[i,j]*c_dagger_c((i,spin),(j,spin)) for i,j,spin in product(range(nb_sites),range(nb_sites),['up','down'])])\
        + U*opesum([n((i,'up'))*n((i,'down')) for i in range(nb_sites)])

def heisenberg(J_matrix):
    nb_sites=J_matrix.shape[0]
    return opesum([J_matrix[i,j]*sisj((i,),(j,)) for i,j in product(range(nb_sites),range(nb_sites))])
    # for p in range(nb_sites):
    #     for q in range(nb_sites):
    #         if np.abs(J_matrix[p, q])>1.e-10:
    #             H += J_matrix[p, q] * (\
    #                 -0.5 * ( c_dagger_c((p,'u'),(q,'u')) * c_dagger_c((q,'d'),(p,'d')) + c_dagger_c((p,'d'),(q,'d')) * c_dagger_c((q,'u'),(p,'u')) )  + \
    #                 ((n((p,'u')) - n((p,'d'))) / 2.) * ((n((q,'u')) - n((q,'d'))) / 2.) \
    #             )
    # return H