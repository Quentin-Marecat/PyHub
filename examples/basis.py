from pyhub.core.basis import Basis 
from time import perf_counter as pc
import numpy as np
import os 

nb_sites = 2
nup = ndown = nb_sites//2

MBBasis = Basis(nb_sites,hilbert=(nup , ndown))
print(f'size : {np.around(MBBasis.file_size/10.**9,3)}Go')
print('Spin ordering : ')
MBBasis.order = "spin"
print(f'   Many-Body basis\n{MBBasis}')
print(f'   full basis :\n{MBBasis[:]}')
print('\nSite ordering : ')
MBBasis.order = "site"
print(f'   Many-Body basis\n{MBBasis}')
print(f'   full basis\n{MBBasis[:]}')
print(f'{MBBasis.bit_repr()}')

print(f'\nFock space : ')
MBBasis = Basis(nb_sites)
print(f'size : {np.around(MBBasis.file_size/10.**9,3)}Go')
print('Spin ordering : ')
MBBasis.order = "spin"
print(f'   Many-Body basis\n   {MBBasis}')
print(f'{MBBasis.bit_repr()}')
print(f'\n   Number of states / 4**N \n{MBBasis.nb_states} / {4**nb_sites}')
print(f'\n   Restricted to Hilbert space in site basis :')
hilbert_index = MBBasis.hilbert_restricted((nup , ndown))
print(MBBasis[hilbert_index])


t0=pc()
print(f'\nGenerate large basis : ')
MBBasis = Basis(30)
print(f'size : {np.around(MBBasis.file_size/10.**9,3)}Go')
print(f'time {pc()-t0} s')
t0=pc()
print(f'From existing file : ')
MBBasis = Basis(30)
print(f'size : {np.around(MBBasis.file_size/10.**9,3)}Go')
print(f'time {pc()-t0} s')

os.remove('basis.h5')