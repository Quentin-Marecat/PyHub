from pyhub.core.basis import Basis 
from time import perf_counter as pc
import numpy as np
import os 

nb_sites = nb_elec = 2
sz = 0.

MBBasis_h = Basis(nb_sites,hilbert=(nb_elec,sz))
print(f'index : {MBBasis_h.index}\nsize : {np.around(MBBasis_h.file_size/10.**9,3)}Go')
print('Spin ordering : ')
MBBasis_h.order = "spin"
print(f'   Many-Body basis\n{MBBasis_h}')
print(f'   full basis\n{MBBasis_h[:]}')
print('\nSite ordering : ')
MBBasis_h.order = "site"
print(f'   Many-Body basis\n{MBBasis_h}')
print(f'   full basis\n{MBBasis_h[:]}')

print(f'\nFock space : ')
MBBasis_f = Basis(nb_sites)
print(f'index : {MBBasis_f.index}\nsize : {np.around(MBBasis_f.file_size/10.**9,3)}Go')
print('Spin ordering : ')
MBBasis_f.order = "spin"
print(f'   Many-Body basis\n{MBBasis_f}')
print('\nSite ordering : ')
MBBasis_f.order = "site"
print(f'   Many-Body basis\n{MBBasis_f}')
print(f'\n   Number of states / 4**N \n{MBBasis_f.nb_states} / {4**nb_sites}')
print(f'\n   Restricted to Hilbert space in site basis :')
hilbert_index = MBBasis_f.hilbert_restricted((nb_elec,sz))
print(MBBasis_f[hilbert_index])


t0=pc()
print(f'\nGenerate large basis : ')
MBBasis_l = Basis(30)
print(f'index : {MBBasis_l.index}\nsize : {np.around(MBBasis_l.file_size/10.**9,3)}Go')
print(f'time {pc()-t0} s')
t0=pc()
print(f'Frome existing file : ')
MBBasis2 = Basis(30)
print(f'index : {MBBasis2.index}\nsize : {np.around(MBBasis2.file_size/10.**9,3)}Go')
print(f'time {pc()-t0} s')
# del MBBasis2
# print(f'size : {np.around(os.path.getsize("basis.h5")/10.**9,3)}Go')