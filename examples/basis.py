from pyhub.core.basis import Basis 
from time import perf_counter as pc

nb_sites = nb_elec = 4
sz = 0.

MBBasis = Basis(nb_sites,nb_elec,sz)
print('Spin ordering : ')
MBBasis.order = "spin"
print(f'   Many-Body basis\n{MBBasis}')
print(f'   full basis\n{MBBasis[:]}')
print('\nSite ordering : ')
MBBasis.order = "site"
print(f'   Many-Body basis\n{MBBasis}')
print(f'   full basis\n{MBBasis[:]}')

print(f'\nFock space : ')
MBBasis = Basis(nb_sites,nb_elec,sz,fock=True)
print('Spin ordering : ')
MBBasis.order = "spin"
print(f'   Many-Body basis\n{MBBasis}')
print('\nSite ordering : ')
MBBasis.order = "site"
print(f'   Many-Body basis\n{MBBasis}')
print(f'\n   Number of states / 4**N \n{MBBasis.nb_states} / {4**nb_sites}')
print(f'\n   Restricted to Hilbert space in site basis :')
MBBasis.hilbert_restricted(nb_elec,sz)
print(MBBasis[MBBasis.hilbert_index])


t0=pc()
print(f'\nGenerate large basis : ')
MBBasis = Basis(30,None,None)
print(f'time {pc()-t0} s')
t0=pc()
print(f'Frome existing file : ')
MBBasis = Basis(30,None,None)
print(f'time {pc()-t0} s')
MBBasis.__basisremove__()