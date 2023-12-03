from pyhub.core.basis import Basis 
from time import perf_counter as pc
import numpy as np
import os 

nb_sites = 2  ## number of sites to initialize basis class
print(f'Set Basis in Fock space : ')
MBBasis = Basis(nb_sites)
print(f'full basis :\n{MBBasis[:]}')  ## Basis class can be used as an np.array
print(f'{MBBasis.bit_repr()}')
# print(f'Restriction to (nup=1,ndown=1) Hilbert space:')
# hilbert_index = MBBasis.hilbert_restricted((nb_sites//2 , nb_sites//2))
# print(MBBasis[hilbert_index])  
print('\n'+"-"*40)

print(f'Set basis directly in Hilbert subspace')
nup = ndown = nb_sites//2  ## Restrict the basis in the Hilbert subspace
MBBasis = Basis(nb_sites,hilbert=(nup , ndown))
print('\nSpin ordering : ')
MBBasis.order = "spin" ## Set the basis in spin ordering, default : True
print(f'{MBBasis}')
print('\nbit string representation in the given ordering, here "spin" : ')
print(f'{MBBasis.bit_repr()}')
print('\nSite ordering : ')
MBBasis.order = "site" ## Set the basis in site ordering, default : False
print(f'{MBBasis}')
print('\n'+"-"*40)

t0=pc()
print(f'\nGenerate large basis : ')
MBBasis = Basis(20)
print(f'size : {np.around(MBBasis.file_size/10.**9,3)}Go')
print(f'time {pc()-t0} s')
t0=pc()
print(f'From existing file : ')
MBBasis = Basis(20)
print(f'size : {np.around(MBBasis.file_size/10.**9,3)}Go')
print(f'time {pc()-t0} s')

os.remove('basis.h5')