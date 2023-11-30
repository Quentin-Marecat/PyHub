from itertools import product
import numpy as np 
from pyhub.tools.operators import n,c_dagger_c,opesum
from pyhub.core.basis import Basis
from time import perf_counter as pc
import os

nb_sites = 2
print('Define counting operator over all sites')
N=0
for i in range(nb_sites):
    N += n((i,'up'))+n((i,'down'))
print(N)
print('\nSame, but using opesum function')
N = opesum([n((i,spin)) for i,spin in product(range(nb_sites),('up','down'))])
print(N)
print('\n'+"-"*40)

print('\nTo use the operator, set the Basis set using Basis class')
nup = ndown = nb_sites//2
hilbert = (nup,ndown)
mbbasis_hilbert = Basis(nb_sites,hilbert=hilbert)
N.set_basis(mbbasis_hilbert)
print('\nCompute the operator in the basis set')
print(N.to_matrix)
print('\nSet a normalized random vector')
np.random.seed(10)
v = np.random.random(N.nstates) - 0.5
v/=np.linalg.norm(v)
print('\nApply N operator over V : N@v or N(v)')
print(N@v)
print('\nCompute average value of N over v, must be 2.')
print(N.avg(v))
print('\n'+"-"*40)

os.remove('basis.h5')
os.remove('operator.h5')