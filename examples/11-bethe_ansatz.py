import numpy as np
from pyhub.solver import bethe_anzatz

U=4
d=bethe_anzatz.solve(U)
print(f'E_gs {d["energy"]}')
print(f'Particle per site {d["particles_per_site"]}')
print(f'From Bessel {bethe_anzatz.E_bessel(U)}')
print(f'Batching : used away from half-filling')
dct = bethe_anzatz.batch(U)
print('fill   energy   mag')
for i,val in enumerate(dct['energy']):
    print(dct['particles_per_site'][i],val,dct['magnetization'][i])
