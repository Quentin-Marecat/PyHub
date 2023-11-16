import numpy as np
from bethe_anzatz import batch
from scipy.interpolate import CubicSpline as CS

delta = 1.e-5
nrange = np.linspace(0.,1.,100)
ba = {key1:{key2:np.zeros(100) for key2 in ['n','etot','ekin','ni']} for key1 in [1,4,8]}
for U in [1,4,8]:
    dct = batch(U)
    e = CS(dct['particles_per_site'],dct['energy'])
    dct = batch(U+delta)
    epdelta = CS(dct['particles_per_site'],dct['energy'])
    dct = batch(U-delta)
    emdelta = CS(dct['particles_per_site'],dct['energy'])
    ba['n'] = nrange
    ba['etot'] = e(nrange)
    ba['ni'] = (epdelta(nrange) - emdelta(nrange))/(2*delta)
    ba['ekin'] = ba['etot'] - U*ba['ni']
    with open(f'dat/BA_U{U}.dat','a') as f:
        for i in range(100):
            f.write(f"{ba['n'][i]} {ba['etot'][i]} {ba['ekin'][i]} {ba['ni'][i]}\n")


