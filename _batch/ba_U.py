import numpy as np
from bethe_anzatz import E_bessel

dx=1.e-6
u_range = np.linspace(0.005,0.995,2000)

for U in 4*u_range/(1-u_range):
    e0,ni = E_bessel(U),(E_bessel(U+dx)-E_bessel(U-dx))/(2*dx)
    with open('dat/BA.dat','a') as f:
        f.write(f'{U} {e0} {e0-U*ni} {ni}\n')