#!/usr/bin/env python3
import numpy as np


def d1(w,t=-1.):
    if isinstance(w,(list,np.ndarray)):
        return np.array([d1(w_,t=t) for w_ in w])
    B = 2*abs(t)
    if B**2>w**2:
        return 2./(np.pi*np.sqrt(B**2 - w**2))
    else:
        return 0

def semicircular(w,t=-1):
    if isinstance(w,(list,np.ndarray)):
        return np.array([semicircular(w_,t=t) for w_ in w])
    B = 2*abs(t) ## t renormalized
    if B**2>w**2:
        return 2./(np.pi*B**2)*np.sqrt(B**2 - w**2)
    else:
        return 0
