#!/usr/bin/env python3
import numpy as np
import h5py
from pyhub.core.core_solver import CoreSolver
import os

__author__ = "Quentin Marécat"
__maintainer__ = "Quentin Marécat"
__email__ = "quentin.marecat@gmail.com"
__date__ = "February, 2023"


class FermiHubbard(CoreSolver):
    r"""
    This Class allows solutions of a general hubbard cluster
    Using exact diagonalization or lanczos solver
    Fortran computation is used for efficiency
    HDF5 file is generated
    """
    def __init__(self,nb_sites,n_up:int,n_down:int,t_matrix:np.ndarray,U:float=0,T=0.,order='spin'): 
        r'''
        H = \sum_{ij\sigma}t_{ij}c^\dagger_{i\sigma}c_{j\sigma} + \sum_{ijkl} U_{ijkl}c^\dagger_{i\sigma}c_{j\sigma}c^\dagger_{k\bar{\sigma}}c_{l\bar{\sigma}}
        '''
        CoreSolver.__init__(self,nb_sites,n_up,n_down,'hubbard.x',T,order=order)
        self.t_matrix=t_matrix 
        if not np.isclose(self.t_matrix,self.t_matrix.T).all():
            raise ValueError(f't_matrix must be symmetric')
        self.ek,self.Vk = np.linalg.eigh(self.t_matrix)
        self.U=U

    def __fileremove__(self):
        try:
            os.remove(self.filename)
        except FileNotFoundError:
            pass
