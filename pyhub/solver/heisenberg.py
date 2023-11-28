#!/usr/bin/env python3
import numpy as np
import h5py
from pyhub.core.core_solver import CoreSolver
import os

__author__ = "Quentin Marécat"
__maintainer__ = "Quentin Marécat"
__email__ = "quentin.marecat@gmail.com"
__date__ = "February, 2023"


class Heisenberg(CoreSolver):
    r"""
    This Class allows solutions of a sz=0 subspace Heisenberg model
    Using exact diagonalization or lanczos solver
    Fortran computation is used for efficiency
    HDF5 file is generated
    """
    def __init__(self,nb_sites,J_matrix:np.ndarray,T=0.): 
        r'''
        H = \sum_{ij\sigma}t_{ij}c^\dagger_{i\sigma}c_{j\sigma} + \sum_{ijkl} U_{ijkl}c^\dagger_{i\sigma}c_{j\sigma}c^\dagger_{k\bar{\sigma}}c_{l\bar{\sigma}}
        '''
        CoreSolver.__init__(self,nb_sites,nb_sites//2,0,'heisenberg.x',T)
        self.J_matrix=J_matrix 
        if not np.isclose(self.J_matrix,self.J_matrix.T).all():
            raise ValueError(f'J_matrix must be symmetric')


    @property
    def filename(self):
        return 'solver.h5'

    def __fileremove__(self):
        try:
            os.remove(self.filename)
        except FileNotFoundError:
            pass

    @property
    def order(self):
        return 'spin'


    def bit_repr(self,index=None):
        if not isinstance(index,np.ndarray):
            index = np.array(range(self.nstates),dtype=int)
        string ='-'*(3*(self.nb_sites+4))+'\n'
        string +='empty        : 0\nuparrow      : 1\ndownarrow    :-1\nupdoawnarrow : 2\n'
        string +='-'*(3*(self.nb_sites+4))+'\n'
        string +='index | bit\n'
        string +='------|--'+'-'*(3*(self.nb_sites+1))+'\n'
        for index in self.basis[index]:
            up = index & 2**self.nb_sites-1
            down = 2**self.nb_sites-1-up
            string += "{:4d}  |".format(index)
            for i in range(self.nb_sites):
                if ((up >> i) & 1)!=1 and ((down >> i) & 1)!=1:
                    string += ' 0 '
                elif ((up >> i) & 1)==1 and ((down >> i) & 1)!=1:
                    string += ' 1 '
                elif ((up >> i) & 1)!=1 and ((down >> i) & 1)==1:
                    string += '-1 '
                else:
                    string += ' 2 '
            string+='\n'
        return string