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
    def __init__(self,nb_sites,nb_elec:int,sz:float,t_matrix:np.ndarray,U:float=0,T=0.,fock=False,order='spin'): 
        r'''
        H = \sum_{ij\sigma}t_{ij}c^\dagger_{i\sigma}c_{j\sigma} + \sum_{ijkl} U_{ijkl}c^\dagger_{i\sigma}c_{j\sigma}c^\dagger_{k\bar{\sigma}}c_{l\bar{\sigma}}
        '''
        CoreSolver.__init__(self,nb_sites,nb_elec,sz,T,fock=fock,order=order)
        self.t_matrix=t_matrix 
        self.ek,self.Vk = np.linalg.eigh(self.t_matrix)
        self.U=U
        self.is_Uloc=True if len(self.U_matrix.shape)==1 else False

    @property
    def filename_hubbard(self):
        return 'hubbard.h5'

    def __hubbardremove__(self):
        try:
            os.remove(self.filename_hubbard)
        except FileNotFoundError:
            pass

    def write_in_hubbard(self):
        file=h5py.File(self.filename_hubbard,'a')
        grp=file.create_group('input')
        grp.attrs['nb_sites'] = self.nb_sites
        grp.attrs['nb_elec'] = self.nb_elec
        grp.attrs['sz'] = self.sz
        grp.attrs['do_solution']=True
        grp.attrs['do_rq']=self.compute_rq_hubbard
        grp.attrs['do_rq_two_body']=self.compute_two_body
        grp.attrs['do_spgf']=self.compute_spgf
        grp.attrs['do_tprf'] = False
        grp.attrs['store_H'] = self.store_H
        grp.attrs['nb_sites_comp'] = self.nb_sites
        grp.create_dataset('U',(self.nb_sites,self.nb_sites,self.nb_sites,self.nb_sites,) if (not self.is_Uloc) else (self.nb_sites,) ,dtype=np.double,data=self.U_matrix)
        grp.attrs['is_Uloc']=self.is_Uloc
        grp.create_dataset('t_matrix',(self.nb_sites,self.nb_sites,),dtype=np.double,data=self.t_matrix)
#        grp.create_dataset('J',(self.nb_sites,self.nb_sites,),dtype=np.double,data=self.J_matrix)
        grp.attrs['T']=self.T
        grp.attrs['nb_comp_states']=self.nb_comp_states_
        grp.attrs['excited_state']=0
        grp_lcz=grp.create_group('lanczos')
        grp_lcz.attrs['max_lcz']=self.max_lcz
        grp_lcz.attrs['acc_lcz']=self.acc_lcz
        grp_lcz.attrs['renorm']=False
        file.close()
