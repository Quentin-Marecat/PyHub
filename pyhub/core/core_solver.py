#!/usr/bin/env python3
import numpy as np
import h5py
import os
from time import perf_counter as pc
from copy import deepcopy as dc
from pyhub.core.reduced_quantities import GreenFunction, StaticQuantities
from pyhub.tools.pole import Pole
from pyhub.core.basis import Basis
from pyhub.tools.tools import find_file


__author__ = "Quentin Marécat"
__maintainer__ = "Quentin Marécat"
__email__ = "quentin.marecat@gmail.com"
__date__ = "January, 2023"

class CoreSolver(Basis,GreenFunction,StaticQuantities):
    def __init__(self,nb_sites,nb_elec:int,sz:float,T=0.,fock=False,order='spin'): 
        Basis.__init__(self,nb_sites,hilbert=(nb_elec,sz) if not fock else None,order=order)
        self.T=T
        self.t_matrix = np.zeros((nb_sites,nb_sites))
        self.J_matrix = np.zeros((nb_sites,nb_sites))
        self.U = 0.
        self.is_solve = False
        self.is_reduced_quantities = False
        self.is_spgf = False
        self.is_tprf = False


    def kernel(self,max_lcz:int=600,acc_lcz:float=1.e-8,nb_comp_states=1,compute_rq_hubbard=True,compute_two_body=False,compute_spgf=False,nb_sites_comp=None,store_H=False,verbose=False):
        self.max_lcz = max_lcz
        self.acc_lcz = acc_lcz
        self.compute_rq_hubbard = compute_rq_hubbard
        self.compute_two_body = compute_two_body
        self.nb_comp_states_ = nb_comp_states
        self.compute_spgf = compute_spgf
        self.compute_tprf = False
        self.store_H = store_H
        self.verbose = verbose

        self.exec = find_file('../../','hubbard.x')
        self.path = self.exec.replace('hubbard.x', '')
        tall=pc()
        self.__hubbardremove__()
        self.write_in_hubbard()
        self.printv(f'start calculation')
        self._exec()
        self.wgt = dc(self.boltzmann)
        self.time=pc()-tall
        if compute_rq_hubbard:
            self.one_rdm,self.ni = self.compute_one_rdm(),self.compute_ni()
        if compute_spgf:
            self.nb_sites_comp = self.nb_sites if nb_sites_comp is None else nb_sites_comp
            self.gf = self.compute_gf()
        self.printv(f'total elapsed time {np.around(self.time,4)} s\n')


    def reduced_quantities(self,excited_state=-1,two_body=False):
        t0=pc()
        self.printv(f'start reduced quantities calculation | set 2RDM {two_body}')
        with h5py.File(self.filename_hubbard,"a") as file:
            if 'reduced_quantities' in file.keys():
                file.__delitem__('reduced_quantities')
            file['input'].attrs['excited_state'] = excited_state+1
            file['input'].attrs['do_rq']=True
            file['input'].attrs['do_rq_two_body']=two_body
        self._exec()
        self.printv(f'end calculation | elapsed time {np.around(pc()-t0,4)} s') 

        self.wgt = dc(self.boltzmann)
        if excited_state != -1:
            self.wgt[:] = 0.
            self.wgt[excited_state] = 1.


    def spgf(self,excited_state=-1,nb_sites_comp=None):
        self.nb_sites_comp = self.nb_sites if nb_sites_comp is None else nb_sites_comp
        t0=pc()
        self.printv(f'start spgf calculation')
        with h5py.File(self.filename_hubbard,"a") as file:
            if 'spgf' in file.keys():
                file.__delitem__('spgf')
            file['input'].attrs['excited_state'] = excited_state+1
            file['input'].attrs['do_spgf']=True
            file['input'].attrs['nb_sites_comp'] = self.nb_sites_comp
        self._exec() 
        self.printv(f'end calculation | elapsed time {np.around(pc()-t0,4)} s')

        self.wgt = dc(self.boltzmann)
        if excited_state != -1:
            self.wgt[:] = 0.
            self.wgt[excited_state] = 1.
        self.gf = self.compute_gf()

    def tprf(self,excited_state=-1,nb_sites_comp=None):
        pass


    def _exec(self):
        os.system(self.exec)
        with h5py.File(self.filename_hubbard,"a") as file:
            file['input'].attrs['do_solution']=False

    def printv(self,*args):
        if self.verbose:
            print(*args)

    @property
    def boltzmann(self):
        with h5py.File(self.filename_hubbard,"r") as file:
            return np.array(file['solve/boltzmann'],dtype=float)
    @property
    def Z(self):
        return np.sum(self.boltzmann)


    @property
    def nb_deg(self):
        with h5py.File(self.filename_hubbard,"r") as file:
            return file['solve'].attrs['deg'][0]
    @property
    def H(self):
        with h5py.File(self.filename_hubbard,"r") as file:
            return np.array(file['solve/H']).T
    @property
    def e(self):
        with h5py.File(self.filename_hubbard,"r") as file:
            return np.array(file['solve/eigenvalues'])
    @property 
    def etot(self):
        return self.e0
    @property 
    def e0(self):
        return self.e[:self.nb_comp_states]@self.boltzmann/self.Z
    @property 
    def lcz_results(self):
        with h5py.File(self.filename_hubbard,"r") as file:
            lcz={'alpha': np.array(file['solve/lanczos/alpha']),\
                        'beta' : np.array(file['solve/lanczos/beta']),\
                        'conv':file['solve/lanczos'].attrs['lcz_conv'][0],\
                        'step':file['solve/lanczos'].attrs['nb_lcz'][0]}
        if lcz['step'] > self.max_lcz:
            raise StopIteration(f'Lanczos not converged, step,', lcz['step'],f' greater than the maximum wanted {self.max_lcz}')
        return lcz
    @property
    def nb_lcz(self):
        with h5py.File(self.filename_hubbard,"r") as file:
            return file['solve/lanczos'].attrs['nb_lcz'][0]
    @property
    def lcz_vectors(self):
        with h5py.File(self.filename_hubbard,"r") as file:
            return np.array(file['solve/lanczos/lanczos_vectors'][:self.nb_lcz,:]).T
    @property
    def nb_comp_states(self):
        return len(self.boltzmann)

    @property
    def psi0(self):
        with h5py.File(self.filename_hubbard,"r") as file:
            if not self.is_lanczos:
                return np.einsum('bi,b->i',file['solve/eigenvectors'][:self.nb_comp_states,:],self.boltzmann/self.Z)
            else:
                return np.einsum('ia,bi,b->a',np.array(file['solve/lanczos/lanczos_vectors'][:self.nb_lcz,:]),file['solve/eigenvectors'][:self.nb_comp_states,:],self.boltzmann/self.Z)

    def psi(self,vector:int):
        with h5py.File(self.filename_hubbard,"r") as file:
            if not self.is_lanczos:
                return np.array(file['solve/eigenvectors'][vector,:],dtype=float)
            else:
                return self.lcz_vectors@np.array(file['solve/eigenvectors'][vector,:],dtype=float)


    def psi_as_dict(self,vector:int):
        psi = self.psi(vector)
        psi_dict = [{'basis':[],'coeff':[]} for i in range(len(psi))]
        for i,vec in enumerate(psi):
            index = np.where(abs(vec)>1.e-12)[0]
            psi_dict[i]['basis'] = self.basis[index]
            psi_dict[i]['coeff'] = psi[index]
        return psi_dict

    @property
    def is_lanczos(self):
        try:
            self.nb_lcz
            return True
        except:
            return False


    @property
    def U_matrix(self):
        if isinstance(self.U,np.ndarray):
            if len(self.U.shape) == 4 or len(self.U.shape) == 1:
                return self.U
            else :
                raise ValueError
        if isinstance(self.U,(float,int)):
            return np.full(self.nb_sites,self.U,dtype=float)

    @property 
    def q_vector(self):
        return np.unique([self.ek-ek for ek in self.ek])
    @property 
    def q_vector_index(self):
        indx = np.zeros((len(self.q_vector),self.nb_sites))
        for i,q in enumerate(self.q_vector):
            for j,k in enumerate(self.ek):
                try:
                    indx[i,j]=np.where(abs(k+q-self.ek) < 1.e-12)[0][0]
                except:
                    indx[i,j]=-1
        return np.array(indx,dtype=int)


    @property 
    def chi0_q(self):
        return np.array([Pole({'positions':np.array([self.ek[j]-self.ek[self.q_vector_index[i,j]] if self.q_vector_index[i,j] != -1 else 0 for j in range(self.nb_sites)]),\
            'weights':-np.array([(self.eta_k0['up'][j]+self.eta_k0['down'][j])*(1-(self.eta_k0['up'][self.q_vector_index[i,j]]+self.eta_k0['down'][self.q_vector_index[i,j]])) if self.q_vector_index[i,j] != -1 else 0. for j in range(self.nb_sites)])},sign_eta=-1.)+\
            Pole({'positions':np.array([self.ek[j]-self.ek[self.q_vector_index[i,j]] if self.q_vector_index[i,j] != -1 else 0. for j in range(self.nb_sites)]),\
            'weights':np.array([(self.eta_k0['up'][self.q_vector_index[i,j]]+self.eta_k0['down'][self.q_vector_index[i,j]])*(1-(self.eta_k0['up'][j]+self.eta_k0['down'][j])) if self.q_vector_index[i,j] != -1 else 0. for j in range(self.nb_sites)])},sign_eta=1.) \
            for i,q in enumerate(self.q_vector)],dtype=Pole)*2/self.nb_sites**2 \

    def chi_RPA_q(self,wgrid,eta=0.05):
        chi0_q = self.chi0_q
        return np.array([chi0_q[i](wgrid,eta=eta)/(1+self.U*chi0_q[i](wgrid,eta=eta)/2)\
            for i in range(len(self.q_vector))]).T
        


