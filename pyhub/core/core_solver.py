#!/usr/bin/env python3
import numpy as np
import h5py
import os
from time import perf_counter as pc
from copy import deepcopy as dc
from pyhub.core.reduced_quantities import GreenFunction, StaticQuantities
from pyhub.core.basis import Basis
from pyhub.tools.tools import find_file


__author__ = "Quentin Marécat"
__maintainer__ = "Quentin Marécat"
__email__ = "quentin.marecat@gmail.com"
__date__ = "January, 2023"

class CoreSolver(Basis,GreenFunction,StaticQuantities):
    def __init__(self,nb_sites,n_up:int,n_down:int,exec_name,T=0.,order='spin'): 
        Basis.__init__(self,nb_sites,hilbert=(n_up,n_down),order=order)
        self.T=T
        self.exec_name = exec_name
        self.t_matrix = np.zeros((nb_sites,nb_sites))
        self.J_matrix = np.zeros((nb_sites,nb_sites))
        self.U = 0.
        self.is_solve = False
        self.is_reduced_quantities = False
        self.is_spgf = False
        self.is_tprf = False
        self.excited_state = 0

    def kernel(self,max_lcz:int=600,acc_lcz:float=1.e-8,nb_comp_states=1,compute_rq=True,compute_spgf=False,nb_sites_comp=None,store_H=False,verbose=False):
        self.max_lcz = max_lcz
        self.acc_lcz = acc_lcz
        self.compute_rq = compute_rq
        self.nb_comp_states_ = nb_comp_states
        self.compute_spgf = compute_spgf
        self.compute_tprf = False
        self.store_H = store_H
        self.verbose = verbose
        self.nb_sites_comp = self.nb_sites if nb_sites_comp is None else nb_sites_comp

        self.exec = find_file('../../',self.exec_name)
        self.path = self.exec.replace(self.exec_name, '')
        tall=pc()
        self.__fileremove__()
        self.write_in()
        self.printv(f'start calculation')
        self._exec()
        with h5py.File(self.filename,'a') as file:
            file['input'].attrs['do_solution'] = False
        self.wgt = dc(self.boltzmann)
        if compute_rq:
            self.reduced_quantities()
        if compute_spgf:
            self.spgf(nb_sites_comp=self.nb_sites_comp)
        self.time=pc()-tall
        self.printv(f'total elapsed time {np.around(self.time,4)} s\n')

    def set_state(self,excited_state=-1):
        self.excited_state = excited_state
        self.wgt = dc(self.boltzmann)
        if excited_state == -1:
            pass
        else:
            self.wgt[:] = 0.
            self.wgt[excited_state] = 1.


    def reduced_quantities(self):
        t0=pc()
        self.printv(f'start reduced quantities calculation')
        with h5py.File(self.filename,"a") as file:
            if 'reduced_quantities' in file.keys():
                file.__delitem__('reduced_quantities')
            file['input'].attrs['excited_state'] = self.excited_state+1
            file['input'].attrs['do_rq']=True
        self._exec()
        self.printv(f'end calculation | elapsed time {np.around(pc()-t0,4)} s') 


    def spgf(self,nb_sites_comp=None):
        for nup in [-1,1]:
            Basis(self.nb_sites,([self.n_up[0]+nup],self.n_down))
        for ndown in [-1,1]:
            Basis(self.nb_sites,(self.n_up,[self.n_down[0]+ndown]))
        with h5py.File(self.filename_basis,'a') as file:
            file.attrs['lenelem2comp'] = 2
            file.attrs['elem2comp'] = [self.n_up,self.n_down]
        self.nb_sites_comp = self.nb_sites if nb_sites_comp is None else nb_sites_comp
        t0=pc()
        self.printv(f'start spgf calculation')
        with h5py.File(self.filename,"a") as file:
            if 'spgf' in file.keys():
                file.__delitem__('spgf')
            file['input'].attrs['excited_state'] = self.excited_state+1
            file['input'].attrs['do_spgf']=True
            file['input'].attrs['nb_sites_comp'] = self.nb_sites_comp
        self._exec() 
        self.printv(f'end calculation | elapsed time {np.around(pc()-t0,4)} s')

        self.gf = self.compute_gf()

    def tprf(self):
        pass


    def _exec(self):
        with h5py.File(self.filename_basis,'a') as file:
            file.attrs['lenelem2comp'] = 2
            file.attrs['elem2comp'] = [self.n_up,self.n_down]
        os.system(self.exec)

    def printv(self,*args):
        if self.verbose:
            print(*args)


    def write_in(self):
        file=h5py.File(self.filename,'a')
        grp=file.create_group('input')
        grp.attrs['nb_sites'] = self.nb_sites
        grp.attrs['n_up'] = self.n_up
        grp.attrs['n_down'] = self.n_down
        grp.attrs['do_solution']=True
        grp.attrs['do_rq']=False
        grp.attrs['do_spgf']=False
        grp.attrs['do_tprf'] = False
        grp.attrs['store_H'] = self.store_H
        grp.attrs['nb_sites_comp'] = self.nb_sites
        # if isinstance(self.U,(float,int)):
        #     grp.attrs['is_Uloc']=1
        #     grp.attrs['U'] = np.float64(self.U)
        if isinstance(self.U,(float,int)):
            grp.attrs['is_Uloc']=2
            grp.create_dataset('U',(self.nb_sites,),dtype=np.float64,data=np.full(self.nb_sites,self.U,dtype=np.float64))
            grp.attrs['do_rq_two_body']=False
        elif isinstance(self.U,(np.ndarray,list)):
            if len(np.shape(self.U))==1:
                grp.attrs['is_Uloc']=2
                grp.create_dataset('U',(self.nb_sites,),dtype=np.float64,data=self.U)
                grp.attrs['do_rq_two_body']=False
            elif len(np.shape(self.U))==4:
                if not np.isclose(self.U,np.einsum('ijkl->jilk',self.U)).all():
                    raise ValueError(f'U tensor must be symmetric')
                grp.attrs['is_Uloc']=3
                grp.create_dataset('U',(self.nb_sites,self.nb_sites,self.nb_sites,self.nb_sites,),dtype=np.float64,data=self.U)
                grp.create_dataset('Upartial',(self.nb_sites,self.nb_sites,),dtype=np.float64,data=np.einsum('ijkl,ijkl->ij',self.U,self.U))
                grp.attrs['do_rq_two_body']=True
            else :
                raise ValueError(f'U must be a float or a np.ndarray of size (nb_sites,) or (self.nb_sites,self.nb_sites,self.nb_sites,self.nb_sites,), not {type(self.U)}')
        else :
            raise ValueError(f'U must be a float or a np.ndarray of size (nb_sites,) or (self.nb_sites,self.nb_sites,self.nb_sites,self.nb_sites,), not {type(self.U)}')
        grp.create_dataset('t_matrix',(self.nb_sites,self.nb_sites,),dtype=np.double,data=self.t_matrix)
        grp.create_dataset('J',(self.nb_sites,self.nb_sites,),dtype=np.double,data=self.J_matrix)
#        grp.create_dataset('J',(self.nb_sites,self.nb_sites,),dtype=np.double,data=self.J_matrix)
        grp.attrs['T']=self.T
        grp.attrs['nb_comp_states']=self.nb_comp_states_
        grp.attrs['excited_state']=0
        grp_lcz=grp.create_group('lanczos')
        grp_lcz.attrs['max_lcz']=self.max_lcz
        grp_lcz.attrs['acc_lcz']=self.acc_lcz
        grp_lcz.attrs['renorm']=False
        file.close()

    @property
    def filename(self):
        return 'solver.h5'

    @property
    def boltzmann(self):
        with h5py.File(self.filename,"r") as file:
            return np.array(file['solve/boltzmann'],dtype=float)
    @property
    def Z(self):
        return np.sum(self.wgt)


    @property
    def nb_deg(self):
        with h5py.File(self.filename,"r") as file:
            return file['solve'].attrs['deg'][0]
    @property
    def H(self):
        with h5py.File(self.filename,"r") as file:
            return np.array(file['solve/H']).T
    @property
    def e(self):
        with h5py.File(self.filename,"r") as file:
            return np.array(file['solve/eigenvalues'])
    @property 
    def etot(self):
        return self.e0
    @property 
    def e0(self):
        return self.e[:self.nb_comp_states]@self.wgt/self.Z
    @property 
    def lcz_results(self):
        with h5py.File(self.filename,"r") as file:
            lcz={'alpha': np.array(file['solve/lanczos/alpha']),\
                        'beta' : np.array(file['solve/lanczos/beta']),\
                        'conv':file['solve/lanczos'].attrs['lcz_conv'][0],\
                        'step':file['solve/lanczos'].attrs['nb_lcz'][0]}
        if lcz['step'] > self.max_lcz:
            raise StopIteration(f'Lanczos not converged, step,', lcz['step'],f' greater than the maximum wanted {self.max_lcz}')
        return lcz
    @property
    def nb_lcz(self):
        with h5py.File(self.filename,"r") as file:
            return file['solve/lanczos'].attrs['nb_lcz'][0]
    @property
    def lcz_vectors(self):
        with h5py.File(self.filename,"r") as file:
            return np.array(file['solve/lanczos/lanczos_vectors'][:self.nb_lcz,:]).T
    @property
    def nb_comp_states(self):
        return len(self.boltzmann)

    @property
    def psi0(self):
        with h5py.File(self.filename,"r") as file:
            if not self.is_lanczos:
                return np.einsum('bi,b->i',file['solve/eigenvectors'][:self.nb_comp_states,:],self.boltzmann/self.Z)
            else:
                return np.einsum('ia,bi,b->a',np.array(file['solve/lanczos/lanczos_vectors'][:self.nb_lcz,:]),file['solve/eigenvectors'][:self.nb_comp_states,:],self.boltzmann/self.Z)

    def psi(self,vector:int):
        with h5py.File(self.filename,"r") as file:
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

    # @property 
    # def q_vector(self):
    #     return np.unique([self.ek-ek for ek in self.ek])
    # @property 
    # def q_vector_index(self):
    #     indx = np.zeros((len(self.q_vector),self.nb_sites))
    #     for i,q in enumerate(self.q_vector):
    #         for j,k in enumerate(self.ek):
    #             try:
    #                 indx[i,j]=np.where(abs(k+q-self.ek) < 1.e-12)[0][0]
    #             except:
    #                 indx[i,j]=-1
    #     return np.array(indx,dtype=int)


    # @property 
    # def chi0_q(self):
    #     return np.array([Pole({'positions':np.array([self.ek[j]-self.ek[self.q_vector_index[i,j]] if self.q_vector_index[i,j] != -1 else 0 for j in range(self.nb_sites)]),\
    #         'weights':-np.array([(self.eta_k0['up'][j]+self.eta_k0['down'][j])*(1-(self.eta_k0['up'][self.q_vector_index[i,j]]+self.eta_k0['down'][self.q_vector_index[i,j]])) if self.q_vector_index[i,j] != -1 else 0. for j in range(self.nb_sites)])},sign_eta=-1.)+\
    #         Pole({'positions':np.array([self.ek[j]-self.ek[self.q_vector_index[i,j]] if self.q_vector_index[i,j] != -1 else 0. for j in range(self.nb_sites)]),\
    #         'weights':np.array([(self.eta_k0['up'][self.q_vector_index[i,j]]+self.eta_k0['down'][self.q_vector_index[i,j]])*(1-(self.eta_k0['up'][j]+self.eta_k0['down'][j])) if self.q_vector_index[i,j] != -1 else 0. for j in range(self.nb_sites)])},sign_eta=1.) \
    #         for i,q in enumerate(self.q_vector)],dtype=Pole)*2/self.nb_sites**2 \

    # def chi_RPA_q(self,wgrid,eta=0.05):
    #     chi0_q = self.chi0_q
    #     return np.array([chi0_q[i](wgrid,eta=eta)/(1+self.U*chi0_q[i](wgrid,eta=eta)/2)\
    #         for i in range(len(self.q_vector))]).T
        


