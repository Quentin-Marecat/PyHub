#!/usr/bin/env python3
import numpy as np
import h5py
import os
from time import perf_counter as pc
from copy import deepcopy as dc
from pyhub.core.reduced_quantities import GreenFunction, StaticQuantities
from pyhub.tools.pole import Pole


__author__ = "Quentin Marécat"
__maintainer__ = "Quentin Marécat"
__email__ = "quentin.marecat@gmail.com"
__date__ = "January, 2023"

def find_file(directory, filename):
    for root, dirs, files in os.walk(directory):
        if filename in files:
            return os.path.join(root, filename)
    raise FileNotFoundError(f"Le fichier {filename} n'a pas été trouvé dans le répertoire {directory} ni dans ses sous-répertoires.")

class CoreSolver(GreenFunction,StaticQuantities):
    def __init__(self,nb_sites,nb_elec:int,sz:float,T=0.): 
        self.nb_sites=nb_sites
        self.nb_elec=nb_elec
        self.sz=sz
        self.T=T
        self.t_matrix = np.zeros((nb_sites,nb_sites))
        self.J_matrix = np.zeros((nb_sites,nb_sites))
        self.U = 0.
        self.compute_basis = True
        self.is_basis = False
        self.is_solve = False
        self.is_reduced_quantities = False
        self.is_spgf = False
        self.is_tprf = False


    def run(self,max_lcz:int=600,acc_lcz:float=1.e-8,nb_comp_states=1,compute_reduced_quantities=True,compute_two_body=False,compute_spgf=False,compute_tprf=False,store_H=False,verbose=False):
        self.max_lcz=max_lcz
        self.acc_lcz=acc_lcz
        self.compute_reduced_quantities=compute_reduced_quantities
        self.nb_comp_states_ = nb_comp_states
        self.compute_spgf=compute_spgf
        self.compute_tprf=compute_tprf
        self.store_H=store_H
        self.verbose = verbose

        self.exec = find_file('../../','main.x')
        self.path = self.exec.replace('main.x', '')
        tall=pc()
        self._check_in()
        self.write_in()
        self.printv(f'start calculation')
        t0=pc()
        self._exec()
        self.printv(f'end calculation | elapsed time {np.around(pc()-t0,4)} s')
        with h5py.File(self.filename,  "a") as file:
            file['input'].attrs['do_basis']=False
            file['input'].attrs['do_solution']=False
        self.is_solve = True
        with h5py.File(self.filename,  "r") as file:
            self.boltzmann = np.array(file['solve/boltzmann'],dtype=float)
        self.Z = np.sum(self.boltzmann)
        if self.compute_reduced_quantities:
            self.reduced_quantities(two_body=compute_two_body)
        if self.compute_spgf:
            self.spgf()
        self.time=pc()-tall
        self.printv(f'total elapsed time {np.around(self.time,4)} s\n')


    def reduced_quantities(self,excited_state=-1,two_body=True):
        if not self.is_reduced_quantities:
            t0=pc()
            self.printv(f'start reduced quantities calculation | set 2RDM {two_body}')
            with h5py.File(self.filename,  "a") as file:
                if 'reduced_quantities' in file.keys():
                    file.__delitem__('reduced_quantities')
                file['input'].attrs['excited_state'] = excited_state+1
                file['input'].attrs['do_rq']=True
                file['input'].attrs['do_rq_two_body']=two_body
            self._exec()
            with h5py.File(self.filename,  "a") as file:
                file['input'].attrs['do_rq']=False
            self.is_reduced_quantities = True  
            self.printv(f'end calculation | elapsed time {np.around(pc()-t0,4)} s') 

        self.wgt = dc(self.boltzmann)
        if excited_state != -1:
            self.wgt[:] = 0.
            self.wgt[excited_state] = 1.
        with h5py.File(self.filename,  "r") as file:
            self.one_rdm = {'up':np.einsum('jki,i->jk',file['reduced_quantities/density_matrix_up'],self.wgt/np.sum(self.wgt)),\
                'down':np.einsum('jki,i->jk',file['reduced_quantities/density_matrix_down'],self.wgt/np.sum(self.wgt))}
            self.ni = np.einsum('ji,i->j',file['reduced_quantities/ni'],self.wgt/np.sum(self.wgt))
        self.compute_rq = True


    def spgf(self,excited_state=-1,nb_sites_comp=None):
        nb_sites_comp = self.nb_sites if nb_sites_comp is None else nb_sites_comp
        if not self.is_spgf:
            t0=pc()
            self.printv(f'start spgf')
            with h5py.File(self.filename,  "a") as file:
                if 'spgf' in file.keys():
                    file.__delitem__('spgf')
                file['input'].attrs['excited_state'] = excited_state+1
                file['input'].attrs['do_spgf']=True
                file['input'].attrs['nb_sites_comp'] = nb_sites_comp
            self._exec()
            with h5py.File(self.filename,  "a") as file:
                file['input'].attrs['do_spgf']=False
            self.is_spgf = True 
            self.printv(f'end calculation | elapsed time {np.around(pc()-t0,4)} s')

        self.wgt = dc(self.boltzmann)
        if excited_state != -1:
            self.wgt[:] = 0.
            self.wgt[excited_state] = 1.
        t1=pc()
        self.printv(f'  set gf as pole')
        with h5py.File(self.filename,  "r") as file:
            self.nb_poles = np.array(file['spgf'].attrs['nb_poles'],dtype=int).T
            self.gf = {spin:np.array([[Pole.sum([\
                Pole({'positions':-self.e[b]+np.array(file[f'spgf/positions_greater_{spin}'],dtype=float), \
                    'weights':np.einsum('k,k->k',file[f'spgf/Q_greater_{spin}'][i,:self.nb_poles[1,0 if spin =='up' else 1],b],file[f'spgf/Q_greater_{spin}'][j,:self.nb_poles[1,0 if spin =='up' else 1],b])*w/np.sum(self.wgt)}).clean(wtol=1.e-8) +
                    Pole({'positions':self.e[b]-np.array(file[f'spgf/positions_lesser_{spin}'],dtype=float), \
                    'weights':np.einsum('k,k->k',file[f'spgf/Q_lesser_{spin}'][i,:self.nb_poles[0,0 if spin =='up' else 1],b],file[f'spgf/Q_lesser_{spin}'][j,:self.nb_poles[0,0 if spin =='up' else 1],b])*w/np.sum(self.wgt)}).clean(wtol=1.e-8)\
                for b,w in enumerate(self.wgt)]) for i in range(nb_sites_comp)] for j in range(nb_sites_comp)],dtype = Pole) for spin in ['up','down']}
            self.ip = {spin:(file[f'spgf/positions_lesser_{spin}'][0]-self.e0) for spin in ['up','down']}
            self.ae = {spin:(self.e0-file[f'spgf/positions_greater_{spin}'][0]) for spin in ['up','down']}
            self.mu = {spin:(self.ip[spin] - self.ae[spin])/2 for spin in ['up','down']}
        self.printv(f'  end gf as pole | elapsed time {np.around(pc()-t1,4)} s')

    def tprf(self,excited_state=-1,nb_sites_comp=None):
        pass


    def _exec(self):
        os.system(self.exec)

    def printv(self,*args):
        if self.verbose:
            print(*args)

    def _check_in(self):
        if abs(self.sz)<1.e-14:
            self.sz=0.
        if not 0<=self.nb_elec<=2*self.nb_sites :
            raise ValueError(f"Number of electrons {self.nb_elec} does not match with the number of orbitals {self.nb_sites}")
        if not (self.nb_elec+int(2*self.sz))%2==0:
            raise ValueError(f"Number of electrons {self.nb_elec} does not match with Sz {self.sz}")
        if os.path.exists(self.filename):
            with h5py.File(self.filename,  "a") as file:
                if file['input'].attrs['nb_sites'] == self.nb_sites and file['input'].attrs['nb_elec'] == self.nb_elec and abs(file['input'].attrs['sz']-self.sz) < 1.e-14:
                    file.__delitem__('input')
                    self.compute_basis = False
                    self.printv(f'basis already compute: skip')
                    file2 = h5py.File('new_'+self.filename,  "w")
                    file.copy('basis',file2)
                os.remove(self.filename)
            if os.path.exists('new_'+self.filename):
                with h5py.File('new_'+self.filename,  "r") as file2:
                    file = h5py.File(self.filename,  "w")
                    file2.copy('basis',file)
                os.remove('new_'+self.filename)
        if os.path.exists(self.path+'/'+self.basis_filename) and self.compute_basis:
            with h5py.File(self.filename,'a') as file:
                grp=file.create_group('input')
                grp.attrs['nb_sites'] = self.nb_sites
                grp.attrs['nb_elec'] = self.nb_elec 
                grp.attrs['sz'] = self.sz 
                grp2=file.create_group('basis')
                with h5py.File(self.path+'/'+self.basis_filename,  "r") as file2:
                    grp2.attrs['nstates']=file2[f'{self.nb_sites}/{self.nb_elec}/{np.around(self.sz,2)}'].attrs['nb_states']
                    grp2.attrs['nup']=file2[f'{self.nb_sites}/{self.nb_elec}/{np.around(self.sz,2)}'].attrs['nup']
                    grp2.attrs['ndown']=file2[f'{self.nb_sites}/{self.nb_elec}/{np.around(self.sz,2)}'].attrs['ndown']
                    grp2.attrs['nsup']=file2[f'{self.nb_sites}/{self.nb_elec}/{np.around(self.sz,2)}'].attrs['nsup']
                    grp2.attrs['nsdown']=file2[f'{self.nb_sites}/{self.nb_elec}/{np.around(self.sz,2)}'].attrs['nsdown']
                    grp2.create_dataset("basis_up", (grp2.attrs['nsup'],), dtype='i',data=file2[f'{self.nb_sites}/{self.nb_elec}/{np.around(self.sz,2)}/basis_up'])
                    grp2.create_dataset("basis_down", (grp2.attrs['nsdown'],), dtype='i',data=file2[f'{self.nb_sites}/{self.nb_elec}/{np.around(self.sz,2)}/basis_down'])
            self._check_in()


    @property
    def nup(self):
        return int(float(self.nb_elec)/2+self.sz)
    @property
    def ndown(self):
        return int(float(self.nb_elec)/2-self.sz)

    @property 
    def basis_up(self):
        with h5py.File(self.filename,  "r") as file:
            return np.array(file['basis/basis_up'],dtype=int)
    @property 
    def basis_down(self):
        with h5py.File(self.filename,  "r") as file:
            return np.array(file['basis/basis_down'],dtype=int)
    @property 
    def basis(self):
        with h5py.File(self.filename,  "r") as file:
            return np.array(file['basis/basis_full'],dtype=int)
    @property
    def nb_states(self):
        with h5py.File(self.filename,  "r") as file:
            return file['basis'].attrs['nstates'][0]
    @property
    def nsup(self):
        with h5py.File(self.filename,  "r") as file:
            return file['basis'].attrs['nsup'][0]
    @property
    def nsdown(self):
        with h5py.File(self.filename,  "r") as file:
            return file['basis'].attrs['nsdown'][0]

    @property
    def nb_deg(self):
        with h5py.File(self.filename,  "r") as file:
            return file['solve'].attrs['deg'][0]
    @property
    def H(self):
        with h5py.File(self.filename,  "r") as file:
            return np.array(file['solve/H']).T
    @property
    def e(self):
        with h5py.File(self.filename,  "r") as file:
            return np.array(file['solve/eigenvalues'])
    @property 
    def etot(self):
        return self.e0
    @property 
    def e0(self):
        return self.e[:self.nb_comp_states]@self.boltzmann/self.Z
    @property 
    def lcz_results(self):
        with h5py.File(self.filename,  "r") as file:
            lcz={'alpha': np.array(file['solve/lanczos/alpha']),\
                        'beta' : np.array(file['solve/lanczos/beta']),\
                        'conv':file['solve/lanczos'].attrs['lcz_conv'][0],\
                        'step':file['solve/lanczos'].attrs['nb_lcz'][0]}
        if lcz['step'] > self.max_lcz:
            raise StopIteration(f'Lanczos not converged, step,', lcz['step'],f' greater than the maximum wanted {self.max_lcz}')
        return lcz
    @property
    def nb_lcz(self):
        with h5py.File(self.filename,  "r") as file:
            return file['solve/lanczos'].attrs['nb_lcz'][0]
    @property
    def lcz_vectors(self):
        with h5py.File(self.filename,  "r") as file:
            return np.array(file['solve/lanczos/lanczos_vectors'][:self.nb_lcz,:]).T
    @property
    def nb_comp_states(self):
        return len(self.boltzmann)

    @property
    def psi0(self):
        with h5py.File(self.filename,  "r") as file:
            if not self.is_lanczos:
                return np.einsum('bi,b->i',file['solve/eigenvectors'][:self.nb_comp_states,:],self.boltzmann/self.Z)
            else:
                return np.einsum('ia,bi,b->a',np.array(file['solve/lanczos/lanczos_vectors'][:self.nb_lcz,:]),file['solve/eigenvectors'][:self.nb_comp_states,:],self.boltzmann/self.Z)

    def psi(self,vector:int):
        with h5py.File(self.filename,  "r") as file:
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
        return set_U_matrix(self.nb_sites,self.U)

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
        




def set_U_matrix(nb_sites:int,l:float or int or list or np.ndarray):
    if isinstance(l,np.ndarray):
        if len(l.shape) == 4:
            return l
        if len(l.shape) == 1:
            U = np.zeros((nb_sites,nb_sites,nb_sites,nb_sites))
            for i in range(nb_sites):
                U[i,i,i,i] = l[i]
            return U
    if isinstance(l,(float,int)):
        U = np.zeros((nb_sites,nb_sites,nb_sites,nb_sites))
        for i in range(nb_sites):
            U[i,i,i,i] = l
        return U

