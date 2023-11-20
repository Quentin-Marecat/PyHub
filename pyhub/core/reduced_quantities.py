import numpy as np 
import h5py
from pyhub.tools.pole import Pole
from scipy.misc import derivative

class StaticQuantities():

    def  __init__(self):
        pass

    def compute_one_rdm(self):
        with h5py.File(self.filename_hubbard,"r") as file:
            return {'up':np.einsum('jki,i->jk',file['reduced_quantities/density_matrix_up'],self.wgt/np.sum(self.wgt)),\
                'down':np.einsum('jki,i->jk',file['reduced_quantities/density_matrix_down'],self.wgt/np.sum(self.wgt))}
    def compute_ni(self):
        with h5py.File(self.filename_hubbard,"r") as file:
            return np.einsum('ji,i->j',file['reduced_quantities/ni'],self.wgt/np.sum(self.wgt))

    @property
    def two_rdm(self):
        with h5py.File(self.filename_hubbard,"r") as file:
            return np.einsum('ijklb,b->lkji',file['reduced_quantities/dbl_occ'],self.wgt/np.sum(self.wgt))

    @property
    def density(self):
        return {spin:np.diag(self.one_rdm[spin]) for spin in ['up','down']}
    @property
    def density_matrix(self):
        return self.one_rdm
    @property
    def one_rdm0(self):
        return {'up':self.Vk@np.diag([1. if i<self.nup else 0. for i in range(self.nb_sites)])@self.Vk.T,\
            'down':self.Vk@np.diag([1. if i<self.ndown else 0. for i in range(self.nb_sites)])@self.Vk.T}
    @property
    def density_matrix0(self):
        return self.one_rdm0
    @property 
    def wi(self):
        return self.ni
    @property
    def spin_cor_matrix(self):
        with h5py.File(self.filename_hubbard,"r") as file:
            return np.einsum('jki,i->jk',file['reduced_quantities/sisj'],self.wgt/np.sum(self.wgt))
    @property
    def s2(self):
        return np.sum(self.spin_cor_matrix)

    @property
    def eta_k(self):
        return {'up':np.diag(self.Vk.T@self.density_matrix['up']@self.Vk),\
            'down':np.diag(self.Vk.T@self.density_matrix['down']@self.Vk)}
    
    @property
    def eta_k0(self):
        return {'up':np.array([1. if i<self.nup else 0. for i in range(self.nb_sites)]),\
                'down':np.array([1. if i<self.ndown else 0. for i in range(self.nb_sites)])}

    @property
    def e_one_body(self):
        return np.einsum('ij,ij->',self.one_rdm['up'],self.t_matrix)+np.einsum('ij,ij->',self.one_rdm['down'],self.t_matrix)

    @property
    def e_two_body(self):
        with h5py.File(self.filename_hubbard,"r") as file:
            if file['input'].attrs['is_Uloc']:
                return np.sum(self.ni*self.U)+ np.einsum('ijb,ij,b',file['reduced_quantities/sisj'],self.J_matrix,self.wgt/np.sum(self.wgt))
            else:
                return np.einsum('lkjib,ijkl,b',file['reduced_quantities/dbl_occ'],self.U_matrix,self.wgt/np.sum(self.wgt)) \
                + np.einsum('ijb,ij,b',file['reduced_quantities/sisj'],self.J_matrix,self.wgt/np.sum(self.wgt))




class GreenFunction():

    def __init__(self):
        pass 

    def compute_gf(self):
        with h5py.File(self.filename_hubbard,"r") as file:
            return {spin:np.array([[Pole.sum([\
                Pole({'positions':-self.e[b]+np.array(file[f'spgf/positions_greater_{spin}'],dtype=float), \
                    'weights':np.einsum('k,k->k',file[f'spgf/Q_greater_{spin}'][i,:self.nb_poles[1,0 if spin =='up' else 1],b],file[f'spgf/Q_greater_{spin}'][j,:self.nb_poles[1,0 if spin =='up' else 1],b])*w/np.sum(self.wgt)}).clean(wtol=1.e-8) +
                    Pole({'positions':self.e[b]-np.array(file[f'spgf/positions_lesser_{spin}'],dtype=float), \
                    'weights':np.einsum('k,k->k',file[f'spgf/Q_lesser_{spin}'][i,:self.nb_poles[0,0 if spin =='up' else 1],b],file[f'spgf/Q_lesser_{spin}'][j,:self.nb_poles[0,0 if spin =='up' else 1],b])*w/np.sum(self.wgt)}).clean(wtol=1.e-8)\
                for b,w in enumerate(self.wgt)]) for i in range(self.nb_sites_comp)] for j in range(self.nb_sites_comp)],dtype = Pole) for spin in ['up','down']}


    @property 
    def nb_poles(self):
        with h5py.File(self.filename_hubbard,"r") as file:
            return np.array(file['spgf'].attrs['nb_poles'],dtype=int).T

    @property 
    def ip(self):
        with h5py.File(self.filename_hubbard,"r") as file:
            return {spin:(file[f'spgf/positions_lesser_{spin}'][0]-self.e0) for spin in ['up','down']}
    @property 
    def ae(self):
        with h5py.File(self.filename_hubbard,"r") as file:
            return {spin:(self.e0-file[f'spgf/positions_greater_{spin}'][0]) for spin in ['up','down']}
    @property 
    def mu(self):
        return {spin:(self.ip[spin] - self.ae[spin])/2 for spin in ['up','down']}

    @property 
    def density_matrix_from_gf(self):
        return self.one_rdm_from_gf

    @property 
    def one_rdm_from_gf(self):
        with h5py.File(self.filename_hubbard,"r") as file:
            if self.nup == self.ndown:
                A = np.einsum('ikb,jkb,b->ij',file[f'spgf/Q_lesser_up'],file[f'spgf/Q_lesser_up'],self.wgt/np.sum(self.wgt))
                return {'up':A,'down':A}
            else:
                return {'up':np.einsum('ikb,jkb,b->ij',file[f'spgf/Q_lesser_up'],file[f'spgf/Q_lesser_up'],self.wgt/np.sum(self.wgt)),\
                'down':np.einsum('ikb,jkb,b->ij',file[f'spgf/Q_lesser_down'],file[f'spgf/Q_lesser_down'],self.wgt/np.sum(self.wgt))}
            

    @property
    def gf0(self):
        return {'up':np.array([[Pole({'positions':self.ek,'weights':self.Vk[i,:]*self.Vk[j,:]}) \
                        for i in range(self.nb_sites)] for j in range(self.nb_sites)]),\
                'down':np.array([[Pole({'positions':self.ek,'weights':self.Vk[i,:]*self.Vk[j,:]}) \
                        for i in range(self.nb_sites)] for j in range(self.nb_sites)])}

    def self_energy(self,w,eta=0.05):
        if isinstance(w,(list,np.ndarray)):
            A = [self.self_energy(w_,eta=eta) for w_ in w]
            with h5py.File(self.filename_hubbard,"r") as file:
                return {'up':np.array([[[A[k]['up'][i,j] for k in range(len(w))] for i in range(file['input'].attrs['nb_sites_comp'])] for j in range(file['input'].attrs['nb_sites_comp'])],dtype=np.complex128),\
                    'down':np.array([[[A[k]['down'][i,j] for k in range(len(w))] for i in range(file['input'].attrs['nb_sites_comp'])] for j in range(file['input'].attrs['nb_sites_comp'])],dtype=np.complex128)}
        with h5py.File(self.filename_hubbard,"r") as file:
            if self.nup == self.ndown:
                A = np.linalg.pinv(np.einsum('ik,jk,k->ij',self.Vk,self.Vk,1./(w-self.ek+1j*eta))) - \
                    np.linalg.pinv(sum(\
                    np.einsum('ik,jk,k->ij',file[f'spgf/Q_greater_up'][:,:self.nb_poles[1,0],b],file[f'spgf/Q_greater_up'][:,:self.nb_poles[1,0],b],1./(w-(-self.e[b]+np.array(file[f'spgf/positions_greater_up'],dtype=float))+1j*eta))*bltz/sum(self.wgt) + \
                    np.einsum('ik,jk,k->ij',file[f'spgf/Q_lesser_up'][:,:self.nb_poles[0,0],b],file[f'spgf/Q_lesser_up'][:,:self.nb_poles[0,0],b],1./(w-(self.e[b]-np.array(file[f'spgf/positions_lesser_up'],dtype=float))+1j*eta))*bltz/sum(self.wgt)\
                    for b,bltz in enumerate(self.wgt)))
                return {'up':A,'down':A}
            else:
                return {'up':\
                    np.linalg.pinv(np.einsum('ik,jk,k->ij',self.Vk,self.Vk,1./(w-self.ek+1j*eta))) - \
                    np.linalg.pinv(sum(\
                    np.einsum('ik,jk,k->ij',file[f'spgf/Q_greater_up'][:,:self.nb_poles[1,0],b],file[f'spgf/Q_greater_up'][:,:self.nb_poles[1,0],b],1./(w-(-self.e[b]+np.array(file[f'spgf/positions_greater_up'],dtype=float))+1j*eta)) + \
                    np.einsum('ik,jk,k->ij',file[f'spgf/Q_lesser_up'][:,:self.nb_poles[0,0],b],file[f'spgf/Q_lesser_up'][:,:self.nb_poles[0,0],b],1./(w-(self.e[b]-np.array(file[f'spgf/positions_lesser_up'],dtype=float))+1j*eta))\
                    for b in range(self.nb_comp_states))),\
                    'down':\
                    np.linalg.pinv(np.einsum('ik,jk,k->ij',self.Vk,self.Vk,1./(w-self.ek+1j*eta))) - \
                    np.linalg.pinv(sum(\
                    np.einsum('ik,jk,k->ij',file[f'spgf/Q_greater_down'][:,:self.nb_poles[1,1],b],file[f'spgf/Q_greater_down'][:,:self.nb_poles[1,1],b],1./(w-(-self.e[b]+np.array(file[f'spgf/positions_greater_down'],dtype=float))+1j*eta)) + \
                    np.einsum('ik,jk,k->ij',file[f'spgf/Q_lesser_down'][:,:self.nb_poles[0,1],b],file[f'spgf/Q_lesser_down'][:,:self.nb_poles[0,1],b],1./(w-(self.e[b]-np.array(file[f'spgf/positions_lesser_down'],dtype=float))+1j*eta))\
                    for b in range(self.nb_comp_states)))\
                        }

    @property
    def galitskii_migdal(self):
#        integrate(lambda w:(-1/np.pi)*(_gf(w)).imag*distribution_function(w,self.pot_chim,T=self.T)*w,bounds=self.bounds,name=self.name,args=self.args)
        pass
    @property 
    def qpw(self)->float:
        return np.linalg.pinv(np.identity(self.nb_sites) - np.real(derivative(lambda w:self.self_energy(w,eta=0.0001)['up'],self.mu['up'],dx=1.e-6,order=3)))