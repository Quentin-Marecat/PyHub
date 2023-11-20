import numpy as np
import h5py
import os
from pyhub.tools.tools import find_file
from itertools import product

class Basis():

    def __init__(self,nb_sites:int,nb_elec:int,sz:float,order:str='spin',fock:bool=False,init=False):
        self.nb_sites = nb_sites
        self.nb_elec = nb_elec
        self.sz = sz 
        self.fock=False
        if fock or (nb_elec == None and sz == None):
            self.fock = True
            self.nb_elec = 0
            self.sz = 0.
        if order != 'site' and order != 'spin':
            raise ValueError(f'basis order must be in site or spin convention, not {order}')
        self.order = order
        self.comp_basis = True
        if not self.fock:
            self._check_in_basis(self.nb_elec,self.sz)
        self.exec = find_file('../','basis.x')
        self.path = self.exec.replace('basis.x', '')
        if not init:
            if os.path.exists(self.filename_basis):
                with h5py.File(self.filename_basis,  "a") as file:
                    if file['input'].attrs['nb_sites'] == self.nb_sites and file['input'].attrs['nb_elec'] == self.nb_elec and abs(file['input'].attrs['sz']-self.sz) < 1.e-14:
                        self.comp_basis = False
        if self.comp_basis :
            self.__basisremove__()
            self.compute_basis()

    def __basisremove__(self):
        try:
            os.remove(self.filename_basis)
        except FileNotFoundError:
            pass

    def __getitem__(self, index):
        if isinstance(index, np.ndarray):
            index = index.astype(int).tolist()
        if index == slice(None, None, None):
            return self.basis[:]
        else:
            return self.basis[index]

    def __iter__(self):
        return iter(self.basis)

    def __repr__(self):
        if self.order=='spin':
            return f'basis up :\n{[up for up in self.basis_up]}\nbasis down :\n{[down<<self.nb_sites for down in self.basis_down]}'
        elif self.order=='site':
            uplist = []
            for up in self.basis_up:
                a_new = 0
                for i in range(self.nb_sites):
                    bit_a = (up >> i) & 1  # Obtenez le bit à la position i
                    a_new |= bit_a << (2*i)  # Placez le bit à la position 2*i+1
                uplist.append(a_new)
            downlist = []
            for down in self.basis_down:
                a_new = 0
                for i in range(self.nb_sites):
                    bit_a = (down >> i) & 1  # Obtenez le bit à la position i
                    a_new |= bit_a << (2*i+1)  # Placez le bit à la position 2*i+1
                downlist.append(a_new)
            return f'basis up :\n{uplist}\nbasis down :\n{downlist}'

    def compute_basis(self):
        file=h5py.File(self.filename_basis,'a')
        grp=file.create_group('input')
        grp.attrs['nb_sites'] = self.nb_sites
        grp.attrs['nb_elec'] = self.nb_elec 
        grp.attrs['sz'] = self.sz
        grp.attrs['fock'] = self.fock
        file.close()
        self.exec_basis()

        
    def _check_in_basis(self,nb_elec,sz):
        if not 0<=nb_elec<=2*self.nb_sites :
            raise ValueError(f"Number of electrons {nb_elec} does not match with the number of orbitals {self.nb_sites}")
        if not (nb_elec+int(np.around(2*sz)))%2==0:
            raise ValueError(f"Number of electrons {nb_elec} does not match with Sz {sz}")

    def exec_basis(self):
        os.system(self.exec)
        with h5py.File(self.filename_basis,"a") as file:
            file['input'].attrs['do_solution']=True

    def spin_ordering(self,a,b):
        return a,(b<<self.nb_sites)

    def site_ordering(self,a,b):
        a_new,b_new = 0,0
        for i in range(self.nb_sites):
            bit_a = (a >> i) & 1  # Obtenez le bit à la position i
            bit_b = (b >> i) & 1  # Obtenez le bit à la position i
            a_new |= bit_a << (2*i)  # Placez le bit à la position 2*i+1
            b_new |= bit_b << (2*i+1)  # Placez le bit à la position 2*i+1
        return a_new,b_new

    def count_bit(self,n):
        '''
        Count the number of ones in bitstring given by integer n
        '''
        count = 0
        while n:
            count += 1
            n &= (n - 1)
        return count

    def hilbert_restricted(self,number_of_elecs=None,sz=None):
        number_of_elecs = self.nb_elec if number_of_elecs is None else int(number_of_elecs)
        sz = self.sz if sz is None else sz
        nup,ndown = int(number_of_elecs/2+sz),int(number_of_elecs/2-sz)
        self.uplist = []
        for up in self.basis_up:
            if self.count_bit(up) == nup:
                self.uplist.append(up)
        self.downlist = []
        for down in self.basis_down:
            if self.count_bit(down) == ndown:
                self.downlist.append(down)

    @property
    def hilbert_index(self):
        with h5py.File(self.filename_basis,  "r") as file:
            if self.order == 'spin':
                elem = np.array([np.sum(self.spin_ordering(i,j)) for i,j in product(self.uplist,self.downlist)],dtype=int)
            else:
                elem = np.array([np.sum(self.site_ordering(i,j)) for i,j in product(self.uplist,self.downlist)],dtype=int)
        return np.where(np.isin(self.basis, np.intersect1d(self.basis,elem)))[0]


    @property
    def filename_basis(self):
        return 'basis.h5'

    @property
    def hilbert(self):
        return True if not self.fock else False

    @property
    def nup(self):
        int(float(self.nb_elec)/2 + self.sz)
    @property
    def ndown(self):
        int(float(self.nb_elec)/2-self.sz)

    @property 
    def basis(self):
        with h5py.File(self.filename_basis,  "r") as file:
            if self.order == 'spin':
                return np.array([np.sum(self.spin_ordering(i,j)) for i,j in product(np.array(file['basis/basis_up'],dtype=int),np.array(file['basis/basis_down'],dtype=int))],dtype=int)
            else:
                return np.array([np.sum(self.site_ordering(i,j)) for i,j in product(np.array(file['basis/basis_up'],dtype=int),np.array(file['basis/basis_down'],dtype=int))],dtype=int)
    @property 
    def basis_up(self):
        with h5py.File(self.filename_basis,  "r") as file:
            return np.array(file['basis/basis_up'],dtype=int)
    @property 
    def basis_down(self):
        with h5py.File(self.filename_basis,  "r") as file:
            return np.array(file['basis/basis_down'],dtype=int)

    @property
    def nb_states(self):
        with h5py.File(self.filename_basis,  "r") as file:
            return file['basis'].attrs['nstates'][0]
    @property
    def nb_sup(self):
        with h5py.File(self.filename_basis,  "r") as file:
            return file['basis'].attrs['nsup'][0]
    @property
    def nb_sdown(self):
        with h5py.File(self.filename_basis,  "r") as file:
            return file['basis'].attrs['nsdown'][0]

    @property
    def nstates(self):
        return self.nb_states
    @property
    def nsup(self):
        return self.nb_sup
    @property
    def nsdown(self):
        return self.nb_sdown
