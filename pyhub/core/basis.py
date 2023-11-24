import numpy as np
import h5py
import os
from pyhub.tools.tools import find_file, count_bit
from itertools import product
import re

class Basis():

    def __init__(self,nb_sites:int,order:str='spin',hilbert:tuple=None,memspace='10Go'):
        self.nb_sites = nb_sites
        self.fock=True
        self.nb_elec = None
        self.sz = None
        self.exec = find_file('../','basis.x')
        self.path = self.exec.replace('basis.x', '')
        size,unit = float(re.sub(r'(\d+)\w?o?',r'\1',memspace)),re.sub(r'\d+(\w?o?)',r'\1',memspace)
        if unit == '':
            self.memspace = size
        elif unit == 'Ko':
            self.memspace = size*10**3
        elif unit == 'Mo':
            self.memspace = size*10**6
        elif unit == 'Go':
            self.memspace = size*10**9
        elif unit == 'To':
            self.memspace = size*10**12
        else:
            raise ValueError(f'Set memspace as a correct string, not {memspace}')
        if isinstance(hilbert,tuple):
            self.fock = False
            self.nb_elec = hilbert[0]
            self.sz = hilbert[1]
            self._check_hilbert_basis(self.nb_elec,self.sz)
        if order != 'site' and order != 'spin':
            raise ValueError(f'basis order must be in site or spin convention, not {order}')
        self.order = order
        self.comp_basis = True
        if os.path.exists(self.filename_basis):
            with h5py.File(self.filename_basis,"a") as file:
                for k in file.keys():
                    if file[f'{k}/input'].attrs['nb_sites'] == self.nb_sites and \
                        file[f'{k}/input'].attrs['nb_elec'] == self.nb_elec and \
                            abs(file[f'{k}/input'].attrs['sz']-self.sz) < 1.e-14:
                        self.comp_basis = False
                        self.index = int(k)
            if self.comp_basis :
                i=0
                with h5py.File(self.filename_basis,"a") as file:
                    list_index = np.array([*file.keys()],dtype=int)
                while self.file_size > self.memspace and i <= len(list_index):
                    self.__delitem__(list_index[i])
                    i+=1
                index = np.sort(list_index)
                if len(index)>0:
                    self.index = index[0]-1 if index[0]>1 else index[-1]+1
                else:
                    self.index=1
        else:
            self.index = 1
        if self.comp_basis :
            self.compute_basis()

    def __basisremove__(self):
        try:
            os.remove(self.filename_basis)
        except:
            pass

    def __delitem__(self,k):
        with h5py.File(self.filename_basis,"a") as file:
            k_str = str(k)
            while len(k_str)<3:
                k_str = '0'+k_str
            file.__delitem__(k_str)
        return

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
        file.attrs['size_index']=len(str(self.index))
        file.attrs['index']=self.index
        grp0 = file.create_group(f'{self.index_str}')
        grp=grp0.create_group(f'input')
        grp.attrs['nb_sites'] = self.nb_sites
        grp.attrs['nb_elec'] = self.nb_elec if self.nb_elec is not None else -1
        grp.attrs['sz'] = self.sz if self.sz is not None else -1
        grp.attrs['fock'] = self.fock
        file.close()
        self.exec_basis()

        
    def _check_hilbert_basis(self,nb_elec,sz):
        if not 0<=nb_elec<=2*self.nb_sites :
            raise ValueError(f"Number of electrons {nb_elec} does not match with the number of orbitals {self.nb_sites}")
        if not (nb_elec+int(np.around(2*sz)))%2==0:
            raise ValueError(f"Number of electrons {nb_elec} does not match with Sz {sz}")

    def exec_basis(self):
        os.system(self.exec)
        with h5py.File(self.filename_basis,"a") as file:
            file[f'{self.index_str}/input'].attrs['do_solution']=True

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


    def hilbert_restricted(self,hilbert):
        number_of_elecs,sz = hilbert
        nup,ndown = int(number_of_elecs/2+sz),int(number_of_elecs/2-sz)
        self.uplist = []
        for up in self.basis_up:
            if count_bit(up) == nup:
                self.uplist.append(up)
        self.downlist = []
        for down in self.basis_down:
            if count_bit(down) == ndown:
                self.downlist.append(down)
        return self.hilbert_index

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
    def index_str(self):
        index_str = str(self.index)
        while len(index_str)<3:
            index_str = '0'+index_str
        return index_str

    @property
    def hilbert(self):
        return True if not self.fock else False

    @property
    def file_size(self):
        return float(os.path.getsize(self.filename_basis))

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
                return np.array([np.sum(self.spin_ordering(i,j)) for i,j in product(np.array(file[f'{self.index_str}/basis/basis_up'],dtype=int),np.array(file[f'{self.index_str}/basis/basis_down'],dtype=int))],dtype=int)
            else:
                return np.array([np.sum(self.site_ordering(i,j)) for i,j in product(np.array(file[f'{self.index_str}/basis/basis_up'],dtype=int),np.array(file[f'{self.index_str}/basis/basis_down'],dtype=int))],dtype=int)
    @property 
    def basis_up(self):
        with h5py.File(self.filename_basis,  "r") as file:
            return np.array(file[f'{self.index_str}/basis/basis_up'],dtype=int)

    @property 
    def basis_down(self):
        with h5py.File(self.filename_basis,  "r") as file:
            return np.array(file[f'{self.index_str}/basis/basis_down'],dtype=int)

    @property
    def nb_states(self):
        with h5py.File(self.filename_basis,  "r") as file:
            return self.nb_sup * self.nb_sdown
    @property
    def nb_sup(self):
        with h5py.File(self.filename_basis,  "r") as file:
            return file[f'{self.index_str}/basis'].attrs['nsup'][0]
    @property
    def nb_sdown(self):
        with h5py.File(self.filename_basis,  "r") as file:
            return file[f'{self.index_str}/basis'].attrs['nsdown'][0]

    @property
    def nstates(self):
        return self.nb_states
    @property
    def nsup(self):
        return self.nb_sup
    @property
    def nsdown(self):
        return self.nb_sdown
