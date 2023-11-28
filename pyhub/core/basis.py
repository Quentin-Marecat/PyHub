import numpy as np
import h5py
import os
from pyhub.tools.tools import find_file, count_bit, int2str
from itertools import product
import re

class Basis():

    def __init__(self,nb_sites:int,hilbert:tuple=None,order:str='spin',memspace=None):
        self.nb_sites = nb_sites
        self.hilbert = hilbert
        self.exec = find_file('../../','basis.x')
        self.path = self.exec.replace('basis.x', '')
        if isinstance(memspace,str):
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
            i=0
            with h5py.File(self.filename_basis,"a") as file:
                list_index = np.array([*file.keys()],dtype=int)
            while self.file_size > self.memspace and i <= len(list_index):
                file=h5py.File(self.filename_basis,"a")
                file.__delitem__(list_index[i])
                file.close()
                i+=1
        elif memspace==None:
            pass
        else:
            raise ValueError(f'Set memspace as a correct string, not {memspace}')
        if isinstance(hilbert,tuple):
            if isinstance(hilbert[0],int):
                self.nup = [hilbert[0]]
            elif isinstance(hilbert[0],(list,np.ndarray)):
                self.nup = hilbert[0]
            else :
                raise ValueError(f'hilbert must be a tuple of n_up and n_down, not {hilbert}')
            if isinstance(hilbert[1],int):
                self.ndown = [hilbert[1]]
            elif isinstance(hilbert[1],(list,np.ndarray)):
                self.ndown = hilbert[1]
            else :
                raise ValueError(f'hilbert must be a tuple of n_up and n_down, not {hilbert}')
        elif hilbert == None : 
            self.nup = np.array([i for i in range(self.nb_sites+1)])
            self.ndown = np.array([i for i in range(self.nb_sites+1)])
        else :
            raise ValueError(f'hilbert must be a tuple of n_up and n_down, not {hilbert}')
        if order != 'site' and order != 'spin':
            raise ValueError(f'basis order must be in site or spin convention, not {order}')
        try:
            self.order = order
        except AttributeError:
            pass
        if os.path.exists(self.filename_basis):
            with h5py.File(self.filename_basis,"a") as file:
                comp_grp = True
                for nsites in file.keys():
                    if nsites == int2str(self.nb_sites,2) :
                        elem2comp = np.array([element for element in np.union1d(self.nup,self.ndown) if element not in np.array(list(file[f'{int2str(self.nb_sites,2)}'].keys()),dtype=int)],dtype=int)
                        comp_grp = False
                        break
                if comp_grp == True:
                    grp = file.create_group(f'{int2str(self.nb_sites,2)}')
                    elem2comp = np.union1d(self.nup,self.ndown)
        else:
            with h5py.File(self.filename_basis,"a") as file:
                grp = file.create_group(f'{int2str(self.nb_sites,2)}')
                elem2comp = np.union1d(self.nup,self.ndown)
        with h5py.File(self.filename_basis,'a') as file:
            file.attrs['index']=self.nb_sites
        if len(elem2comp)>0:
            self.exec_basis(elem2comp)

    def __basisremove__(self):
        try:
            os.remove(self.filename_basis)
        except:
            pass

    # def __delitem__(self,k):
    #     with h5py.File(self.filename_basis,"a") as file:
    #         k_str = str(k)
    #         while len(k_str)<3:
    #             k_str = '0'+k_str
    #         file.__delitem__(k_str)
    #     return

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


    def bit_repr(self,index=None):
        if not isinstance(index,np.ndarray):
            index = np.array(range(self.nstates),dtype=int)
        string ='-'*(3*(self.nb_sites+4))+'\n'
        string +='empty        : 0\nuparrow      : 1\ndownarrow    :-1\nupdoawnarrow : 2\n'
        string +='-'*(3*(self.nb_sites+4))+'\n'
        string +='index | bit\n'
        string +='------|--'+'-'*(3*(self.nb_sites+1))+'\n'
        for index in self.basis[index]:
            if self.order=='spin':
                up,down = index & 2**self.nb_sites-1,index >> self.nb_sites
            else:
                up,down=0,0
                for k in range(0,2*self.nb_sites,2):
                    if (index&2**k == 2**k):
                        up += 2**(k//2)
                    if (index&2**(k+1)==2**(k+1)):
                        down +=2**(k//2)
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


    def exec_basis(self,elem2comp):
        with h5py.File(self.filename_basis,'a') as file:
            file.attrs['lenelem2comp'] = len(elem2comp)
            file.attrs['elem2comp'] = elem2comp
        os.system(self.exec)

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
        nup,ndown = hilbert
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
    def file_size(self):
        return float(os.path.getsize(self.filename_basis))

    @property 
    def basis(self):
        with h5py.File(self.filename_basis,  "r") as file:
            if self.order == 'spin':
                return np.concatenate(([[np.sum(self.spin_ordering(i,j)) \
                    for i,j in product(np.array(file[f'{int2str(self.nb_sites,2)}/{int2str(u,2)}/basis'],dtype=int),np.array(file[f'{int2str(self.nb_sites,2)}/{int2str(d,2)}/basis'],dtype=int))] \
                    for u,d in product(self.nup,self.ndown)]))
            elif self.order == 'site':
                return np.concatenate(([[np.sum(self.site_ordering(i,j)) \
                    for i,j in product(np.array(file[f'{int2str(self.nb_sites,2)}/{int2str(u,2)}/basis'],dtype=int),np.array(file[f'{int2str(self.nb_sites,2)}/{int2str(d,2)}/basis'],dtype=int))] \
                    for u,d in product(self.nup,self.ndown)]))
    @property 
    def basis_up(self):
        with h5py.File(self.filename_basis,  "r") as file:
            return np.concatenate(([np.array(file[f'{int2str(self.nb_sites,2)}/{int2str(u,2)}/basis'],dtype=int) for u in self.nup]))
    @property 
    def basis_down(self):
        with h5py.File(self.filename_basis,  "r") as file:
            return np.concatenate(([np.array(file[f'{int2str(self.nb_sites,2)}/{int2str(d,2)}/basis'],dtype=int) for d in self.ndown]))


    @property
    def nb_states(self):
        with h5py.File(self.filename_basis,  "r") as file:
            return self.nb_sup * self.nb_sdown
    @property
    def nb_sup(self):
        with h5py.File(self.filename_basis,  "r") as file:
            return np.sum([np.array(file[f'{int2str(self.nb_sites,2)}/{int2str(u,2)}'].attrs['nstates'][0],dtype=int) for u in self.nup])
    @property
    def nb_sdown(self):
        with h5py.File(self.filename_basis,  "r") as file:
            return np.sum([np.array(file[f'{int2str(self.nb_sites,2)}/{int2str(d,2)}'].attrs['nstates'][0],dtype=int) for d in self.ndown])


    @property
    def nstates(self):
        return self.nb_states
    @property
    def nsup(self):
        return self.nb_sup
    @property
    def nsdown(self):
        return self.nb_sdown

    @property
    def n_up(self):
        return self.nup

    @property
    def n_down(self):
        return self.ndown
