import numpy as np 
import re
from itertools import product
import os
import h5py
from pyhub.tools.tools import find_file
from pyhub.core.basis import Basis
from scipy.linalg import eigh_tridiagonal
from scipy.sparse import csc_matrix 

class Operators():

    def __init__(self,*args):
        self.coeff = [1.]
        self.elem_index = [[(args)]]
        self.op2write = True

    def __add__(self,other):
        add = Operators.empty_operator()
        if isinstance(other,Operators):
            add.coeff = self.coeff + other.coeff
            add.elem_index = self.elem_index + other.elem_index
            add.str_index = self.str_index + other.str_index
            self.op2write = True
        elif isinstance(other,(float,int)):
            return self + other*Operators._idt()
        return add

    def __radd__(self,other):
        return self.__add__(other)

    def __sub__(self,other):
        sub = Operators.empty_operator()
        sub.coeff = self.coeff + [-1.*c for c in other.coeff]
        sub.elem_index = self.elem_index + other.elem_index
        sub.str_index = self.str_index + other.str_index
        self.op2write = True
        return sub

    def __rsub__(self,other):
        return self.__sub__(other)

    def __mul__(self,other):
        mul = Operators.empty_operator()
        self.op2write = True
        if isinstance(other,Operators):
            mul.coeff = [c1*c2 for c1,c2 in product(self.coeff,other.coeff)]
            mul.elem_index = [[*c1,*c2] for c1,c2 in product(self.elem_index,other.elem_index)]
            mul.str_index =  [[*c1,*c2] for c1,c2 in product(self.str_index,other.str_index)]
        elif isinstance(other,(float,int)):
            mul.coeff = [other*c for c in self.coeff]
            mul.elem_index = self.elem_index
            mul.str_index = self.str_index
        return mul

    def __rmul__(self,other):
        return self.__mul__(other)

    def __truediv__(self,other):
        if isinstance(other,(float,int)):
            return self.__mul__(1/float(other))

    def __rtruediv__(self,other):
        return self.__rdiv__(other)

    def __matmul__(self,other):
        return self.__call__(other)

    def __rmatmul__(self,other):
        return self.__call__(other.conj())

    def __eq__(self,other):
        if isinstance(other,(float,int)):
            bl = np.zeros(self.nb_selected)
            for i in range(self.nb_selected):
                psi = np.array([1. if j==i else 0. for j in range(self.nstates)])
                bl[i] = np.isclose(self.avg(psi),other)
        return bl 

    def __req__(self,other):
        return self.__eq__(other)


    def __call__(self,psi,max_len_ope=40,avg=False):
        if isinstance(psi[0],np.complex128):
            if not avg:
                return self.__call__(psi.real,max_len_ope=20,avg=avg) + 1j*self.__call__(psi.imag,max_len_ope=20,avg=avg)
            else:
                return self.__call__(psi.real,max_len_ope=20,avg=avg) + self.__call__(psi.imag,max_len_ope=20,avg=avg)
        if self.op2write:
            self._write_operator(max_len_ope=max_len_ope)
        with h5py.File(self.filename_operators,'a') as f:
            f[f'{self.index_str}/psi'].attrs['avg'] = avg
            if (self.selected is None) or (isinstance(self.selected,np.ndarray) and len(psi)==self.nstates):
                f[f'{self.index_str}/psi/input'][:] = psi
            else:
                f[f'{self.index_str}/psi/input'][self.selected] = psi
        self.exec_operators()
        with h5py.File(self.filename_operators,'a') as f:
            if not avg:
                return f[f'{self.index_str}/psi']['output'][:] if self.selected is None else f[f'{self.index_str}/psi']['output'][self.selected] 
            else:
                return f[f'{self.index_str}/psi'].attrs['avg_value']


    def exec_operators(self):
        with h5py.File(self.filename_basis,'a') as file:
            file.attrs['size_index']=len(str(self._basis_index))
            file.attrs['index']=self._basis_index
        with h5py.File(self.filename_operators,'a') as file:
            file.attrs['size_index']=len(str(self._index))
            file.attrs['index']=self._index
        os.system(self.exec)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        full_string = ''
        for i,coeff in enumerate(self.coeff):
            string = f'{np.around(coeff,3)}*'
            for j,elem in enumerate(self.elem_index[i]):
                string += _print_from_string(self.str_index[i][j],elem)
            string += ' + '
            full_string += string
        full_string = re.sub(r'\*\s',r' ',full_string)
        full_string = re.sub(r'(.*) \+ $',r'\1',full_string)
        full_string = re.sub(r'\+\s-',r'- ',full_string)
        return full_string

    def __iter__(self):
        iter_op = []
        for i,coeff in enumerate(self.coeff):
            iter_op.append(coeff*self.__multiplemul__([_operator_from_string(string,elem) for string,elem in zip(self.str_index[i],self.elem_index[i])]))
        return iter_op.__iter__()

    def __multiplemul__(self,other):
        mul = Operators._idt()
        for o in other : 
            mul *= o 
        return mul


    def set_basis(self,B:Basis,selected:np.ndarray=None,memspace='10Go'):
        self._basis_index = B.index
        self.nstates = B.nstates 
        self.filename_basis = B.filename_basis
        self.memspace = memspace
        self.nb_selected = len(selected) if isinstance(selected,np.ndarray) else self.nstates
        self.selected = selected
        self._set_basis = True

    def avg(self,psi):
        return self.__call__(psi,avg=True)
            

    def __delitem__(self,k):
        with h5py.File(self.filename_operators,"a") as file:
            k_str = str(k)
            while len(k_str)<3:
                k_str = '0'+k_str
            file.__delitem__(k_str)
        return


    @property
    def index_str(self):
        index_str = str(self._index)
        while len(index_str)<3:
            index_str = '0'+index_str
        return index_str


    @property
    def file_size(self):
        return float(os.path.getsize(self.filename_operators))

    @staticmethod
    def empty_operator():
        empty = Operators()
        empty.coeff = []
        empty.elem_index = []
        empty.str_index = []
        return empty

    @property
    def T(self):
        pass

    @property
    def dagger(self):
        return self.T 
    
    @property
    def to_matrix(self):
        return np.array([self@np.array([1. if j==i else 0. for j in range(self.nb_selected)]) for i in range(self.nb_selected)])

    @property
    def to_sparse_matrix(self):
        return csc_matrix(self.to_matrix,dtype=np.float64)

    def lanczos(self,maxstep=200,acc_lcz=1.e-8,neigen=1,init_states = None):
        if self.op2write:
            self._write_operator(max_len_ope=40)
        with h5py.File('operators.h5','a') as f:
            if 'lanczos' in f[f'{self.index_str}'].keys():
                raise NotImplementedError
            grp = f[f'{self.index_str}'].create_group('lanczos')
            grp.create_dataset('lcz_vectors', data=np.zeros((maxstep,self.nb_selected)))
        #v = np.full(self.nb_selected,1.)
        v = np.random.random(self.nb_selected)-.5  if not isinstance(init_states,np.ndarray) else init_states
        alpha,beta = [],[]
        beta.append(np.linalg.norm(v))
        v = v/beta[0]
        v0 = np.zeros(self.nb_selected)
        nstep,ek,ekm1 = 0,np.full(neigen,100.),np.full(neigen,0.)
        nadd = 5
        while nstep < maxstep-nadd and not np.linalg.norm(ek[:neigen]-ekm1)/neigen<acc_lcz:
            ekm1 = np.copy(ek[:neigen])
            for i in range(nadd):
                with h5py.File(self.filename_operators,'a') as f:
                    f[f'{self.index_str}/lanczos/lcz_vectors'][nstep,:] = v
                v2 = self@v 
                alpha.append(v2@v)
                v2 += -alpha[-1]*v - beta[-1]*v0
                beta.append(np.linalg.norm(v2))
                v0 = v
                v = v2/beta[-1] 
                nstep+=1
            matrix = np.diag(alpha) + np.diag(beta[1:nstep],k=1) + np.diag(beta[1:nstep],k=-1)
#            ek,Vk = np.linalg.eigh(matrix)
            ek,Vk = eigh_tridiagonal(alpha,beta[1:nstep],eigvals_only=False, select=f'a')
#            print(ek[:5])
#            print(alpha[-nadd:])
#            print(beta[-nadd:])
        with h5py.File(self.filename_operators,'a') as f:
            f[f'{self.index_str}/lanczos'].create_dataset('alpha', data=alpha)
            f[f'{self.index_str}/lanczos'].create_dataset('beta', data=alpha)
            f[f'{self.index_str}/lanczos'].create_dataset('ek', data=ek[:neigen])
            f[f'{self.index_str}/lanczos'].create_dataset('Vk', data=Vk[:,:neigen])
            psi = f[f'{self.index_str}/lanczos/lcz_vectors'][:len(alpha),:].T@Vk[:,:neigen]
            return ek[:neigen], psi

    @property 
    def alpha(self):
        with h5py.File(self.filename_operators,'a') as f:
            return f[f'{self.index_str}/lanczos/alpha'][:]

    @property 
    def beta(self):
        with h5py.File(self.filename_operators,'a') as f:
            return f[f'{self.index_str}/lanczos/beta'][:]



    @property
    def filename_operators(self):
        return 'operators.h5'

    @property 
    def nb_ope(self):
        return len(self.coeff)

    @staticmethod
    def _idt():
        idt = Operators()
        idt.coeff = [1.]
        idt.elem_index = [[()]]
        idt.str_index = [['id.']]
        return idt


    def _write_operator(self,max_len_ope=40):
        self.op2write = False
        try:
            self._set_basis
        except AttributeError:
            raise AttributeError(f'Set the basis set using self.set_basis(YourBasis:Basis)')
        self.exec = find_file('../','operators.x')
        self.path = self.exec.replace('operators.x', '')

        size,unit = float(re.sub(r'(\d+)\w?o?',r'\1',self.memspace)),re.sub(r'\d+(\w?o?)',r'\1',self.memspace)
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
            raise ValueError(f'Set memspace as a correct string, not {self.memspace}')
        ## index attribtion
        if os.path.exists(self.filename_operators):
            i=0
            with h5py.File(self.filename_operators,"a") as file:
                list_index = np.array([*file.keys()],dtype=int)
            while self.file_size > self.memspace and i <= len(list_index):
                self.__delitem__(list_index[i])
                i+=1
            index = np.sort(list_index)
            if len(index)>0:
                self._index = index[0]-1 if index[0]>1 else index[-1]+1
        else:
            self._index = 1
        with h5py.File(self.filename_operators,'a') as f:
            grp = f.create_group(f'{self.index_str}')
            grp1 = grp.create_group(f'psi')
            grp1.attrs['avg_value'] = 0.
            grp1.create_dataset('input', data=np.zeros(self.nstates),dtype = np.float64)
            grp1.create_dataset('output', data=np.zeros(self.nstates),dtype = np.float64)
        spin1 ,site1 ,spin2 ,site2 ,string = [np.zeros((max_len_ope,self.nb_ope),dtype=int) for _ in range(5)]
        nprod = np.zeros(self.nb_ope,dtype=int) 
        for j,op in enumerate(self):
            nprod[j] = len(op.elem_index[0])
            for i,index in enumerate(op.elem_index[0]):
                try:
                    site1[i,j] = index[0][0]+1
                    site2[i,j] = index[1][0]+1
                except :
                    pass
                try : 
                    spin1[i,j] = 1 if (index[0][1]==1 or index[0][1]=='u' or index[0][1]=='up' or index[0][1]=='+') else 2
                    spin2[i,j] = 1 if (index[1][1]==1 or index[1][1]=='u' or index[1][1]=='up' or index[1][1]=='+') else 2
                except :
                    pass
                if op.str_index[0][i] == 'c_dagger_c':
                    string[i][j] = 1 if index[0][1] == index[1][1] else 9
                elif op.str_index[0][i] == 'n':
                    string[i][j] = 2
                elif op.str_index[0][i] == 'id.':
                    string[i][j] = 3
                elif op.str_index[0][i] == '1_n':
                    string[i][j] = 4
                elif op.str_index[0][i] == 'sz':
                    string[i][j] = 5
                elif op.str_index[0][i] == 'splus':
                    string[i][j] = 6
                elif op.str_index[0][i] == 'sminus':
                    string[i][j] = 7
                elif op.str_index[0][i] == 'sisj':
                    string[i][j] = 8
                else :
                    raise NotImplementedError
        with h5py.File('operators.h5','a') as f:
            grp = f.create_group(f'{self.index_str}/operators')
            grp.create_dataset('coeff', data=self.coeff)
            grp.create_dataset('str_index', data=string)
            grp.create_dataset('site1_index', data=site1)
            grp.create_dataset('spin1_index', data=spin1)
            grp.create_dataset('site2_index', data=site2)
            grp.create_dataset('spin2_index', data=spin2)
            grp.attrs['nb_ope'] = self.nb_ope
            grp.attrs['max_nb_ope'] = max_len_ope
            try:
                grp.attrs['basis_index'] = self._basis_index
            except:
                raise AttributeError(f'Set the basis set using self.set_basis(YourBasis:Basis)')




def _operator_from_string(string,args):
    op = Operators(*args)
    if string == 'c_dagger_c':
        op.str_index = [['c_dagger_c']]
    elif string == 'n':
        op.str_index = [['n']]
    elif string == '1_n':
        op.str_index = [['1_n']]
    elif string == 'sz':
        op.str_index = [['sz']]
    elif string == 'sisj':
        op.str_index = [['sisj']]
    elif string == 'splus':
        op.str_index = [['splus']]
    elif string == 'sminus':
        op.str_index = [['sminus']]
    elif string == 'id.':
        op.str_index = [['id.']]
    else :
        raise ValueError(f'{string} is not an implemented operator')
    return op

def _print_from_string(string,args):
    if string == 'c_dagger_c':
        return f'c^\dagger_{{{args[0][0]},{args[0][1]}}}c_{{{args[1][0]},{args[1][1]}}}'
    elif string == 'n':
        return f'n_{{{args[0][0]},{args[0][1]}}}'
    elif string == '1_n':
        return f'(1 - n_{{{args[0][0]},{args[0][1]}}})'
    elif string == 'sz':
        return f'S^z_{{{args[0][0]}}}'
    elif string == 'sisj':
        return f'S_{{{args[0][0]}}}S_{{{args[1][0]}}}'
    elif string == 'splus':
        return f'S^+_{{{args[0][0]}}}'
    elif string == 'sminus':
        return f'S^-_{{{args[0][0]}}}'
    elif string == 'id.':
        return ''
    else :
        raise ValueError(f'{string} is not an implemented operator')