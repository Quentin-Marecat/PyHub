import numpy as np 
import re
from itertools import product
import os
import h5py
from pyhub.tools.tools import find_file
from pyhub.core.basis import Basis
from scipy.linalg import eigh_tridiagonal

class Operators(Basis):

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
            # index = np.where(np.abs(add.coeff)<1.e-14)[0]
            # k=0
            # for i in index:
            #     del add.coeff[i-k]
            #     del add.elem_index[i-k]
            #     del add.str_index[i-k]
            #     k+=1
            #     if (len(add.elem_index)==1):
            #         break
        elif isinstance(other,(float,int)):
            add = self + other*Operators._idt()
            return add
        return add

    def __radd__(self,other):
        return self.__add__(other)

    def __sub__(self,other):
        sub = Operators.empty_operator()
        if isinstance(other,Operators):
            sub.coeff = self.coeff + [-1.*c for c in other.coeff]
            sub.elem_index = self.elem_index + other.elem_index
            sub.str_index = self.str_index + other.str_index
            # index = np.where(np.abs(sub.coeff)<1.e-14)[0]
            # k=0
            # for i in index:
            #     del sub.coeff[i-k]
            #     del sub.elem_index[i-k]
            #     del sub.str_index[i-k]
            #     k+=1
            #     if (len(sub.elem_index)==1):
            #         break
        elif isinstance(other,(float,int)):
            sub = self - other*Operators._idt()
            return sub
        return sub

    def __rsub__(self,other):
        return -1*self.__sub__(other)

    def __mul__(self,other):
        mul = Operators.empty_operator()
        if isinstance(other,Operators):
            mul.coeff = [c1*c2 for c1,c2 in product(self.coeff,other.coeff)]
            mul.elem_index = [[*c1,*c2] for c1,c2 in product(self.elem_index,other.elem_index)]
            mul.str_index =  [[*c1,*c2] for c1,c2 in product(self.str_index,other.str_index)]
        elif isinstance(other,(float,int)):
            mul.coeff = [other*c for c in self.coeff]
            mul.elem_index = self.elem_index
            mul.str_index = self.str_index
            # index = np.where(np.abs(mul.coeff)<1.e-14)[0]
            # k=0
            # for i in index:
            #     del mul.coeff[i-k]
            #     del mul.elem_index[i-k]
            #     del mul.str_index[i-k]
            #     k+=1
            #     if (len(mul.elem_index)==1):
            #         break
            if not self.op2write:
                mul.set_basis(Basis(self.nb_sites,self.hilbert,self.order))
                mul.selected = self.selected
                mul.nb_selected = self.nb_selected
                mul._write_operator()
        return mul

    def __neg__(self):
        neg = Operators.empty_operator()
        neg.coeff = [-1.*c for c in self.coeff]
        neg.elem_index = self.elem_index
        neg.str_index = self.str_index
        if not self.op2write:
            neg.set_basis(Basis(self.nb_sites,self.hilbert,self.order))
            neg.selected = self.selected
            neg.nb_selected = self.nb_selected
            neg._write_operator()
        return neg

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

    def __getitem__(self, index):
        if isinstance(index, (list,np.ndarray)):
            if isinstance(index,np.ndarray):
                index = index.astype(int).tolist()
            other = Operators.empty_operator()
            other.set_basis(Basis(self.nb_sites,self.hilbert,self.order))
            other.coeff = self.coeff
            other.elem_index = self.coeff
            other.str_index = self.str_index
            other.op2write = False
            other.exec = self.exec
            other._index = self._index
            other.selected = index 
            other.nb_selected = len(index)
            other.max_len_ope = self.max_len_ope
            return other

    def _set_basis_elem(self,posu,posd,up,down):
        with h5py.File(self.filename_operators,'a') as f:
            f[f'{self.index_str}/psi'].attrs['posu'] = posu+1
            f[f'{self.index_str}/psi'].attrs['posd'] = posd+1
            f[f'{self.index_str}/psi'].attrs['up'] = up
            f[f'{self.index_str}/psi'].attrs['down'] = down

    def __req__(self,other):
        return self.__eq__(other)


    def __call__(self,psi,avg=False):
        if self.op2write:
            raise AttributeError(f'Set the basis set using self.set_basis(YourBasis:Basis)')
        if isinstance(psi[0],np.complex128):
            if not avg:
                return self.__call__(psi.real,avg=avg) + 1j*self.__call__(psi.imag,avg=avg)
            else:
                return self.__call__(psi.real,avg=avg) + self.__call__(psi.imag,avg=avg)
        with h5py.File(self.filename_operators,'a') as f:
            f[f'{self.index_str}/psi'].attrs['avg'] = avg
            f[f'{self.index_str}/psi'].attrs['avg_value'] = 0.
            if not isinstance(self.selected,(list,np.ndarray)):
                f[f'{self.index_str}/psi/input'][:] = psi
            else:
                print(psi)
                f[f'{self.index_str}/psi/input'][:] = np.zeros(self.nstates)
                f[f'{self.index_str}/psi/input'][self.selected] = psi
        if self.stable :
            self.exec_operators()
            with h5py.File(self.filename_operators,'a') as f:
                if not avg:
                    return f[f'{self.index_str}/psi']['output'][:] if self.selected is None else f[f'{self.index_str}/psi']['output'][self.selected] 
                else:
                    return f[f'{self.index_str}/psi'].attrs['avg_value']
        else:
            with h5py.File(self.filename_operators,'a') as f:
                f[f'{self.index_str}'].attrs['basis_elem'] = True
                f[f'{self.index_str}/psi'].attrs['avg'] = False
            if isinstance(self.selected,(list,np.ndarray)):
                psi_out = np.zeros(self.nb_selected)
                for i,index in enumerate(self.basis[self.selected]):
                    up,down = self.decompose(index)
                    posu = np.where(up==self.basis_up)[0][0]
                    posd = np.where(down==self.basis_up)[0][0]
                    self._set_basis_elem(posu,posd,up,down)
                    self.exec_operators()
                    with h5py.File(self.filename_operators,'a') as f:
                        val = f[f'{self.index_str}/psi']['input'][i]
                        psi_out += val*f[f'{self.index_str}/psi']['output'][self.selected]
                if avg:
                    with h5py.File(self.filename_operators,'a') as f:
                        f[f'{self.index_str}/psi'].attrs['avg_value'] = f[f'{self.index_str}/psi']['input'][:]@psi_out
            else:
                psi_out = np.zeros(self.nstates)
                for posd,down in enumerate(self.basis_down):
                    for posu,up in enumerate(self.basis_up):
                        self._set_basis_elem(posu,posd,up,down)
                        self.exec_operators()
                        with h5py.File(self.filename_operators,'a') as f:
                            val = f[f'{self.index_str}/psi']['input'][posd*self.nsup+posu]
                            psi_out += val*f[f'{self.index_str}/psi']['output'][:]
                if avg:
                    with h5py.File(self.filename_operators,'a') as f:
                        f[f'{self.index_str}/psi'].attrs['avg_value'] = f[f'{self.index_str}/psi']['input'][:]@psi_out
            with h5py.File(self.filename_operators,'a') as f:
                f[f'{self.index_str}'].attrs['basis_elem'] = False
                f[f'{self.index_str}/psi'].attrs['avg'] = avg
            if not avg:
                return psi_out
            else :
                with h5py.File(self.filename_operators,'a') as f:
                    return f[f'{self.index_str}/psi'].attrs['avg_value']

    def exec_operators(self):
        if self.op2write:
            raise AttributeError(f'Set the basis set using self.set_basis(YourBasis:Basis)')
        with h5py.File(self.filename_basis,'a') as file:
            file.attrs['lenelem2comp'] = len(self.n_up)+len(self.n_down)
            file.attrs['elem2comp'] = [*self.n_up,*self.n_down]
        with h5py.File(self.filename_operators,'a') as file:
            file.attrs['size_index']=len(str(self._index))
            file.attrs['index']=self._index
            if 'output' in file[f'{self.index_str}/psi'].keys():
                file[f'{self.index_str}/psi'].__delitem__('output')
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
        Basis.__init__(self,B.nb_sites,hilbert=B.hilbert,order=B.order)
        self.nb_selected = len(selected) if isinstance(selected,np.ndarray) else self.nstates
        self.selected = selected if isinstance(selected,np.ndarray) else None
        self._set_basis = True
        if self.op2write:
            self._write_operator(max_len_ope=40)

    def avg(self,psi):
        return self.__call__(psi,avg=True)
            

    def __delitem__(self):
        with h5py.File(self.filename_operators,"a") as file:
            file.__delitem__(self.index_str)


    @property
    def index_str(self):
        index_str = str(self._index)
        while len(index_str)<6:
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
        if self.op2write:
            raise AttributeError(f'Set the basis set using self.set_basis(YourBasis:Basis)')
        with h5py.File(self.filename_operators,'a') as f:
            f[f'{self.index_str}'].attrs['basis_elem'] = True
            f[f'{self.index_str}/psi'].attrs['avg'] = False
        if isinstance(self.selected,(list,np.ndarray)):
            H = np.zeros((self.nb_selected,self.nb_selected))
            for i,index in enumerate(self.basis[self.selected]):
                up,down = self.decompose(index)
                posu = np.where(up==self.basis_up)[0][0]
                posd = np.where(down==self.basis_up)[0][0]
                self._set_basis_elem(posu,posd,up,down)
                self.exec_operators()
                with h5py.File(self.filename_operators,'a') as f:
                    H[i,:] = f[f'{self.index_str}/psi']['output'][self.selected]
        else:
            H = np.zeros((self.nstates,self.nstates))
            for posd,down in enumerate(self.basis_down):
                for posu,up in enumerate(self.basis_up):
                    self._set_basis_elem(posu,posd,up,down)
                    self.exec_operators()
                    with h5py.File(self.filename_operators,'a') as f:
                        H[posd*self.nsup+posu,:] = f[f'{self.index_str}/psi']['output'][:]
        with h5py.File(self.filename_operators,'a') as f:
            f[f'{self.index_str}'].attrs['basis_elem'] = False
        return H
#        return np.array([self@np.array([1. if j==i else 0. for j in range(self.nb_selected)]) for i in range(self.nb_selected)])

    def lanczos(self,maxstep=200,acc_lcz=1.e-8,neigen=1,init_states = None):
        if self.op2write:
            raise AttributeError(f'Set the basis set using self.set_basis(YourBasis:Basis)')
        with h5py.File('operators.h5','a') as f:
            if 'lanczos' in f[f'{self.index_str}'].keys():
                f[f'{self.index_str}'].__delitem__('lanczos')
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
#            matrix = np.diag(alpha) + np.diag(beta[1:nstep],k=1) + np.diag(beta[1:nstep],k=-1)
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
        self.stable = True
        self.max_len_ope = max_len_ope
        self.exec = find_file('../../','operators.x')
        self.path = self.exec.replace('operators.x', '')

        ## index attribtion
        if os.path.exists(self.filename_operators):
            with h5py.File(self.filename_operators,"a") as file:
                list_index = np.array([*file.keys()],dtype=int)
            self._index = list_index[-1]+1
        else:
            self._index = 1
        with h5py.File(self.filename_operators,'a') as f:
            grp = f.create_group(f'{self.index_str}')
            grp.attrs['basis_elem'] = False
            grp1 = grp.create_group(f'psi')
            grp1.attrs['avg_value'] = 0.
            grp1.attrs['nstates'] = self.nstates
            grp1.attrs['nsup'] = self.nsup
            grp1.attrs['nsdown'] = self.nsdown
            grp1.attrs['nelemup'] = len(self.n_up)
            grp1.attrs['nelemdown'] = len(self.n_down)
            grp1.create_dataset('input', data=np.zeros(self.nstates),dtype = np.float64)
            grp1.create_dataset('output', data=np.zeros(self.nstates),dtype = np.float64)
        spin1 ,site1 ,spin2 ,site2 ,string = [np.zeros((max_len_ope,self.nb_ope),dtype=int) for _ in range(5)]
        coeff,nprod = np.zeros(self.nb_ope,dtype=np.float64) ,np.zeros(self.nb_ope,dtype=int) 
        if len(self.coeff)>0:
            coeff = self.coeff
            for j,op in enumerate(self):
                nprod[j] = len(op.elem_index[0])
                if nprod[j]>0:
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
                        elif op.str_index[0][i] == 'c_dagger':
                            string[i][j] = 10
                        elif op.str_index[0][i] == 'c':
                            string[i][j] = 11
                        else :
                            raise NotImplementedError
        with h5py.File('operators.h5','a') as f:
            grp = f.create_group(f'{self.index_str}/operators')
            grp.create_dataset('coeff', data=coeff)
            grp.create_dataset('str_index', data=string)
            grp.create_dataset('site1_index', data=site1)
            grp.create_dataset('spin1_index', data=spin1)
            grp.create_dataset('site2_index', data=site2)
            grp.create_dataset('spin2_index', data=spin2)
            grp.attrs['nb_ope'] = self.nb_ope if self.nb_ope>1 else 1
            grp.attrs['max_nb_ope'] = max_len_ope


    def bit_repr(self):
        with h5py.File(self.filename_operators,'a') as f:
            f[f'{self.index_str}'].attrs['basis_elem'] = True
            f[f'{self.index_str}/psi'].attrs['avg'] = True
        select = self.selected if isinstance(self.selected,(list,np.ndarray)) else list(range(self.nstates))
        string ='-'*(3*(self.nb_sites+4))+'\n'
        string +='empty        : 0\nuparrow      : 1\ndownarrow    :-1\nupdoawnarrow : 2\n'
        string +='-'*(3*(self.nb_sites+4))+'\n'
        string +='index | bit'+' '*(3*(self.nb_sites-1)-1)+' | <H> '+'\n'
        string +='------|----'+'-'*(3*(self.nb_sites-1)-1)+'-|-----'+'\n'
        for index in self.basis[select]:
            up,down = self.decompose(index)
            posu = np.where(up==self.basis_up)[0][0]
            posd = np.where(down==self.basis_down)[0][0]
            self._set_basis_elem(posu,posd,up,down)
            self.exec_operators()
            with h5py.File(self.filename_operators,'a') as f:
                avg=f[f'{self.index_str}/psi'].attrs['avg_value']
            up,down = self.decompose(index)
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
            string+=f' | {np.around(avg,4)}'+'\n'
        with h5py.File(self.filename_operators,'a') as f:
            f[f'{self.index_str}'].attrs['basis_elem'] = False
            f[f'{self.index_str}/psi'].attrs['avg'] = False
        return string


    def find_stable_space(self):
        print('Experimental function, do not use!')
        val = 0
        with h5py.File(self.filename_operators,'a') as f:
            f[f'{self.index_str}'].attrs['basis_elem'] = True
            f[f'{self.index_str}/psi'].attrs['avg'] = False
        if isinstance(self.selected,(list,np.ndarray)):
            table = [self.basis[self.selected],[None for i in range(self.nb_selected)]]
            for i,index in enumerate(self.basis[self.selected]):
                if table[1][i] is None:
                    up,down = self.decompose(index)
                    posu = np.where(up==self.basis_up)[0]
                    posd = np.where(down==self.basis_up)[0]
                    self._set_basis_elem(posu,posd,up,down)
                    self.exec_operators()
                    with h5py.File(self.filename_operators,'a') as f:
                        vec = f[f'{self.index_str}/psi']['output'][self.selected]
                    pos = np.concatenate(([i],np.where(abs(vec)>1.e-14)[0]))
                    for j in pos:
                        if table[1][j] is None:
                            table[1][j] = val 
                    val+=1 
        else:
            table = [self.basis,[None for i in range(self.nb_selected)]]
            for posd,down in enumerate(self.basis_down):
                for posu,up in enumerate(self.basis_up):
                    if table[1][posd*self.nsup+posu] is None:
                        self._set_basis_elem(posu,posd,up,down)
                        self.exec_operators()
                        with h5py.File(self.filename_operators,'a') as f:
                            vec = f[f'{self.index_str}/psi']['output'][:]
                        pos = np.concatenate(([posd*self.nsup+posu],np.where(abs(vec)>1.e-14)[0]))
                        for j in pos:
                            if table[1][j] is None:
                                table[1][j] = val 
                        val+=1 
        table_ = []
        for i in range(val):
            index = np.where(np.array(table[1][:]) == i)[0]
            table_.append(index)
        with h5py.File(self.filename_operators,'a') as f:
            f[f'{self.index_str}'].attrs['basis_elem'] = False
        return table_


    def __eq__(self,other):
        if isinstance(other,(float,int)):
            avg = []
            with h5py.File(self.filename_operators,'a') as f:
                f[f'{self.index_str}'].attrs['basis_elem'] = True
                f[f'{self.index_str}/psi'].attrs['avg'] = True
            if isinstance(self.selected,(list,np.ndarray)):
                for i,index in enumerate(self.basis[self.selected]):
                    up,down = self.decompose(index)
                    posu = np.where(up==self.basis_up)[0][0]
                    posd = np.where(down==self.basis_down)[0][0]
                    self._set_basis_elem(posu,posd,up,down)
                    self.exec_operators()
                    with h5py.File(self.filename_operators,'a') as f:
                        avg.append(np.isclose(f[f'{self.index_str}/psi'].attrs['avg_value'],other))
                        f[f'{self.index_str}/psi'].attrs['avg_value'] = 0.
            else:
                for posu,up in enumerate(self.basis_up):
                    for posd,down in enumerate(self.basis_down):
                        self._set_basis_elem(posu,posd,up,down)
                        self.exec_operators()
                        with h5py.File(self.filename_operators,'a') as f:
                            avg.append(np.isclose(f[f'{self.index_str}/psi'].attrs['avg_value'],other))
                            f[f'{self.index_str}/psi'].attrs['avg_value'] = 0.
            with h5py.File(self.filename_operators,'a') as f:
                f[f'{self.index_str}'].attrs['basis_elem'] = False
                f[f'{self.index_str}/psi'].attrs['avg'] = False
        return avg

    def __ge__(self,other):
        if isinstance(other,(float,int)):
            avg = []
            with h5py.File(self.filename_operators,'a') as f:
                f[f'{self.index_str}'].attrs['basis_elem'] = True
                f[f'{self.index_str}/psi'].attrs['avg'] = True
            if isinstance(self.selected,(list,np.ndarray)):
                for i,index in enumerate(self.basis[self.selected]):
                    up,down = self.decompose(index)
                    posu = np.where(up==self.basis_up)[0][0]
                    posd = np.where(down==self.basis_down)[0][0]
                    self._set_basis_elem(posu,posd,up,down)
                    self.exec_operators()
                    with h5py.File(self.filename_operators,'a') as f:
                        avg.append(f[f'{self.index_str}/psi'].attrs['avg_value']>=other)
                        f[f'{self.index_str}/psi'].attrs['avg_value'] = 0.
            else:
                for posu,up in enumerate(self.basis_up):
                    for posd,down in enumerate(self.basis_down):
                        self._set_basis_elem(posu,posd,up,down)
                        self.exec_operators()
                        with h5py.File(self.filename_operators,'a') as f:
                            avg.append(f[f'{self.index_str}/psi'].attrs['avg_value']>=other)
                            f[f'{self.index_str}/psi'].attrs['avg_value'] = 0.
            with h5py.File(self.filename_operators,'a') as f:
                f[f'{self.index_str}'].attrs['basis_elem'] = False
                f[f'{self.index_str}/psi'].attrs['avg'] = False
        return avg

    def __gt__(self,other):
        if isinstance(other,(float,int)):
            avg = []
            with h5py.File(self.filename_operators,'a') as f:
                f[f'{self.index_str}'].attrs['basis_elem'] = True
                f[f'{self.index_str}/psi'].attrs['avg'] = True
            if isinstance(self.selected,(list,np.ndarray)):
                for i,index in enumerate(self.basis[self.selected]):
                    up,down = self.decompose(index)
                    posu = np.where(up==self.basis_up)[0][0]
                    posd = np.where(down==self.basis_down)[0][0]
                    self._set_basis_elem(posu,posd,up,down)
                    self.exec_operators()
                    with h5py.File(self.filename_operators,'a') as f:
                        avg.append(f[f'{self.index_str}/psi'].attrs['avg_value']>other)
                        f[f'{self.index_str}/psi'].attrs['avg_value'] = 0.
            else:
                for posu,up in enumerate(self.basis_up):
                    for posd,down in enumerate(self.basis_down):
                        self._set_basis_elem(posu,posd,up,down)
                        self.exec_operators()
                        with h5py.File(self.filename_operators,'a') as f:
                            avg.append(f[f'{self.index_str}/psi'].attrs['avg_value']>other)
                            f[f'{self.index_str}/psi'].attrs['avg_value'] = 0.
            with h5py.File(self.filename_operators,'a') as f:
                f[f'{self.index_str}'].attrs['basis_elem'] = False
                f[f'{self.index_str}/psi'].attrs['avg'] = False
        return avg


    def __le__(self,other):
        if isinstance(other,(float,int)):
            avg = []
            with h5py.File(self.filename_operators,'a') as f:
                f[f'{self.index_str}'].attrs['basis_elem'] = True
                f[f'{self.index_str}/psi'].attrs['avg'] = True
            if isinstance(self.selected,(list,np.ndarray)):
                for i,index in enumerate(self.basis[self.selected]):
                    up,down = self.decompose(index)
                    posu = np.where(up==self.basis_up)[0][0]
                    posd = np.where(down==self.basis_down)[0][0]
                    self._set_basis_elem(posu,posd,up,down)
                    self.exec_operators()
                    with h5py.File(self.filename_operators,'a') as f:
                        avg.append(f[f'{self.index_str}/psi'].attrs['avg_value']<=other)
                        f[f'{self.index_str}/psi'].attrs['avg_value'] = 0.
            else:
                for posu,up in enumerate(self.basis_up):
                    for posd,down in enumerate(self.basis_down):
                        self._set_basis_elem(posu,posd,up,down)
                        self.exec_operators()
                        with h5py.File(self.filename_operators,'a') as f:
                            avg.append(f[f'{self.index_str}/psi'].attrs['avg_value']<=other)
                            f[f'{self.index_str}/psi'].attrs['avg_value'] = 0.
            with h5py.File(self.filename_operators,'a') as f:
                f[f'{self.index_str}'].attrs['basis_elem'] = False
                f[f'{self.index_str}/psi'].attrs['avg'] = False
        return avg


    def __lt__(self,other):
        if isinstance(other,(float,int)):
            avg = []
            with h5py.File(self.filename_operators,'a') as f:
                f[f'{self.index_str}'].attrs['basis_elem'] = True
                f[f'{self.index_str}/psi'].attrs['avg'] = True
            if isinstance(self.selected,(list,np.ndarray)):
                for i,index in enumerate(self.basis[self.selected]):
                    up,down = self.decompose(index)
                    posu = np.where(up==self.basis_up)[0][0]
                    posd = np.where(down==self.basis_down)[0][0]
                    self._set_basis_elem(posu,posd,up,down)
                    self.exec_operators()
                    with h5py.File(self.filename_operators,'a') as f:
                        avg.append(f[f'{self.index_str}/psi'].attrs['avg_value']<other)
                        f[f'{self.index_str}/psi'].attrs['avg_value'] = 0.
            else:
                for posu,up in enumerate(self.basis_up):
                    for posd,down in enumerate(self.basis_down):
                        self._set_basis_elem(posu,posd,up,down)
                        self.exec_operators()
                        with h5py.File(self.filename_operators,'a') as f:
                            avg.append(f[f'{self.index_str}/psi'].attrs['avg_value']<other)
                            f[f'{self.index_str}/psi'].attrs['avg_value'] = 0.
            with h5py.File(self.filename_operators,'a') as f:
                f[f'{self.index_str}'].attrs['basis_elem'] = False
                f[f'{self.index_str}/psi'].attrs['avg'] = False
        return avg

    
    def __ne__(self,other):
        eq = self.__eq__(other)
        return [not e for e in eq]


def _operator_from_string(string,args):
    op = Operators(*args)
    if string == 'c_dagger_c':
        op.str_index = [['c_dagger_c']]
    elif string == 'c_dagger':
        op.str_index = [['c_dagger']]
    elif string == 'c':
        op.str_index = [['c']]
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
    elif string == 'c_dagger':
        return f'c^\dagger_{{{args[0][0]},{args[0][1]}}}'
    elif string == 'c':
        return f'c_{{{args[0][0]},{args[0][1]}}}'
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


        