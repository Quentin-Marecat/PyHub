import numpy as np 
from pyhub.core._operators import Operators


class c_dagger_c(Operators):

    def __init__(self,*args):
        Operators.__init__(self,*args)
        self.str_index = [['c_dagger_c']]

class n(Operators):

    def __init__(self,*args):
        Operators.__init__(self,*args)
        self.str_index = [['n']]

class _n(Operators):

    def __init__(self,*args):
        Operators.__init__(self,*args)
        self.str_index = [['1_n']]


class sz(Operators):

    def __init__(self,*args):
        Operators.__init__(self,*args)
        self.str_index = [['sz']]

class sisj(Operators):

    def __init__(self,*args):
        Operators.__init__(self,*args)
        self.str_index = [['sisj']]

class splus(Operators):

    def __init__(self,*args):
        Operators.__init__(self,*args)
        self.str_index = [['splus']]

class sminus(Operators):

    def __init__(self,*args):
        Operators.__init__(self,*args)
        self.str_index = [['sminus']]

class idt(Operators):

    def __init__(self):
        Operators.__init__(self)
        self.coeff = [1.]
        self.elem_index = [[()]]
        self.str_index = [['id.']]


class empty_operator(Operators):

    def __init__(self):
        Operators.__init__(self)
        self.coeff = []
        self.elem_index = []
        self.str_index = []

def opeexp(op:Operators,psi,ftol=1.e-6,unitary=False):
    conv=False
    new_psi,full_psi,i,val = np.copy(psi),np.copy(psi),1,1.e3
    nmr0 = np.linalg.norm(psi)
    while not conv and i<100:
        new_psi = op@new_psi/i
        full_psi += new_psi
        i+=1
        nmr = full_psi.conj()@full_psi
        if unitary:
            if np.abs(nmr0-nmr)<ftol:
                conv = True
        else:
            if np.abs(val-nmr)<ftol:
                conv = True
            else:
                val = nmr
    if i==100:
        raise ValueError('exponential operator does not converge')
    return full_psi


def opesum(lst):
    H = empty_operator()
    for op in lst:
        H += op 
    return H




if __name__ == '__main__':
    op1 = c_dagger_c((0,"u"),(1,"u"))
    op2 = n((0,"u")) + 3*n((1,"u"))
    op3 = 2.*op1+op2+op2 - op1
    op4 = 2.*idt()
    op5 = op1 + 2*idt()
    op6 = op2*op5 + 2.
    op7 = sz((0,))
    op8 = sisj((0,),(1,))
    print(op8)
    # print(op2)
    # print(op5)
    # print(op3)
    # print(op4)
    # print(op1+op4)
    # print(op1+2)
    # print()
    # print('op2\n',op2)
    # print('op5\n',op5)
    # print(op4*op2)
    # print('op2*op5\n',op6)
    # for op in op6:
    #     print(op)
    print(op6)
    print(op6.is_spin_restricted)
    op6([1,2,3])