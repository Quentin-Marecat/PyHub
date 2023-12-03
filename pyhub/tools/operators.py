import numpy as np 
from pyhub.core._operators import Operators


class c_dagger_c(Operators):

    def __init__(self,*args):
        Operators.__init__(self,*args)
        self.str_index = [['c_dagger_c']]

class c_dagger(Operators):

    def __init__(self,*args):
        Operators.__init__(self,*args)
        self.str_index = [['c_dagger']]

class _c(Operators):

    def __init__(self,*args):
        Operators.__init__(self,*args)
        self.str_index = [['c']]

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

# def sz(*args):
#     return (n((args[0][0],'u')) - n((args[0][0],'d'))) / 2.

# class sisj(Operators):

#     def __init__(self,*args):
#         Operators.__init__(self,*args)
#         self.str_index = [['sisj']]

def sisj(*args):
    op = empty_operator()
    if args[0][0]==args[1][0]:
        op += n((args[0][0],'up'))*_n((args[0][0],'down'))/2 + n((args[0][0],'up'))*_n((args[0][0],'down'))/2
    else:
        op -= ( c_dagger_c((args[0][0],'u'),(args[1][0],'u')) * c_dagger_c((args[1][0],'d'),(args[0][0],'d'))\
                 + c_dagger_c((args[0][0],'d'),(args[1][0],'d')) * c_dagger_c((args[1][0],'u'),(args[0][0],'u')) )/2
#    op += ((n((args[0][0],'u')) - n((args[0][0],'d'))) / 2.) * ((n((args[1][0],'u')) - n((args[1][0],'d'))) / 2.)
    op += sz((args[0][0],))*sz((args[1][0],))
    return op

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
    if isinstance(psi,(list,np.ndarray)):
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
    else:
        matrix = np.zeros((op.nb_selected,op.nb_selected))
        psi0 = np.zeros(op.nb_selected)
        for i in range(op.nb_selected):
            psi0[i] = 1.
            matrix[i,:] = opeexp(op,psi0)
            psi0[i] = 0.
            
        return matrix



def opesum(lst):
    H = empty_operator()
    for op in lst:
        H += op 
    return H




# if __name__ == '__main__':
#     from pyhub.core.basis import Basis
#     op1 = c_dagger_c((0,"u"),(1,"u"))
#     op2 = n((0,"u")) + 3*n((1,"u"))
#     op3 = 2.*op1+op2+op2 - op1
#     op4 = 2.*idt()
#     op5 = op1 + 2*idt()
#     op6 = op2*op5 + 2.
#     op7 = sz((0,))
#     op8 = sisj((0,),(1,))
#     print(op8)
#     # print(op2)
#     # print(op5)
#     # print(op3)
#     # print(op4)
#     # print(op1+op4)
#     # print(op1+2)
#     # print()
#     # print('op2\n',op2)
#     # print('op5\n',op5)
#     # print(op4*op2)
#     # print('op2*op5\n',op6)
#     # for op in op6:
#     #     print(op)
#     mbbasis = Basis(4,hilbert=(2,2))
#     print(op6)
#     op6.set_basis(mbbasis)
# #    op6([1,2,3])