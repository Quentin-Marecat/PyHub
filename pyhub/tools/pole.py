#!/usr/bin/env python3
import numpy as np
from typing import Callable

class Pole():
    r"""
    This class allow manipulation of pole
    A pole is defined by 2 parameters within a dict :
        'positions'->list
        'weights'->list
    each poles are ordered within a list
    """
    def __init__(self,pole:dict,atol=1.e-10):
        if list(pole.keys()) != (['positions','weights'] or ['weights','positions']):
            raise ValueError(f"Pole asked is not correct, use correct notation -> {str({'positions':[],'weights':[]})} not {pole.keys()}")
        self.pole={'positions':np.array(pole['positions'],copy=True,dtype=float),'weights':np.array(pole['weights'],copy=True,dtype=float)}
        self.atol=atol
        if len(pole['weights']) != len(pole['positions']):
            raise ValueError(f"number of poles {len(pole['positions'])} should be equal as len of weights {len(pole['weights'])}")

    def __call__(self,w,eta:float=0.05):
        return self.retarded(w,eta=eta)

    def __str__(self):
        return f"{self.pole}"


    def __add__(self,other):
        return Pole({'positions':np.concatenate((self.pole['positions'],other.pole['positions'])),\
            'weights':np.concatenate((self.pole['weights'],other.pole['weights']))})


    def __sub__(self,other):
        return Pole({'positions':np.concatenate((self.pole['positions'],other.pole['positions'])),\
            'weights':np.concatenate((self.pole['weights'],-other.pole['weights']))})


    def __mul__(self,other):
        if isinstance(other,(float,int,complex,np.complex128)):
            return Pole({'positions':self.pole['positions'],'weights':self.pole['weights']*other})
        else:
            raise TypeError(f"Other must be float complex or int, not {type(other)}")

    def __truediv__(self,other):
        if isinstance(other,(float,int)):
            return Pole({'positions':self.pole['positions'],'weights':self.pole['weights']/other})
        else:
            raise TypeError(f"Other must be float or int, not {type(other)}")


    def shift(self,s:float):
        return Pole({'positions':s+self.pole['positions'],'weights':self.pole['weights']})

    def raise_deg(self,value:float):
        for k in range(2):
            index=np.where(self.atol>=abs(self.pole[k]['positions']-value))[0]
            if len(index)>0 :
                s=np.sum(self.pole[k]['weights'][index])
                for i in range(len(index)):
                    self.pole[k]['positions']=np.delete(self.pole[k]['positions'],index[0])
                    self.pole[k]['weights']=np.delete(self.pole[k]['weights'],index[0])
                self.pole[k]['positions']=np.insert(self.pole[k]['positions'],index[0],value+self.atol*(np.random.random()+1)*1.e2)
                self.pole[k]['weights']=np.insert(self.pole[k]['weights'],index[0],s/(len(index)+1))
                for i in range(1,len(index)+1):
                    self.pole[k]['positions']=np.insert(self.pole[k]['positions'],index[0]+(-1)**i*(i//2+1)+i%2,value+(-1)**i*(i//2+1)*self.atol*(np.random.random()+1)*1.e2)
                    self.pole[k]['weights']=np.insert(self.pole[k]['weights'],index[0]+(-1)**i*(i//2+1)+i%2,s/(len(index)+1))


    def clean(self,ptol=1.e-5,wtol=1.e-5):
        index = np.where(abs(self.pole['weights'])<wtol)[0]
        self.pole['weights'] = np.delete(self.pole['weights'],index)
        self.pole['positions'] = np.delete(self.pole['positions'],index)
        self.pole['positions'],index = np.unique(np.around(self.pole['positions'],int(-np.log10(ptol))),return_inverse=True)
        new_weights = np.zeros(self.nb_poles)
        for j,w in enumerate(index):
            new_weights[w] += self.pole['weights'][j]
        self.pole['weights'] = new_weights
        return self


    def retarded(self,w,eta:float=0.05)->complex:
        r"""
        This function compute the retarded GF where poles are given in arguments.
        """
        if isinstance(w,(float,complex)):
            return sum([self.pole['weights'][i]/(w-self.pole['positions'][i]+1j*eta) for i in range(self.nb_poles)])
        elif isinstance(w,(list,np.ndarray)):
            return np.array([sum([self.pole['weights'][i]/(w_-self.pole['positions'][i]+1j*eta) for i in range(self.nb_poles)]) for w_ in w])

    def advanced(self,w,eta:float=0.05)->complex:
        r"""
        This function compute the advanced GF where poles are given in arguments.
        """
        if isinstance(w,(float,complex)):
            return sum([self.pole['weights'][i]/(w-self.pole['positions'][i]-1j*eta) for i in range(self.nb_poles)])
        elif isinstance(w,(list,np.ndarray)):
            return np.array([sum([self.pole['weights'][i]/(w_-self.pole['positions'][i]-1j*eta) for i in range(self.nb_poles)]) for w_ in w])

    def greater(self,w,ef,eta:float=0.05)->complex:
        r"""
        This function compute the greater GF where poles are given in arguments.
        """
        if isinstance(w,(float,complex)):
            return sum([self.pole['weights'][i]/(w-p-1j*eta) if p>ef else 0. for i,p in enumerate(self.pole['positions'])])
        elif isinstance(w,(list,np.ndarray)):
            return np.array([sum([self.pole['weights'][i]/(w-p-1j*eta) if p>ef else 0. for i,p in enumerate(self.pole['positions'])]) for w_ in w])

    def lesser(self,w,ef,eta:float=0.05)->complex:
        r"""
        This function compute the lesser GF where poles are given in arguments.
        """
        if isinstance(w,(float,complex)):
            return sum([self.pole['weights'][i]/(w-p+1j*eta) if p<ef else 0. for i,p in enumerate(self.pole['positions'])])
        elif isinstance(w,(list,np.ndarray)):
            return np.array([sum([self.pole['weights'][i]/(w-p+1j*eta) if p<ef else 0. for i,p in enumerate(self.pole['positions'])]) for w_ in w])

    def causal(self,w,ef,eta:float=0.05)->complex:
        if isinstance(w,(float,complex)):
            return sum([self.pole['weights'][i]/(w-self.pole['positions'][i]+np.sign(self.pole['positions'][i]-ef)*1j*eta) for i in range(self.nb_poles)])
        elif isinstance(w,(list,np.ndarray)):
            return np.array([sum([self.pole['weights'][i]/(w_-self.pole['positions'][i]+np.sign(self.pole['positions'][i]-ef)*1j*eta) for i in range(self.nb_poles)]) for w_ in w]) 


    def spectral_density(self,w,eta:float=0.05):
        return -1.*np.imag(self.retarded(w,eta))/np.pi


    @classmethod
    def empty(cls):
        return cls({'positions':np.array([]),'weights':np.array([])})

    @classmethod
    def zero(cls):
        return cls({'positions':np.array([0.]),'weights':np.array([0.])})


    @classmethod
    def pole_research(cls,gf:Callable[[float],complex],borne:np.ndarray,eta:float=0.01,crit:float=1.,eta_is_var:bool=True,constant=True):
        from scipy.optimize import minimize_scalar as ms
        r"""
        This function research Lorentzian parameters of the gf
        """
        low,sup,nb_steps,dw=borne[0],borne[-1],len(borne),np.average(borne[1:]-borne[:-1])
        pole={'positions':np.array([]),'weights':np.array([])}
        elem_list=[-gf(low-dw,eta=eta).imag,-gf(low,eta=eta).imag]
        new_eta=eta/10 if eta_is_var else eta
        for k in range(nb_steps) :
            elem_list.append(-gf(low + (k+1)*dw,eta=eta).imag)
            if abs(elem_list[1]) >= crit*abs(elem_list[0]) and\
                abs(elem_list[1]) >= crit*abs(elem_list[2]):
                w_max=ms(lambda w:-abs(gf(w,eta=new_eta).imag),bounds=(low + (k-1)*dw,low + (k+1)*dw), method='bounded').x
                pole['positions']=np.append(pole['positions'],w_max)
                pole['weights']=np.append(pole['weights'],-gf(w_max,eta=new_eta).imag*(new_eta))
            del elem_list[0]
        if constant:
            if abs(gf(1.e5,eta=eta).real) >=1/1.e4:
                pole['positions']=np.append(pole['positions'],1.e8)
                pole['weights']=np.append(pole['weights'],-(1.e8)*gf(1.e5,eta=eta).real)
        return cls(pole)

    @staticmethod
    def sum(l:np.ndarray):
        pole=Pole.empty()
        for p in l:
            pole=pole+p 
        return pole

    @property 
    def nb_poles(self):
        return len(self.pole['positions'])



if __name__=='__main__':
    pole1=Pole({'positions':[1,2,4],'weights':[0,1,2]})
    pole2=Pole({'positions':[3,0],'weights':[7,2]})
    print(pole1);print(pole2)
    print(pole1 + pole2)
    print((pole1 + pole2)/2.)
    print(Pole.shift(pole1,2.))
    print(pole1-pole2)
    print(pole1*2)

    np.random.seed(10)
    positions,weights = np.random.random(20),10*np.random.random(20)
    gf = lambda w,eta=0.01: np.sum([weights[i]/(w-positions[i]+1j*eta) for i in range(20)])
    wgrid = np.linspace(-1,2,2000)
    pole3 = Pole.empty()
    pole3 = pole3.pole_research(gf,wgrid) 
    pole4 = Pole({'positions':positions,'weights':weights})

    import matplotlib.pyplot as plt
    plt.plot(wgrid,pole4.spectral_density(wgrid,eta=0.02),color='blue',label='exact')
    plt.plot(wgrid,pole3.spectral_density(wgrid,eta=0.02),color='red',linestyle='--',label='fit')
    plt.legend()
    plt.show()
#    plt.savefig('test')
