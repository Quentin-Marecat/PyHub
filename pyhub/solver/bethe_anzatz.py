#!/usr/bin/env python3
import numpy
import scipy.integrate as integrate
from scipy.special import jv as bessel
from scipy.special import zeta
from scipy.optimize import root_scalar
import warnings
warnings.filterwarnings("ignore",category = RuntimeWarning)

__author__ = "Quentin Marécat"
__maintainer__ = "Quentin Marécat"
__email__ = "quentin.marecat@gmail.com"
__date__ = "April, 2022"


def solve(u,q=numpy.pi, b=100, n_k=2000, n_lambda=2000, singlet=True):
    """
    Solves Bethe ansatz for given parameters Q, B, see Lieb, Wu, Physica A: Statistical Mechanics and its Applications,
    vol. 321, no. 1-2, pp. 1-27, Apr. 2003, section 4.

    Args:
        q (float): the wave number range;
        b (float): the lambda range;
        u (float): the Hubbard interaction term assuming the nearest heighbour hopping is equal to -1;
        n_k (int): number of points to discretize k-integrations;
        n_lambda (int): number of points to discretize lambda-integrations;

    Returns:
        A dict with all calculated quantities.
    """

    def k_func(x):
        """The K function, Eq. 23"""
        return (1/(2*numpy.pi))*(8 * u / (u ** 2 + 16 * x ** 2))

    def k2_func(x):
        """The K-squared function, Eq. 23. Note that is is not the same as the previous function squared."""
        return (1/(2*numpy.pi))*(4 * u / (u ** 2 + 4 * x ** 2))

    def R(x,nc=10):
        '''Eq 2.7 + truncated'''
        x2=numpy.power(x,2)
        x4=numpy.power(x,4)
        s=numpy.array([((-1)**(i+1)/(2*i))*(1/(1+(x2/(2*i)**2)) - 1 + (x2/(2*i)**2) - (x4/(2*i)**4)) for i in range(1,nc+1)])
        return (1/numpy.pi)*((0.5*numpy.log(2)-(3*x2/32)*zeta(3) + (15*x4/512)*zeta(5)) + numpy.sum(s,axis=0))

    # Discretize k 
    k = numpy.linspace(-q, q, n_k)
    dk = k[1]-k[0]
    k_i = k[:, numpy.newaxis]
    k_j = k[numpy.newaxis, :]

    # Discretize l(ambda)
    l = numpy.linspace(-b, b, n_lambda)
    dl = l[1]-l[0]
    l_i = l[:, numpy.newaxis]
    l_j = l[numpy.newaxis, :]

    if not singlet:
        # Calculate discretized integrands and right-hand sides as a function of indexes i,j
        # Eq. 22
        integrand11 = numpy.eye(len(k))
        integrand12 = - numpy.cos(k_i) * k_func(numpy.sin(k_i) - l_j) * dl
        rhs1 = numpy.ones(len(k))/(2*numpy.pi)

        # Eq. 23
        integrand21 = -k_func(numpy.sin(k_j) - l_i) * dk
        integrand22 = numpy.eye(len(l)) + k2_func(l_i - l_j) * dl
        rhs2 = numpy.zeros(len(l))

        # Compose the linear system and solve it
        A = numpy.block([[integrand11, integrand12], [integrand21, integrand22]])
        b = numpy.concatenate((rhs1, rhs2))
        x = numpy.linalg.solve(A, b)

        # Decompose the solution into unknowns rho, sigma
        rho = x[:len(k)]
        sigma = x[len(k):]
    else:
        # Calculate discretized integrands and right-hand sides as a function of indexes i,j
        # Eq. 2.8 Shiba
        integrand11 = numpy.eye(len(k))
#        integrand12 = - numpy.cos(k_i)*(4/u)*numpy.array([[R(4*(numpy.sin(x_i)-numpy.sin(x_j))/u) for x_j in k] for x_i in k])*dk
        integrand12 = - numpy.cos(k_i)*(4/u)* R(4*(numpy.sin(k_i)-numpy.sin(k_j))/u) *dk
        rhs1 = numpy.ones(len(k))/(2*numpy.pi)
        # Compose the linear system and solve it
        A = integrand11 + integrand12
        b = rhs1
        x = numpy.linalg.solve(A, b)
        rho = x
        sigma = numpy.zeros(n_lambda)
    # Calculate quantities
    # Eq. 19
    particles_per_site = sum(rho)*dk
    spin_downs_per_site = sum(sigma)*dl
    spin_ups_per_site = particles_per_site - spin_downs_per_site
    magnetization = (spin_ups_per_site - spin_downs_per_site) / particles_per_site
    # Eq. 25
    energy = - 2*sum(rho*numpy.cos(k)) * dk

    return dict(
        energy=energy,
        magnetization=magnetization,
        particles_per_site=particles_per_site,
    )


def E_bessel(U:float,epsabs=1.e-8):
    r'''
    Calculate GS energy of half-filled BA solution
    '''
    return -4*integrate.quad(lambda w: bessel(0.,w)*bessel(1.,w)/(w*(1+numpy.exp(w*U/2))), 0, numpy.inf,limit=200,epsabs=epsabs)[0]

def beta(U:float):
    return root_scalar(lambda beta:-2.*beta*numpy.sin(numpy.pi/beta)/numpy.pi + 4*integrate.quad(lambda w: bessel(0.,w)*bessel(1.,w)/(w*(1+numpy.exp(w*U/2))), 0, numpy.inf)[0],bracket=(1.,2.)).root

def mu_BA(U:float):
    r'''
    eq (65)
    '''
    return 2 - 4*integrate.quad(lambda w: bessel(1.,w)/(w*(1+numpy.exp(w*U/2))), 0, numpy.inf)[0]

def delta_BA(U:float):
    return U+4*numpy.cos(numpy.pi/beta(U))


def batch(u,nq=20,b=200, n_k=1000, n_lambda=1000,singlet=True):
    q_list=numpy.linspace(0,numpy.pi,nq)
    dbatch={'energy':numpy.zeros(nq),\
            'magnetization':numpy.zeros(nq),\
            'particles_per_site':numpy.zeros(nq)}
    for i,q in enumerate(q_list):
        d=solve(u,q=q,b=b,n_k=n_k,n_lambda=n_lambda,singlet=singlet)
        dbatch['energy'][i]=d['energy']
        dbatch['magnetization'][i]=d['magnetization']
        dbatch['particles_per_site'][i]=d['particles_per_site']
    return dbatch

def E_asymp(U):
    r'''
    Return strongly correlated limit energy
    E\prop t^2/U -> W/Ekin = -0.5 (Helmann feynmann)
    '''
    return -4*numpy.log(2)/U


