import numpy as np
from scipy.integrate import quad
from scipy.optimize import root_scalar

__author__ = "Quentin Marécat"
__maintainer__ = "Quentin Marécat"
__email__ = "quentin.marecat@gmail.com"
__date__ = "March, 2023"
__cite__ = r'''
        @article{gutzwiller_correlation_1965,
                title = {Correlation of Electrons in a Narrow s Band},
                volume = {137},
                issn = {0031-899X},
                url = {https://link.aps.org/doi/10.1103/PhysRev.137.A1726},
                doi = {10.1103/PhysRev.137.A1726},
                pages = {A1726--A1735},
                number = {6},
                journal = {Phys. Rev.},
                author = {Gutzwiller, Martin C.},
                date = {1965-03-15},
        }
'''

r'''
L = number of lattice sites
eta = weight to each configuration depending on the amount of crowding <= 1 (uncorrelated)
m = number of lattoce sites to be occupied by \uparrow particules
mu = number of lattoce sites to be occupied by \downarrow particules
non-ferro state -> m=mu
nu = number of identical liattice sites among the sets (g_1, ..., g_\m) and (\gamma_1, ..., \gamma_\mu)
'''
def q(m,nu,mu=None):
    r'''
    Eq (B6)
    '''
    if mu == None:
        return q_special(m,nu)
    else:
        eta=np.sqrt(1./((m-nu)*(mu-nu)/(nu*(1-m-mu+nu))))
        return ((m-nu)*(1-m-mu+nu))*(1+((mu-nu)*eta)/(1-m-mu+nu))**2/(m*(1-m))

def dq(m,nu,mu=None):
    r'''
    Eq (B6)
    '''
    if mu == None:
        return dq_special(m,nu)
    else:
        eta=np.sqrt(1./((m-nu)*(mu-nu)/(nu*(1-m-mu+nu))))
        deta = ((m-nu)*(mu-nu)*(1-m-mu+2*nu)-(nu*(1-m-mu+nu)*(-m-mu+2*nu)))/(2*np.sqrt(eta)*(m-nu)**2*(mu-nu)**2)
        a,b=((m-nu)*(1-m-mu+nu))/(m*(1-m)),(1+((mu-nu)*eta)/(1-m-mu+nu))
        return b**2*(mu+2*m-1-2*nu)/(m*(1-m))+2*a*b*((1-m-mu+nu)*(-eta+(mu-nu)*deta) - (mu-nu)*eta)/((1-m-mu+nu)**2)

def q_special(m,nu):
    r'''
    m=\mu -> non-ferro state
    Eq (B8)
    '''
    return (m-nu)*(np.sqrt(1-2*m+nu)+np.sqrt(nu))**2/(m*(1-m))

def dq_special(m,nu):
#    return (q(m,nu+1.e-5)-q(m,nu-1.e-5))/2.e-5
    a,b = np.sqrt(1-2*m+nu),np.sqrt(nu)
    return (-(a+b)**2 + (m-nu)*(2+(a**2+b**2)/(a*b)))/(m*(1-m))


def E(e0,U,m,nu,mu=None):
    return 2*e0*q(m,nu,mu=mu)*m + nu*U

def E_min(e0,U,m,mu=None):
    try:
        nu = root_scalar(lambda nu:dq(m,nu)+U/(2*m*e0),bracket=(0.,0.25)).root
        return 2*e0*q(m,nu,mu=mu)*m + nu*U
    except ValueError:
        return 0.

def E_saddle_point(e0,m,nu,mu=None):
    return 2*e0*(q_special(m,nu,mu=mu)*m-nu*dq(m,nu,mu=mu))

if __name__=='__main__':
    from matplotlib import colors
    import matplotlib.pyplot as plt
    import numpy as np

    ICGMmarine=(0.168,0.168,0.525)
    ICGMblue=(0,0.549,0.714)
    ICGMorange=(0.968,0.647,0)
    ICGMyellow=(1,0.804,0)
    gray=(0.985,0.985,0.985)
    map1=tuple(np.array([255, 195, 15])/256)
    map2=tuple(np.array([255, 87, 51])/256)
    map3=tuple(np.array([199, 0, 57])/256)
    map4=tuple(np.array([144, 12, 63])/256)
    map5=tuple(np.array([84, 24, 69])/256)
    clr_map=[ICGMyellow,ICGMorange,ICGMblue,ICGMmarine]
    clr = [ICGMmarine,ICGMorange,ICGMblue,ICGMyellow]
    cmap = colors.LinearSegmentedColormap.from_list('my_list', clr_map, N=100)
    default_color='k'

    plt.matplotlib.rcParams.update({'figure.figsize': (12, 10),'figure.autolayout': False,\
        'font.size':30,'lines.linewidth':2.5,'lines.markersize':0.01,'lines.marker':'*','image.cmap':cmap,\
        'savefig.format':'png',\
        'grid.color':'k','grid.linestyle':'--',\
        'text.color':default_color,'axes.labelcolor':default_color,'xtick.color':default_color,'ytick.color':default_color,'grid.color':default_color})

    def e_1d(m):
        r'''
        e^0 = \sum_k e_k/m
        '''
        return np.sum([-2*np.cos(k) for k in np.linspace(0.,np.pi,402)[:int(m*402)]])/(m*402)
        
    from dos import semicircular
    def e_dinf(m):
        mu = root_scalar(lambda b:m-quad(lambda w:semicircular(w),-np.inf,b)[0],bracket=(-2,2)).root
        return quad(lambda w:w*semicircular(w),-np.inf,mu)[0]/m

    nu = np.linspace(0.,0.25,100)
    for i,m in enumerate([0.1,0.2,0.3,0.4,0.5]):
        q_array=np.array([q(m,n) for n in nu])
        index = np.where((q_array[1:]-q_array[:-1])/(nu[1]-nu[0])>0)[0]
        plt.plot(nu[index],q_array[index],color=cmap(i/5),label=f'm = {np.around(m,1)}')
    plt.xlim(0.,0.26)
    plt.ylim(0.,1.05)
    plt.xticks([0.05,0.1,0.15,0.2,0.25])
    plt.yticks([0.,0.2,0.4,0.6,0.8,1.])
    plt.ylabel(r'$q$')
    plt.xlabel(r'$\nu$')
    plt.title('Fig. 4')
    plt.legend(loc='lower right')
    plt.grid(color='k',linestyle='--')
    plt.tight_layout()
    plt.savefig('gutz_occ')

    plt.close()
    nu = np.linspace(0.0001,0.25,100)
    for i,m in enumerate([0.1,0.2,0.3,0.4,0.5]):
        plt.plot(nu,[q(m,n)-n*dq(m,n) for n in nu],color=cmap(i/5),label=f'm = {np.around(m,1)}')
    plt.xlim(0.,0.25)
    plt.ylim(0.,1.)
    plt.xticks([0.05,0.1,0.15,0.2,0.25])
    plt.yticks([0.,0.2,0.4,0.6,0.8,1.])
    plt.ylabel(r'$q - \nu \partial q / \partial \nu$')
    plt.xlabel(r'$\nu$')
    plt.title('Fig. 2')
    plt.legend(loc='lower right')
    plt.grid(color='k',linestyle='--')
    plt.tight_layout()
    plt.savefig('gutz_renorm')

    plt.close()
    nu = np.linspace(0.0001,0.25,100)
    for i,m in enumerate([0.1,0.2,0.3,0.4,0.5]):
        plt.plot(nu,[m*dq(m,n) for n in nu],color=cmap(i/5),label=f'm = {np.around(m,1)}')
    plt.xlim(0.,0.25)
    plt.ylim(0.,6.)
    plt.xticks([0.05,0.1,0.15,0.2,0.25])
    plt.yticks([0.,1.,2.,3.,4.,5.])
    plt.ylabel(r'$m \partial q / \partial \nu$')
    plt.xlabel(r'$\nu$')
    plt.title('Fig. 3')
    plt.legend(loc='upper right')
    plt.grid(color='k',linestyle='--')
    plt.tight_layout()
    plt.savefig('gutz_ratio')

    plt.close()
    u_range = np.linspace(0.0,1.,100)
    for i,m in enumerate([0.1,0.2,0.3,0.4,0.499]):
        #e0 = e_dinf(m)
        e0 = e_1d(m)
        plt.plot(u_range,[root_scalar(lambda nu:dq(m,nu)+U/(2*m*e0),bracket=(0.,0.25)).root for U in 4*u_range/(1-u_range)],color=cmap(i/5),label=f'm = {np.around(m,1)}')
    plt.xlim(0.,1.)
    plt.ylim(0.,0.25)
    plt.xticks([0.,0.2,0.4,0.6,0.8,1.])
    plt.yticks([0.05,0.1,0.15,0.2,0.25])
    plt.ylabel(r'$\nu$')
    plt.xlabel(r'$U/(U+2B)$')
    plt.legend(loc='upper right')
    plt.grid(color='k',linestyle='--')
    plt.tight_layout()
    plt.savefig('gutz_ni')

    plt.close()
    u_range = np.linspace(0.0,1.,100)
    for i,m in enumerate([0.1,0.2,0.3,0.4,0.499]):
        #e0 = e_dinf(m)
        e0 = e_1d(m)
        nu = [root_scalar(lambda nu:dq(m,nu)+U/(2*m*e0),bracket=(0.,0.25)).root for U in 4*u_range/(1-u_range)]
        plt.plot(u_range,[q(m,nu[i]) for i,U in enumerate(4*u_range/(1-u_range))],color=cmap(i/5),label=f'm = {np.around(m,1)}')
    plt.xlim(0.,1.)
    plt.ylim(0.,1.05)
    plt.xticks([0.,0.2,0.4,0.6,0.8,1.])
    plt.yticks([0.,0.2,0.4,0.6,0.8,1.])
    plt.ylabel(r'$q$')
    plt.xlabel(r'$U/(U+2B)$')
    plt.legend(loc='lower left')
    plt.grid(color='k',linestyle='--')
    plt.tight_layout()
    plt.savefig('gutz_q')

    plt.close()
    plt.plot(np.linspace(-4,4,500)/2,semicircular(np.linspace(-4,4,500))*np.pi,color=ICGMmarine)
    plt.xlabel(r'$w/B$')
    plt.grid(color='k',linestyle='--')
    plt.tight_layout()
    plt.savefig('semi')
