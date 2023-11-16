from pyhub.solver.gutzwiller import *
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np
from pyhub.tools.dos import semicircular
import warnings
warnings.filterwarnings("ignore",category = RuntimeWarning)

ICGMmarine=(0.168,0.168,0.525)
ICGMblue=(0,0.549,0.714)
ICGMorange=(0.968,0.647,0)
ICGMyellow=(1,0.804,0)
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
