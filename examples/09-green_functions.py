from pyhub.solver.fermi_hubbard import FermiHubbard
import numpy as np
import matplotlib.pyplot as plt
import os 

ICGMmarine=(0.168,0.168,0.525)
ICGMblue=(0,0.549,0.714)
ICGMorange=(0.968,0.647,0)
ICGMyellow=(1,0.804,0)
#clr=[map1,map2,map3,map4,map5]
clr=[ICGMmarine,ICGMorange,ICGMblue,ICGMyellow]
plt.matplotlib.rcParams.update({'figure.figsize': (12, 10),'figure.autolayout': False,\
    'font.size':30,'lines.linewidth':2.5,'lines.markersize':0.01,'lines.marker':'*'})

np.set_printoptions(precision=4)
nb_sites=6
nup = ndown = nb_sites//2
T = 0.

t_matrix = np.diag(np.full(nb_sites-1,-1.),k=1) + np.diag(np.full(nb_sites-1,-1.),k=-1)
t_matrix[0,-1],t_matrix[-1,0] = -1.,-1.  

U = 4.

FH = FermiHubbard(nb_sites,nup,ndown,t_matrix,U,T=T)
FH.kernel(max_lcz=1000,acc_lcz = 1.e-8,nb_comp_states=1,\
    compute_rq=False,verbose=True)
FH.spgf()

print(f'mu :{FH.mu["up"]}')
print(f'ae :{FH.ae["up"]}')
print(f'ip :{FH.ip["up"]}')
    

w_grid = np.linspace(-7,7,501)+FH.mu['up']
wticks = np.delete(np.around(np.linspace(w_grid[0],w_grid[-1],5),1),2)
plt.plot(w_grid-FH.mu['up'], FH.gf['up'][0,0].spectral_density(w_grid,eta=0.03), color = clr[0],label=r'$\mathcal{S}(G^R(w))$')
plt.plot(w_grid-FH.mu['up'], FH.gf0['up'][0,0].spectral_density(w_grid-FH.mu['up'],eta=0.03), color = clr[1],label=r'$\mathcal{S}(G^{0R}(w+\mu))$',linestyle='-')
plt.ylabel(r'$A(w)$',fontsize='30') 
plt.xlabel(r'$w-\mu$',fontsize='30')
plt.xticks(np.concatenate((wticks,[-FH.ip['up']-FH.mu['up'],0.,-FH.ae['up']-FH.mu['up']])),np.concatenate((wticks,[r'$-IP$',r'',r'$-AE$'])))
plt.grid(True,which="both", linestyle='--')
plt.legend(loc='upper right',fontsize='25') 
plt.tight_layout()
plt.savefig('spec_dens')

w_grid = np.linspace(-9,9,500)+FH.mu['up']
plt.clf()
reim = FH.self_energy(w_grid,eta=0.03)['up'][0,0]
re,im = np.real(reim),np.imag(reim)
plt.plot(w_grid-FH.mu['up'], re-FH.mu['up'], color = clr[0],label=r'$\mathcal{R}(\Sigma^R(w))-\mu$')
plt.plot(w_grid-FH.mu['up'], -im, color = clr[1],label=r'$-\mathcal{I}(\Sigma^R(w))$')
#    plt.ylabel(r'$A(w)$',fontsize='30') 
plt.xlabel(r'$w-\mu$',fontsize='30')
plt.grid(True,which="both", linestyle='--')
plt.legend(loc='upper right',fontsize='25') 
plt.tight_layout()
plt.savefig('self_energy')

w_grid = np.linspace(-7,7,801)+FH.mu['up']
plt.clf()
eta=0.03
A = np.imag(FH.gf['up'][0,0].retarded(w_grid-1j*(eta-1.e-2),eta=eta))
A/=max(abs(A))
plt.plot(w_grid-FH.mu['up'], A-eta, color = clr[0],label=r'$\mathcal{I}(G^R(iw))$')
A = np.imag(FH.gf['up'][0,0].advanced(w_grid+1j*(eta-1.e-2),eta=eta))
A/=max(abs(A))
plt.plot(w_grid-FH.mu['up'], A+eta, color = clr[1],label=r'$\mathcal{I}(G^A(iw))$')
A = np.imag(FH.gf['up'][0,0].causal(w_grid+1j*(eta-1.e-2)*np.sign(FH.mu['up']-w_grid),FH.mu['up'],eta=eta))
A/=max(abs(A))
plt.plot(w_grid-FH.mu['up'], A+eta*np.sign(FH.mu['up']-w_grid), color ='black',label=r'$\mathcal{I}(G^C(iw))$',linestyle='--')
plt.text(-5.,-0.3,'<',color=clr[0])
plt.text(5.,0.3,'>',color=clr[1])
plt.plot(w_grid-FH.mu['up'],np.full(len(w_grid),1.e-3),color='black',linestyle='--',linewidth=1.)
r,x = (w_grid[-1] - w_grid[0])/2,w_grid-FH.mu['up']
plt.plot(x,np.full(len(w_grid),1.e-3)+0.98*np.sqrt(r**2-x**2)/r,color='black',linestyle='--',linewidth=1.)
plt.text(-5.,0.85,r'$\mathcal{C}$',color='black')
plt.ylabel(r'$iw$',fontsize='30') 
plt.ylim(-1,1.)
plt.xlim(w_grid[0]-FH.mu['up'],w_grid[-1]-FH.mu['up'])
plt.xlabel(r'$w - \mu$',fontsize='30')
plt.grid(True,which="both", linestyle='--')
plt.legend(loc='upper right',fontsize='25') 
plt.tight_layout()
plt.savefig('propagator')

os.remove('solver.h5')