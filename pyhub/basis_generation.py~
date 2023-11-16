import numpy as np
import h5py
from hubbard import Solve

for nb_sites in range(15):
    with h5py.File('basis_hubbard.h5',  "a") as file:
        grp=file.create_group(f'{nb_sites}')
    for nb_elec in range(2*(nb_sites)):
        with h5py.File('basis_hubbard.h5',  "a") as file:
            grp=file[f'{nb_sites}'].create_group(f'{nb_elec}')
        for sz in np.linspace(-nb_elec/2,nb_elec/2,nb_elec+1,endpoint=True):
            print(nb_sites,nb_elec,sz)
            H = Solve(nb_sites,nb_elec,sz,t_matrix=np.diag(np.random.random(nb_sites)))
            H.run(max_lcz=300,acc_lcz = 1.e-2,nb_comp_states=1,compute_reduced_quantities=False)
            with h5py.File('basis_hubbard.h5',  "a") as file:
                grp=file[f'{nb_sites}'][f'{nb_elec}'].create_group(f'{np.around(sz,1)}')
                grp.attrs['nb_states'] = H.nb_states
                grp.attrs['nstates'] = H.nb_states
                grp.attrs['nup'] = int(float(nb_elec)/2+sz)
                grp.attrs['ndown'] = int(float(nb_elec)/2-sz)
                grp.attrs['nsup'] = H.nsup
                grp.attrs['nsdown'] = H.nsdown
                grp.create_dataset("basis_up", (H.nsup,), dtype='i',data=H.basis_up)
                grp.create_dataset("basis_down", (H.nsdown,), dtype='i',data=H.basis_down)
