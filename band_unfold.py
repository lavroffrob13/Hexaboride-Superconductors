import pyprocar
import numpy as np
import matplotlib.pyplot as plt

pyprocar.unfold(
    fname='PROCAR',
    poscar='POSCAR',
    outcar='OUTCAR',
    supercell_matrix=np.array([[2,0,0],[0,2,0],[0,0,1]]),
    ispin=1, # None for non-spin polarized calculation. For spin polarized case, ispin=1: up, ispin=2: down
    efermi=None,
    shift_efermi=True,
    elimit=(-2,4),
    kticks=[0, 49, 99],
    knames=['X','R', 'M', 'G', 'R'],
    print_kpts=False,
    #show_band=True,
    width=5,
    color='red',
    savetab='unfolding.csv',#    savefig='unfolded_band.png',
    exportplt=True,
    repair=False)
plt.show()

#KPOINTS for a 2x2x1 Pm3m:

#0.0 0.0 0.5 ! X
#1.0 1.0 0.5 ! R
#
#1.0 1.0 0.5 ! R
#1.0 1.0 0.0 ! M
#
#1.0 1.0 0.0 ! M
#0.0 0.0 0.0 ! G
#
#0.0 0.0 0.0 ! G
#1.0 1.0 0.5 ! R
