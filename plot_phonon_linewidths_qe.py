#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import re


f1 = open('YB6_EPC_bands.dat', 'r').read()
d1 = re.findall(r'[-.\d]+', f1)
d1 = np.array(d1, dtype=float)
d1.shape = (-1,3)
Nk = len(set(d1[:,0]))
print('Nk', Nk)
print('d1.shape', d1.shape)
d1.shape = (-1,Nk,3)
Nbands1 = d1.shape[0]
vmax1 = np.max(d1[:,:,2])


f2 = open('LaB6_EPC_bands.dat', 'r').read()
d2 = re.findall(r'[-.\d]+', f2)
d2 = np.array(d2, dtype=float)
d2.shape = (-1,3)
Nk = len(set(d2[:,0]))
print('Nk', Nk)
print('d2.shape', d2.shape)
d2.shape = (-1,Nk,3)
Nbands2 = d2.shape[0]
vmax2 = np.max(d2[:,:,2])

vmax = max(vmax1, vmax2)


plt.figure()
plt.subplot(1, 2, 1)
for i in range(Nbands1):
    plt.plot(d1[i,:,0], d1[i,:,1]/33.356683, '-k')
    plt.scatter(d1[i,:,0], d1[i,:,1]/33.356683, s=d1[i,:,2], c=d1[i,:,2], cmap='jet', vmax=vmax1 )
plt.ylim(bottom=0, top=43)
plt.xlim([0,2.7802])
plt.grid()
plt.xticks([0, 0.7071, 1.2071, 1.9142, 2.7802], ['X', 'R', 'M', 'G', 'R'])
plt.ylabel('Frequency (THz)')

plt.subplot(1, 2, 2)
for i in range(Nbands2):
    plt.plot(d2[i,:,0], d2[i,:,1]/33.356683, '-k')
    plt.scatter(d2[i,:,0], d2[i,:,1]/33.356683, s=d2[i,:,2], c=d2[i,:,2], cmap='jet', vmax=vmax2 )
# the 33.356683 is to convert from cm-1 to THz

plt.ylim(bottom=0, top=43)
plt.xlim([0,2.7802])
plt.grid()
plt.xticks([0, 0.7071, 1.2071, 1.9142, 2.7802], ['X', 'R', 'M', 'G', 'R'])
plt.colorbar(label='EPC')

plt.show()
