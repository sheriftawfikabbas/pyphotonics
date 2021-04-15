
import numpy as np
import sys
import matplotlib.pyplot as plt
import importlib
import photoluminescence
importlib.reload(photoluminescence)
from photoluminescence import Photoluminescence


# if __name__ == "__main__":
m = np.zeros(63)
for i in range(63):
    m[i] = 12.011 * 1.660539040e-27

m[62] = 14.007 * 1.660539040e-27

if True:
    sys = "NV+"
    path = "/mnt/c/MyJava/CMTLab_V3/Work_VII/QuantumEmission/Diamond/VASP/PBE/"
    path_phonopy = "/mnt/c/MyJava/CMTLab_V3/Work_VII/QuantumEmission/Diamond/VASP/PBE/Phonons_diamond_2x2x2_NV+_Download/"
    p = Photoluminescence(path_phonopy,
                            path + "diamond_2x2x2_NV+/positions.xyz",
                            path + "diamond_2x2x2_NV+_HOMO_LUMO/positions.xyz",
                            189, "phonopy", m,1000)
else:
    sys = "NV-"
    path = "/mnt/c/MyJava/CMTLab_V3/Work_VII/QuantumEmission/Diamond/VASP/PBE/"
    path_phonopy = "/mnt/c/MyJava/CMTLab_V3/Work_VII/QuantumEmission/Diamond/VASP/PBE/diamond_2x2x2_NV-_PHON/Results/"
    p = Photoluminescence(path_phonopy,
                        path + "diamond_2x2x2_NV-_PHON/Results/POSITIONS.xyz",
                        path + "diamond_2x2x2_NV-_HOMO_LUMO_PHON/Results/POSITIONS.xyz",
                        189, "phonopy", m, 1000)


print("Delta_R=", p.Delta_R)
print("Delta_Q=", p.Delta_Q)
print("HuangRhyes=", p.HuangRhyes)

plt.figure(figsize=(10, 5))
plt.plot(p.S_omega)
plt.ylabel('S')
plt.xlabel('Energy (meV)')
plt.xlim(0, 200)
# plt.ylim(0, 0.01)
plt.savefig('S_omega', bbox_inches='tight')

p.write_S('S')

A, I = p.PL(2, 2, 1.95)

plt.figure(figsize=(10, 10))
plt.plot(I.__abs__())
plt.ylabel('I')
plt.xlabel('Points')
plt.xlim(1200, 2000)
x_values, labels = plt.xticks()
labels = [float(x)/p.resolution for x in x_values]
plt.xticks(x_values, labels)
plt.ylim(0,600)
plt.savefig('I_'+sys, bbox_inches='tight')


plt.figure(figsize=(10, 10))
plt.plot(A.__abs__())
plt.ylabel('A')
plt.xlabel('Points')
plt.xlim(1200, 2000)
x_values, labels = plt.xticks()
labels = [float(x)/p.resolution for x in x_values]
plt.xticks(x_values, labels)
plt.ylim(0,600)
plt.savefig('A_'+sys, bbox_inches='tight')
plt.close()


from ase.build import bulk
from ase.io import read
from ase.lattice import CUB
diamond = read('C_mp-66_conventional_standard.cif')
lat = diamond.cell.get_bravais_lattice()
print(list(lat.get_special_points()))
path = diamond.cell.bandpath('GMRX', npoints=100)