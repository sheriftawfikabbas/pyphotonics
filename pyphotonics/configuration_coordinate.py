import sys
import os

import argparse
import numpy as np

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.inputs import Poscar, Kpoints
# from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.structure import Structure
import matplotlib.pyplot as plt


class ConfigurationCoordinate:
    def read_poscar(self,i_path, l_get_sorted_symbols=False):
        poscar = Poscar.from_file("{}".format(i_path))
        struct = poscar.structure
        if l_get_sorted_symbols:
            return struct, poscar.site_symbols
        else:
            return struct


    def Delta_Q(self,i_file, f_file, disp_range=None):
        struct_i, sorted_symbols = self.read_poscar(i_file, True)
        struct_f, sorted_symbols = self.read_poscar(f_file, True)
        delta_R = struct_f.frac_coords - struct_i.frac_coords
        delta_R = (delta_R + 0.5) % 1 - 0.5
        
        lattice = struct_i.lattice.matrix
        delta_R = np.dot(delta_R, lattice)

        masses = np.array([spc.atomic_mass for spc in struct_i.species])
        delta_Q2 = masses[:,None] * delta_R ** 2
        return np.sqrt(delta_Q2.sum())

    def get_init_fin(self,i_file, f_file, disp_range=np.linspace(-1, 1, 11), output_dir='disp_dir'):
        '''
        '''
        # A. Alkauskas, Q. Yan, and C. G. Van de Walle, Physical Review B 90, 27 (2014)
        struct_i, sorted_symbols = self.read_poscar(i_file, True)
        struct_f, sorted_symbols = self.read_poscar(f_file, True)
        delta_R = struct_f.frac_coords - struct_i.frac_coords
        delta_R = (delta_R + 0.5) % 1 - 0.5
        
        lattice = struct_i.lattice.matrix #[None,:,:]
        delta_R = np.dot(delta_R, lattice)

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        masses = np.array([spc.atomic_mass for spc in struct_i.species])
        delta_Q2 = masses[:,None] * delta_R ** 2

        print('Delta_Q^2', np.sqrt(delta_Q2.sum()))

        for frac in disp_range:
            disp = frac * delta_R
            print(disp[0][0])
            struct = Structure(struct_i.lattice, struct_i.species, \
                            struct_i.cart_coords + disp, \
                            coords_are_cartesian=True)
            Poscar(struct).write_file('{0}/POSCAR_{1:03d}'.format(output_dir, int(np.rint(frac*10))))

# path = "/mnt/c/MyJava/CMTLab_V3/Work_VII/QuantumEmission/Diamond/VASP/PBE/"
# s_i = path + "diamond_2x2x2_NV-/CONTCAR"
# s_f = path + "diamond_2x2x2_NV-_HOMO_LUMO/CONTCAR"

# cc = ConfigurationCoordinate()
# # cc.get_init_fin(s_i,s_f)

# f = cc.read_poscar(s_f) 
# i = cc.read_poscar(s_i) 

# i.translate_sites(range(len(i.frac_coords)),[0,0,0.1], frac_coords=False)
# f.translate_sites(range(len(f.frac_coords)),[0,0,0.1], frac_coords=False)

# f.frac_coords-i.frac_coords

# from ase.geometry import wrap_positions
# f_w=wrap_positions(f.frac_coords,f.lattice.matrix,pretty_translation=True,center=(0,0,0))
# i_w=wrap_positions(i.frac_coords,i.lattice.matrix,pretty_translation=True,center=(0,0,0))
# f_w-i_w
# f.frac_coords-i.frac_coords
