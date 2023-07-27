# pyphotonics

PyPhotonics is a post-processing python code that calculates photonic properties of materials. Based on the outcome of DFT, constrained DFT and vibrational calculations using DFT performed using VASP for a defect system, PyPhotonics uses the results in the output files and calculates the Huang-Rhys factor of the defect and the photoluminescence line-shape. Soon, the code will calculate the carrier capture coefficient and carrier lifetimes for defects, which are essential quantities for assessing the photovoltaic efficiency of materials.

If you wish to use PyPhotonics, please cite our paper:

- [Sherif AbdulkaderTawfik, Salvy P.Russo, PyPhotonics: A python package for the evaluation of luminescence properties of defects, Computer Physics Communications, 2022, 273, 108222.](https://www.sciencedirect.com/science/article/pii/S0010465521003349)


# Installation

You can install PyPhotonics using the pip command: `pip install pyphotonics`. The following python packages are required:

- scipy
- numpy
- pandas
- matplotlib
- oganesson

# Directory structure

The `pyphotonics` module is composed of the following files:

```
/pyphotonics
    /__init__.py
    /constants.py: A list of physical constants.
    /photoluminescence.py: Contains the Photoluminescence class, which drives the photonics calculations.

# How to use

To calculate the photonics properties of a crystal using pyphotonics, a number of DFT calculations should be performed with VASP first:
- The ground state structure of the crystal should be optimized. Let's call the output file CONTCAR_GS.
- The excited state structure of the crystal should be optimized. Let's call the output file CONTCAR_ES. For a tutorial on how to setup such calculation, see this [tutorial](./vasp_constrained_occupations.md)
- The phonon modes of the ground state system should be calculated with VASP, and the `bands.yaml` file should be produced using the `phonopy` code.

Once all of the above is done, you can calculate the HR factor and PL line-shape as follows (the complete example is in the `test/` directory):

```
from pyphotonics.photoluminescence import Photoluminescence
modes = 189 #number of modes
m = <<list of atomic masses>>
path_phonopy = './' #path to the bands.yaml file
path = './' #path to the structure files
p = Photoluminescence(path_phonopy,
                        path + "CONTCAR_GS",
                        path + "CONTCAR_ES",
                        189, "phonopy", m, 1000)

```

# How to setup INCAR for CDFT

Instructions and examples are [here](vasp_constrained_occupations.md).

