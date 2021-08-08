# pyphotonics

PyPhotonics is a post-processing python code that calculates photonic properties of materials. Based on the outcome of DFT, constrained DFT and vibrational calculations using DFT performed using VASP for a defect system, PyPhotonics uses the results in the output files and calculates the Huang-Rhys factor of the defect and the photoluminescence line-shape. Soon, the code will calculate the carrier capture coefficient and carrier lifetimes for defects, which are essential quantities for assessing the photovoltaic efficiency of materials.

# Installation

You can install PyPhotonics using the pip command: `pip install pyphotonics`. The following python packages are required:

- scipy
- numpy
- pandas
- matplotlib

# Directory structure

The `pyphotonics` module is composed of the following files:

```
/pyphotonics
    /__init__.py
    /configuration_coordinate.py: Contains the class ConfigurationCoordinate which calculates the Huang-Rhys factor and other quantities.
    /constants.py: A list of physical constants.
    /photoluminescence.py: Contains the Photoluminescence class, which drives the photonics calculations.
    /schrodinger.py: Contains the Schrodinger class which solves the 1-dimensional Schrodinger equation for an arbitrary potential. Will be used in a future version of the code.
    /xyz.py: Contains the XYZ class which stores coordinate information of the crystals.
```

# How to use

To calculate the photonics properties of a crystal using pyphotonics, a number of DFT calculations should be performed with VASP first:
- The ground state structure of the crystal should be optimized. Let's call the output file CONTCAR_GS.
- The excited state structure of the crystal should be optimized. Let's call the output file CONTCAR_ES. For a tutorial on how to setup such calculation, see: 
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
                        189, "phonopy", m, 1000, shift_vector=[0.0, 0, 0.1])

```

