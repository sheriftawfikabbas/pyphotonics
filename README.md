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
```

# How to use

## PL line-shape and the Huang-Rhys factor

To calculate the PL line-shape and the Huang-Rhys of a crystal structure using pyphotonics, a number of DFT calculations should be performed with VASP first:
- The ground state structure of the crystal should be optimized. Let's call the output file CONTCAR_GS.
- The excited state structure of the crystal should be optimized. Let's call the output file CONTCAR_ES. For a tutorial on how to setup such calculation, see this [tutorial](./vasp_constrained_occupations.md)
- The phonon modes of the ground state system should be calculated with VASP, and the `bands.yaml` file should be produced using the `phonopy` code.

Once all of the above is done, you can calculate the HR factor and PL line-shape as follows (the complete example is in the `test/` directory):


```
from pyphotonics.photoluminescence import Photoluminescence
p = Photoluminescence("CONTCAR_GS", "CONTCAR_ES", 189, method="phonopy", resolution=1000)

```

You can also run `pyphotonics` from the command line by typing the following in the directory `test/photoluminscence/`:

```
pyphotonics -cgs CONTCAR_GS -ces CONTCAR_ES -m 189 -M phonopy -r 1000
```


## Setting up the INCAR for CDFT

VASP allows you to excite a crystal structure by constraining the occupations of the electronic bands. Let's say you want to excite a spin-up (spin component 1) or spin-down (spin component 2) electron from the ground state to the excited state. To do this, run the following command:

```
pyphotonics-incar <path to your OUTCAR file>
```


`pyphotonics` will write to the terminal the values of the two INCAR tags, `FERWE` and `FERDO`, which you can then copy-paste into your INCAR file to perform the excited state optimisation.

There is a sample `OUTCAR` file in the directory `test/ferwe_ferdo` which you can test with. Go to the directory and type: `pyphotonics-incar`.

### What are the `FERWE` and `FERDO`?

These are INCAR tags that specify the electronic occupation. Here is how you can construct those tags for any structure:

- Find out how many bands do you have. You can obtain that by search for `NBANDS` in the OUTCAR file
- Find out how many k points do you have. You can look that up in the OUTCAR by searching for the section that starts with the string spin component 1 in the OUTCAR. Let's call that number `N_K`.
- Find out how many electrons are occupying spin components 1 and 2 by searching the **last** such section in the OUTCAR that starts with spin component 1. Let's say there are `O_UP` electrons that occupy the spin up, `O_DN` electrons that occupy the spin down.
- Next, calculate how many states are empty in both spins: empty states in the spin up = `U_UP` where `U_UP = NBANDS - O_UP`. Same applies for `U_DN`.

Finally, add the `FERWE` (spin up, or spin component 1) and `FERDO` (spin down, or spin component 2) as follows:

- `FERWE = [O_UP-1]*1.0 1*0.0 1*1.0 [U_UP-1]*0.0 <repeated N_K times>`
- `FERDO = [O_DN]*1.0 [U_DN]*0.0 <repeated N_K times>`

Note the `<repeated N_K times>` above. It means: repeat `[O_UP-1]*1.0 1*0.0 1*1.0 [U_UP-1]*0.0` for `N_K` times.

Let's have an example. Let's say we have 8 kpoints in both spins, 216 bands, where the spin up electrons occupy 144 bands and the spin down occupy 142 bands. Then, here are the two tags:

- `FERWE = 143*1.0 1*0.0 1*1.0 71*0.0 143*1.0 1*0.0 1*1.0 71*0.0 143*1.0 1*0.0 1*1.0 71*0.0 143*1.0 1*0.0 1*1.0 71*0.0 143*1.0 1*0.0 1*1.0 71*0.0 143*1.0 1*0.0 1*1.0 71*0.0 143*1.0 1*0.0 1*1.0 71*0.0 143*1.0 1*0.0 1*1.0 71*0.0`
- `FERDO = 142*1.0 74*0.0 142*1.0 74*0.0 142*1.0 74*0.0 142*1.0 74*0.0 142*1.0 74*0.0 142*1.0 74*0.0 142*1.0 74*0.0 142*1.0 74*0.0`

