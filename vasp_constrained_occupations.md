# Constrained occupations in VASP

I am assuming here that the user is familiar with VASP.

For a spin-polarized calculation (where `ISPIN=2` in the INCAR), let's say you want to excite a spin-up (spin component 1) or spin-down (spin component 2) electron from the ground state to the excited state.

To do this, you will need to do the following:

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

I hope that helps!