# Deck amplitude with advanced pi-p interaction.

For Deck amplitude, requires GSL integration routines. Download gsl library at https://www.gnu.org/software/gsl/
To compile, do e.g.
```bash
g++ -std=c++11 deck.cpp -lgsl
```

## pion-proton amplitude
The model for the pion-nucleon scatterings is published in [Mat15b](http://cgl.soic.indiana.edu/jpac/PiN.php#mat15b).
It uses SAID low energy parametrization of the partial waves and Regge model for the high energy behavior.
 - A short description of the model can be found in [JPAC webpage](http://cgl.soic.indiana.edu/jpac/PiN.php#mat15b).
 - The function referred as `A`, `B`, and `C` in the code (see `pion_nucleon` folder) are `A`, `B`, and `A'` in the description

## Deck model
The code incorporates `pi-n` amplitude to the Deck model following the [Ascoli paper](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.8.3894).
 - The matrix element in the GJ frame, i.e. functions `M++`, `M--`, `M+-`, and `M-+` can be called (see `deck/deck_pi_n`)
 - Also, M-odd and M-even combinations (see `M_pm_m_mp` which is `M+- - M-+`)
