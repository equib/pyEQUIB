"""Tests for pyequib.cel."""

import pyequib

ion='o_iii'
temperature=10000.0
density=5000.0
levels5007='3,4/'
iobs5007=1200.0
Abb5007=pyequib.cel.calc_abundance(temperature=temperature, density=density, line_flux=iobs5007, ion=ion, atomic_levels=levels5007)

ion='s_ii'
upper_levels='1,2,1,3/'
lower_levels='1,5/'
density = 2550.0
niiTratio=10.753
temp=pyequib.cel.calc_temperature(line_flux_ratio=niiTratio, density=density, ion=ion, upper_levels=upper_levels, lower_levels=lower_levels)

ion='s_ii'
upper_levels='1,2/'
lower_levels='1,3/'
temperature = 7000.0
siiNratio=1.506
dens=pyequib.cel.calc_density(line_flux_ratio=siiNratio, temperature=temperature, ion=ion, upper_levels=upper_levels, lower_levels=lower_levels)
