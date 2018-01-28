"""Tests for pyequib.cel."""

import pyequib

ion='o_iii'
tempi=10000.0
densi=5000.0
levels5007='3,4/'
iobs5007=1200.0
Abb5007=pyequib.cel.calc_abundance(temperature=tempi, density=densi, line_flux=iobs5007, ion=ion, atomic_levels=levels5007)

ion='s_ii'
levu='1,2,1,3/'
levl='1,5/'
dens = 2550.0
niiTratio=10.753
temp=pyequib.cel.calc_temperature(line_flux_ratio=niiTratio, density=dens, ion=ion, upper_levels=levu, lower_levels=levl)

ion='s_ii'
levu='1,2/'
levl='1,3/'
temp = 7000.0
siiNratio=1.506
dens=pyequib.cel.calc_density(line_flux_ratio=siiNratio, temperature=temp, ion=ion, upper_levels=levu, lower_levels=levl)
