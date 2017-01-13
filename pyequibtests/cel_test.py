"""Tests for pyequib.cel."""

import pyequib

ion='o_iii'
tempi=10000.0
densi=5000.0
levels5007='3,4/'
iobs5007=1200.0
Abb5007=pyequib.cel.calc_abundance(ion, levels5007, tempi, densi, iobs5007)

ion='s_ii'
levu='1,2,1,3/'
levl='1,5/'
diagtype='T'
dens = 2550.0
niiTratio=10.753
temp=pyequib.cel.calc_temp_dens(ion, levu, levl, niiTratio, diagtype, dens)

ion='s_ii'
levu='1,2/'
levl='1,3/'
diagtype='D'
temp = 7000.0
siiNratio=1.506
dens=pyequib.cel.calc_temp_dens(ion, levu, levl, siiNratio, diagtype, temp)

