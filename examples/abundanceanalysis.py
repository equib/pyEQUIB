# Example: calc_abundance()
#     determine ionic abundance from observed 
#     flux intensity for gievn electron density 
#     and temperature using  calc_abundance function
#     from pyEQUIB
# 
import pyequib

ion='o_iii'
temperature=10000.0
density=5000.0
levels5007='3,4/'
iobs5007=1200.0
emis=pyequib.cel.calc_emissivity(temperature=temperature, density=density, ion=ion, levels=levels5007)
print(emis)

ion='o_iii'
temperature=10000.0
density=5000.0
levels5007='3,4/'
iobs5007=1200.0
Abb5007=pyequib.cel.calc_abundance(temperature=temperature, density=density, line_flux=iobs5007, ion=ion, atomic_levels=levels5007)
print(Abb5007)
