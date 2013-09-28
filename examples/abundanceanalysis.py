# Example: calc_abundance()
#     determine ionic abundance from observed 
#     flux intensity for gievn electron density 
#     and temperature using  calc_abundance function
#     from pyEQUIB
# 
import pyequib

ion='o_iii'
tempi=10000.0
densi=5000.0
levels5007='3,4/'
iobs5007=1200.0
Abb5007=pyequib.cel.calc_abundance(ion, levels5007, tempi, densi, iobs5007)
print(Abb5007)
