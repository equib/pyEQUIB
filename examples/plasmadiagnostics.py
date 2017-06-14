# Example: calc_temperature() and calc_density()
#     determine electron density or temperature from given 
#     flux intensity ratio for a fixed electron density 
#     or temperature using  calc_temperature and 
#     calc_density functions from pyEQUIB
# 
import pyequib

ion='s_ii'
levu='1,2,1,3/'
levl='1,5/'
diagtype='T'
density = 2550.0
siiTratio=10.753
temperature=pyequib.cel.calc_temperature(line_flux_ratio=siiTratio, density=density, ion=ion, upper_levels=levu, lower_levels=levl)
print(temperature)

ion='s_ii'
levu=levu='1,2/'
levl='1,3/'
diagtype='D'
temperature = 7000.0
siiNratio=1.506
density=pyequib.cel.calc_density(line_flux_ratio=siiNratio, temperature=temperature, ion=ion, upper_levels=levu, lower_levels=levl)
print(density)

