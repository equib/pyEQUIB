# Example: calc_temp_dens()
#     determine electron density or temperature from given 
#     flux intensity ratio for a fixed electron density 
#     or temperature using  calc_temp_dens function
#     from pyEQUIB
# 
import pyequib

ion='s_ii'
levu='1,2,1,3/'
levl='1,5/'
diagtype='T'
dens = 2550.0
niiTratio=10.753
temp=pyequib.cel.calc_temp_dens(ion, levu, levl, niiTratio, diagtype, dens)
print(temp)

ion='s_ii'
levu='1,2/'
levl='1,3/'
diagtype='D'
temp = 7000.0
siiNratio=1.506
dens=pyequib.cel.calc_temp_dens(ion, levu, levl, siiNratio, diagtype, temp)
print(dens)
