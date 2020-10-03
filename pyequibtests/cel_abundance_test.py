# Example: calc_abundance()
#     determine ionic abundance from observed
#     flux intensity for gievn electron density
#     and temperature using  calc_abundance function
#     from proEQUIB
#
# --- Begin MAIN program. ---------------
#
#
import pyequib
import atomneb
import os
import numpy as np

# Locate datasets
base_dir = os.getcwd()
data_dir = os.path.join('externals', 'atomneb', 'atomic-data', 'chianti70')
data_rc_dir = os.path.join('externals', 'atomneb', 'atomic-data-rc')
atom_elj_file = os.path.join(base_dir,data_dir, 'AtomElj.fits')
atom_omij_file = os.path.join(base_dir,data_dir, 'AtomOmij.fits')
atom_aij_file = os.path.join(base_dir,data_dir, 'AtomAij.fits')
atom_rc_sh95_file = os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')

atom = 'h'
ion = 'ii' # H I Rec
hi_rc_data = atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)

atom = 'o'
ion = 'iii' # [O III]
o_iii_elj = atomneb.read_elj(atom_elj_file, atom, ion, level_num=5) # read Energy Levels (Ej)
o_iii_omij = atomneb.read_omij(atom_omij_file, atom, ion) # read Collision Strengths (Omegaij)
o_iii_aij = atomneb.read_aij(atom_aij_file, atom, ion) # read Transition Probabilities (Aij)

levels5007 = '3,4/'
temperature = np.float64(10000.0,)
density = np.float64(5000.0)
iobs5007 = np.float64(1200.0)
abb5007 = np.float64(0.0)

emis = pyequib.calc_emissivity(temperature=temperature, density=density, atomic_levels=levels5007, elj_data=o_iii_elj, omij_data=o_iii_omij, aij_data=o_iii_aij)
print('Emissivity(O III 5007):', emis)

abb5007 = pyequib.calc_abundance(temperature=temperature, density=density, line_flux=iobs5007, atomic_levels=levels5007,
                         elj_data=o_iii_elj, omij_data=o_iii_omij, aij_data=o_iii_aij, h_i_aeff_data=hi_rc_data.aeff)
print('N(O^2+)/N(H+):', abb5007)

nlj = pyequib.calc_populations(temperature=temperature, density=density, elj_data=o_iii_elj, omij_data=o_iii_omij, aij_data=o_iii_aij)
print('Atomic Level Populations:', nlj)


n_crit = pyequib.calc_crit_density(temperature=temperature, elj_data=o_iii_elj, omij_data=o_iii_omij, aij_data=o_iii_aij)
print('Critical Densities:', n_crit)

temperature = np.float64(10000.0)
omij_t = pyequib.get_omij_temp(temperature=temperature, omij_data=o_iii_omij, level_num=8)
print('Effective Collision Strengths: ')
print(omij_t)

pyequib.print_ionic(temperature=temperature, density=density,
            elj_data=o_iii_elj, omij_data=o_iii_omij, aij_data=o_iii_aij,
            h_i_aeff_data=hi_rc_data.aeff)
#
# --- End MAIN program. ---------------

