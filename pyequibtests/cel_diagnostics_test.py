# Example: calc_temperature() and calc_density()
#     determine electron density or temperature from given
#     flux intensity ratio for a fixed electron density
#     or temperature using calc_temperature function
#     calc_density function from proEQUIB
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

atom = 's'
ion = 'ii'
s_ii_elj = atomneb.read_elj(atom_elj_file, atom, ion, level_num=5) # read Energy Levels (Ej)
s_ii_omij = atomneb.read_omij(atom_omij_file, atom, ion) # read Collision Strengths (Omegaij)
s_ii_aij = atomneb.read_aij(atom_aij_file, atom, ion) # read Transition Probabilities (Aij)\

upper_levels = '1,2,1,3/'
lower_levels = '1,5/'
density = np.float64(2550)
line_flux_ratio = np.float64(10.753)
temperature = pyequib.calc_temperature(line_flux_ratio=line_flux_ratio, density=density, upper_levels=upper_levels, lower_levels=lower_levels, elj_data=s_ii_elj, omij_data=s_ii_omij, aij_data=s_ii_aij)
print("Electron Temperature:", temperature)

upper_levels = '1,2/'
lower_levels = '1,3/'
temperature = np.float64(7000.0)
line_flux_ratio = np.float64(1.506)
density = pyequib.calc_density(line_flux_ratio=line_flux_ratio, temperature=temperature, upper_levels=upper_levels, lower_levels=lower_levels, elj_data=s_ii_elj, omij_data=s_ii_omij, aij_data=s_ii_aij)
print("Electron Density:", density)

density = np.float64(1000)
temperature = np.float64(10000.0)
nlj = pyequib.calc_populations(temperature=temperature, density=density, elj_data=s_ii_elj, omij_data=s_ii_omij, aij_data=s_ii_aij)
print('Atomic Level Populations:', nlj)

temperature = np.float64(10000.0)
n_crit = pyequib.calc_crit_density(temperature=temperature, elj_data=s_ii_elj, omij_data=s_ii_omij, aij_data=s_ii_aij)
print('Critical Densities:', n_crit)

temperature = np.float64(10000.0)
omij_t = pyequib.get_omij_temp(temperature=temperature, omij_data=s_ii_omij)
print('Effective Collision Strengths: ')
print(omij_t)

pyequib.print_ionic(temperature=temperature, density=density,
            elj_data=s_ii_elj, omij_data=s_ii_omij, aij_data=s_ii_aij,
            h_i_aeff_data=hi_rc_data.aeff)

# --- End MAIN program. ---------------
