# Example: calc_abund_he_i_rl(), calc_abund_he_ii_rl()
#          calc_abund_c_ii_rl(), calc_abund_c_iii_rl()
#          calc_abund_n_ii_rl(), calc_abund_n_iii_rl()
#          calc_abund_o_ii_rl(), calc_abund_ne_ii_rl()
#     determine ionic abundance from observed
#     flux intensity for gievn electron density
#     and temperature using calc_abund_he_i_rl,
#     calc_abund_he_ii_rl, calc_abund_c_ii_rl, calc_abund_c_iii_rl
#     calc_abund_n_ii_rl, calc_abund_n_iii_rl
#     calc_abund_o_ii_rl, and calc_abund_ne_ii_rl from proEQUIB
#
# --- Begin MAIN program. ---------------
#
#
import pyequib
import atomneb
import os
import numpy as np

# Locate datasets
base_dir = '../externals/atomneb/'
data_rc_dir = os.path.join('atomic-data-rc')
atom_rc_all_file = os.path.join(base_dir,data_rc_dir, 'rc_collection.fits')
atom_rc_he_i_file = os.path.join(base_dir,data_rc_dir, 'rc_he_ii_PFSD12.fits')
atom_rc_ppb91_file = os.path.join(base_dir,data_rc_dir, 'rc_PPB91.fits')
atom_rc_sh95_file = os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')

atom = 'h'
ion = 'ii' # H I
h_i_rc_data = atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)

atom = 'he'
ion = 'ii' # He I
he_i_rc_data = atomneb.read_aeff_he_i_pfsd12(atom_rc_he_i_file, atom, ion)

atom = 'he'
ion = 'iii' # He II
he_ii_rc_data = atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)

atom = 'c'
ion = 'iii' # C II
c_ii_rc_data = atomneb.read_aeff_collection(atom_rc_all_file, atom, ion)

atom = 'c'
ion = 'iv' # C III
c_iii_rc_data = atomneb.read_aeff_ppb91(atom_rc_ppb91_file, atom, ion)

atom = 'n'
ion = 'iii' # N II
n_ii_rc_data = atomneb.read_aeff_collection(atom_rc_all_file, atom, ion)
n_ii_rc_data_br = atomneb.read_aeff_collection(atom_rc_all_file, atom, ion, br=True)

atom = 'n'
ion = 'iv' # N III
n_iii_rc_data = atomneb.read_aeff_ppb91(atom_rc_ppb91_file, atom, ion)

atom = 'o'
ion = 'iii' # O II
o_ii_rc_data = atomneb.read_aeff_collection(atom_rc_all_file, atom, ion)
o_ii_rc_data_br = atomneb.read_aeff_collection(atom_rc_all_file, atom, ion, br=True)

atom = 'ne'
ion = 'iii' # Ne II
ne_ii_rc_data = atomneb.read_aeff_collection(atom_rc_all_file, atom, ion)

h_i_aeff_data = h_i_rc_data['aeff'][0]
he_i_aeff_data = he_i_rc_data['aeff'][0]
he_ii_aeff_data = he_ii_rc_data['aeff'][0]

temperature = np.float64(10000.0)
density = np.float64(5000.0)

# 4120.84: linenum=7
# 4387.93: linenum=8
# 4437.55: linenum=9
# 4471.50: linenum=10
# 4921.93: linenum=12
# 5015.68: linenum=13
# 5047.74: linenum=14
# 5875.66: linenum=15
# 6678.16: linenum=16
# 7065.25: linenum=17
# 7281.35: linenum=18
linenum = 10# 4471.50
emiss_he_i = pyequib.calc_emiss_he_i_rl(temperature=temperature, density=density,
                                linenum=linenum, he_i_aeff_data=he_i_aeff_data)
print('He I Emissivity:', emiss_he_i)
he_i_4471_flux = 2.104
abund_he_i = pyequib.calc_abund_he_i_rl(temperature=temperature, density=density,
                                linenum=linenum, line_flux=he_i_4471_flux,
                                he_i_aeff_data=he_i_aeff_data, h_i_aeff_data=h_i_aeff_data)
print('N(He^+)/N(H^+):', abund_he_i)

emiss_he_ii = pyequib.calc_emiss_he_ii_rl(temperature=temperature, density=density,
                                  he_ii_aeff_data=he_ii_aeff_data)
print('He II Emissivity:', emiss_he_ii)
he_ii_4686_flux = 135.833
abund_he_ii = pyequib.calc_abund_he_ii_rl(temperature=temperature, density=density,
                                  line_flux=he_ii_4686_flux,
                                  he_ii_aeff_data=he_ii_aeff_data, h_i_aeff_data=h_i_aeff_data)
print('N(He^2+)/N(H^+):', abund_he_ii)

wavelength = 6151.43
emiss_c_ii = pyequib.calc_emiss_c_ii_rl(temperature=temperature, density=density,
                                wavelength=wavelength, c_ii_rc_data=c_ii_rc_data)
print('C II Emissivity:', emiss_c_ii)
c_ii_6151_flux = 0.028
abund_c_ii = pyequib.calc_abund_c_ii_rl(temperature=temperature, density=density,
                                wavelength=wavelength, line_flux=c_ii_6151_flux,
                                c_ii_rc_data=c_ii_rc_data, h_i_aeff_data=h_i_aeff_data)
print('N(C^2+)/N(H+):', abund_c_ii)

wavelength = 4647.42
emiss_c_iii = pyequib.calc_emiss_c_iii_rl(temperature=temperature, density=density,
                                  wavelength=wavelength,
                                  c_iii_rc_data=c_iii_rc_data)
print('C III Emissivity:', emiss_c_iii)
c_iii_4647_flux = 0.107
abund_c_iii = pyequib.calc_abund_c_iii_rl(temperature=temperature, density=density,
                                  wavelength=wavelength,
                                  line_flux=c_iii_4647_flux, c_iii_rc_data=c_iii_rc_data,
                                  h_i_aeff_data=h_i_aeff_data)
print('N(C^3+)/N(H+):', abund_c_iii)

wavelength = 4442.02
emiss_n_ii = pyequib.calc_emiss_n_ii_rl(temperature=temperature, density=density,
                                wavelength=wavelength,
                                n_ii_rc_br=n_ii_rc_data_br, n_ii_rc_data=n_ii_rc_data)
print('N II Emissivity:', emiss_n_ii)
n_ii_4442_flux = 0.017
abund_n_ii = pyequib.calc_abund_n_ii_rl(temperature=temperature, density=density,
                                wavelength=wavelength, line_flux=n_ii_4442_flux,
                                n_ii_rc_br=n_ii_rc_data_br, n_ii_rc_data=n_ii_rc_data,
                                h_i_aeff_data=h_i_aeff_data)
print('N(N^2+)/N(H+):', abund_n_ii)


wavelength = 4640.64
emiss_n_iii = pyequib.calc_emiss_n_iii_rl(temperature=temperature, density=density,
                                  wavelength=wavelength, n_iii_rc_data=n_iii_rc_data)
print('N III Emissivity:', emiss_n_iii)
n_iii_4641_flux = 0.245
abund_n_iii = pyequib.calc_abund_n_iii_rl(temperature=temperature, density=density,
                                  wavelength=wavelength, line_flux=n_iii_4641_flux,
                                  n_iii_rc_data=n_iii_rc_data, h_i_aeff_data=h_i_aeff_data)
print('N(N^3+)/N(H+):', abund_n_iii)


wavelength = 4613.68
emiss_o_ii = pyequib.calc_emiss_o_ii_rl(temperature=temperature, density=density,
                                wavelength=wavelength,
                                o_ii_rc_br=o_ii_rc_data_br, o_ii_rc_data=o_ii_rc_data)
print('O II Emissivity:', emiss_o_ii)
o_ii_4614_flux = 0.009
abund_o_ii = pyequib.calc_abund_o_ii_rl(temperature=temperature, density=density,
                                wavelength=wavelength, line_flux=o_ii_4614_flux,
                                o_ii_rc_br=o_ii_rc_data_br,
                                o_ii_rc_data=o_ii_rc_data,
                                h_i_aeff_data=h_i_aeff_data)
print('N(O^2+)/N(H+):', abund_o_ii)

wavelength = 3777.14
emiss_ne_ii = pyequib.calc_emiss_ne_ii_rl(temperature=temperature, density=density,
                                  wavelength=wavelength, ne_ii_rc_data=ne_ii_rc_data)
print('Ne II Emissivity:', emiss_ne_ii)
ne_ii_3777_flux = 0.056
abund_ne_ii = pyequib.calc_abund_ne_ii_rl(temperature=temperature, density=density,
                                  wavelength=wavelength, line_flux=ne_ii_3777_flux,
                                  ne_ii_rc_data=ne_ii_rc_data, h_i_aeff_data=h_i_aeff_data)
print('N(Ne^2+)/N(H+):', abund_ne_ii)

# --- End MAIN program. ---------------

