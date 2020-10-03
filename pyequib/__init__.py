"""
pyequib - Python Package for Plasma Diagnostics and Abundance Analysis
"""

__all__ =  ["calc_temperature", "calc_density",
            "calc_populations", "calc_crit_density",
            "calc_emissivity", "calc_abundance",
            "print_ionic", "get_omij_temp",
            "calc_emiss_h_beta",
            "calc_emiss_he_i_rl", "calc_emiss_he_ii_rl",
            "calc_emiss_c_ii_rl", "calc_emiss_c_iii_rl",
            "calc_emiss_n_ii_rl", "calc_emiss_n_iii_rl",
            "calc_emiss_o_ii_rl", "calc_emiss_ne_ii_rl",
            "calc_abund_he_i_rl", "calc_abund_he_ii_rl",
            "calc_abund_c_ii_rl", "calc_abund_c_iii_rl",
            "calc_abund_n_ii_rl", "calc_abund_n_iii_rl",
            "calc_abund_o_ii_rl", "calc_abund_ne_ii_rl",
            "redlaw", "redlaw_gal", "redlaw_gal2",
            "redlaw_ccm", "redlaw_jbk", "redlaw_fm",
            "redlaw_smc", "redlaw_lmc",
            "deredden_flux", "deredden_relflux"]

from .pyequib import calc_temperature, calc_density, \
             calc_populations, calc_crit_density, \
             calc_emissivity, calc_abundance, \
             print_ionic, get_omij_temp, \
             calc_emiss_h_beta, \
             calc_emiss_he_i_rl, calc_emiss_he_ii_rl, \
             calc_emiss_c_ii_rl, calc_emiss_c_iii_rl, \
             calc_emiss_n_ii_rl, calc_emiss_n_iii_rl, \
             calc_emiss_o_ii_rl, calc_emiss_ne_ii_rl, \
             calc_abund_he_i_rl, calc_abund_he_ii_rl, \
             calc_abund_c_ii_rl, calc_abund_c_iii_rl, \
             calc_abund_n_ii_rl, calc_abund_n_iii_rl, \
             calc_abund_o_ii_rl, calc_abund_ne_ii_rl, \
             redlaw, redlaw_gal, redlaw_gal2, \
             redlaw_ccm, redlaw_jbk, redlaw_fm, \
             redlaw_smc, redlaw_lmc, \
             deredden_flux, deredden_relflux

from .version import __version__

import sys
from numpy.version import version as numpy_version

if sys.version_info[0:2] < (2, 6):
    log_.warn('pyequib requires Python version >= 2.6, but it is version {0}'.format(sys.version_info), calling='pyequib')
try:
    if [int(n) for n in (numpy_version.split('.')[:3])] < [1, 5, 1] :
        log_.warn('pyequib Numpy version >= 1.5.1, but it is version {0}'.format(numpy_version), calling='pyequib')
except:
    log_.warn('Cannot find Numpy version {0}, report the bug'.format(numpy_version), calling='pyemcee')
    


