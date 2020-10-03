# Example: redlaw() and deredden_flux ()
#     determine reddening law function for
#     given wavelength and deredden flux
#     intensity (relative to Hb=100) based on
#     the reddening law
#
# --- Begin MAIN program. ---------------
#
#
import pyequib
import numpy as np

wavelength = 6563.0
m_ext = 1.0
flux = 1.0
r_v = 3.1

fl = pyequib.redlaw_gal(wavelength, rv=r_v)
print('fl(6563):', fl)

fl = pyequib.redlaw_gal2(wavelength)
print('fl(6563):', fl)

fl = pyequib.redlaw_ccm(wavelength, rv=r_v)
print('fl(6563):', fl)

fl = pyequib.redlaw_jbk(wavelength)
print('fl(6563):', fl)

#fmlaw='LMC2'
fmlaw = 'AVGLMC'
fl = pyequib.redlaw_fm(wavelength, fmlaw=fmlaw, rv=r_v)
print('fl(6563):', fl)

fl = pyequib.redlaw_smc(wavelength)
print('fl(6563):', fl)

fl = pyequib.redlaw_lmc(wavelength)
print('fl(6563):', fl)

fl = pyequib.redlaw(wavelength, rv=r_v)
print('fl(6563):', fl)

ext_law = 'GAL'
r_v = 3.1
# deredden flux intensity relative to Hb=100
flux_deredden = pyequib.deredden_relflux(wavelength, flux, m_ext, ext_law=ext_law, rv=r_v)
print('dereddened flux(6563):', flux_deredden)

# deredden absolute flux intensity
flux_deredden = pyequib.deredden_flux(wavelength, flux, m_ext, ext_law=ext_law, rv=r_v)
print('dereddened flux(6563):', flux_deredden)

# --- End MAIN program. ---------------
