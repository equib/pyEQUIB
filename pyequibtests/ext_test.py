
"""Tests for pyequib.ext."""

import pyequib

wavelength=6563.0
m_ext=1.0
flux=1.0
fl=pyequib.ext.redlaw_gal(wavelength)
fl=pyequib.ext.redlaw_gal2(wavelength)
fl=pyequib.ext.redlaw_ccm(wavelength)
fl=pyequib.ext.redlaw_jbk(wavelength)
fl=pyequib.ext.redlaw_fm(wavelength)
fl=pyequib.ext.redlaw_smc(wavelength)
fl=pyequib.ext.redlaw_lmc(wavelength)
fl=pyequib.ext.redlaw(wavelength)
flux_deredden=pyequib.ext.deredden_relflux(wavelength, flux, m_ext)
flux_deredden=pyequib.ext.deredden_flux(wavelength, flux, m_ext)

wavelength=1250.0
m_ext=1.0
flux=1.0
fl=pyequib.ext.redlaw_gal(wavelength)
fl=pyequib.ext.redlaw_gal2(wavelength)
fl=pyequib.ext.redlaw_ccm(wavelength)
fl=pyequib.ext.redlaw_jbk(wavelength)
# fl=pyequib.ext.redlaw_fm(wavelength)
# fl=pyequib.ext.redlaw_smc(wavelength)
fl=pyequib.ext.redlaw_lmc(wavelength)
fl=pyequib.ext.redlaw(wavelength)
flux_deredden=pyequib.ext.deredden_relflux(wavelength, flux, m_ext)
flux_deredden=pyequib.ext.deredden_flux(wavelength, flux, m_ext)

wavelength=2000.0
m_ext=1.0
flux=1.0
fl=pyequib.ext.redlaw_gal(wavelength)
fl=pyequib.ext.redlaw_gal2(wavelength)
fl=pyequib.ext.redlaw_ccm(wavelength)
fl=pyequib.ext.redlaw_jbk(wavelength)
# fl=pyequib.ext.redlaw_fm(wavelength)
fl=pyequib.ext.redlaw_smc(wavelength)
fl=pyequib.ext.redlaw_lmc(wavelength)
fl=pyequib.ext.redlaw(wavelength)
flux_deredden=pyequib.ext.deredden_relflux(wavelength, flux, m_ext)
flux_deredden=pyequib.ext.deredden_flux(wavelength, flux, m_ext)

wavelength=3333.0
m_ext=1.0
flux=1.0
fl=pyequib.ext.redlaw_gal(wavelength)
fl=pyequib.ext.redlaw_gal2(wavelength)
fl=pyequib.ext.redlaw_ccm(wavelength)
fl=pyequib.ext.redlaw_jbk(wavelength)
fl=pyequib.ext.redlaw_fm(wavelength)
fl=pyequib.ext.redlaw_smc(wavelength)
fl=pyequib.ext.redlaw_lmc(wavelength)
fl=pyequib.ext.redlaw(wavelength)
flux_deredden=pyequib.ext.deredden_relflux(wavelength, flux, m_ext)
flux_deredden=pyequib.ext.deredden_flux(wavelength, flux, m_ext)

wavelength=4000.0
m_ext=1.0
flux=1.0
fl=pyequib.ext.redlaw_gal(wavelength)
fl=pyequib.ext.redlaw_gal2(wavelength)
fl=pyequib.ext.redlaw_ccm(wavelength)
fl=pyequib.ext.redlaw_jbk(wavelength)
fl=pyequib.ext.redlaw_fm(wavelength)
fl=pyequib.ext.redlaw_smc(wavelength)
fl=pyequib.ext.redlaw_lmc(wavelength)
fl=pyequib.ext.redlaw(wavelength)
flux_deredden=pyequib.ext.deredden_relflux(wavelength, flux, m_ext)
flux_deredden=pyequib.ext.deredden_flux(wavelength, flux, m_ext)
