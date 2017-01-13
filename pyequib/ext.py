"""
This module contains functions for extinction corrections using different laws
"""

# A. Danehkar
#
# Version 0.1, 15/08/2016
# First Release
#

import numpy, os
import array, math
from scipy import interpolate

"""class _Extinction(object):
    
    def __init__(self):
        pass
""" 

def deredden_relflux(wavelength, relflux, m_ext, ext_law=None, rv=None, fmlaw=None):
   """
    NAME:
        deredden_relflux
    PURPOSE:
        determine deredden flux intensity relative to Hb=100,
        depending on the law used
    EXPLANATION:
   
    CALLING SEQUENCE:
        relflux_deredden = deredden_relflux(wavelengths, relflux, m_ext)
   
    INPUTS:
        wavelength -     Wavelength in Angstrom
        relflux -      flux intensity relative to Hb=100
        m_ext -      logarithmic extinction
        ext_law -     extinction law
           ext_law='GAL' ; Howarth Galactic
           ext_law='GAL2' ; Savage and Mathis
           ext_law='CCM' ; CCM galactic
           ext_law='JBK' ; Whitford, Seaton, Kaler
           ext_law='FM' ; Fitxpatrick
           ext_law='SMC' ; Prevot SMC
           ext_law='LMC' ; Howarth LMC
    RETURN:  deredden relative intensity
   
    REVISION HISTORY:
        Python code by A. Danehkar, 31/08/2012
   """
   if (ext_law is not None):   
      extlaw = ext_law
   else:   
      extlaw = 'GAL'
   if (rv is not None):   
      r_v = rv
   else:   
      r_v = 3.1
   if (fmlaw is not None):   
      fm_law = fmlaw
   else:   
      fm_law = 'STANDARD'
   fl = redlaw(wavelength, ext_law=extlaw, rv=r_v, fmlaw=fm_law)
   int_dered = relflux * 10.0 ** (m_ext * fl)
   return int_dered

def deredden_flux(wavelength, flux, m_ext, ext_law=None, rv=None, fmlaw=None):
   """
    NAME:
        deredden
    PURPOSE:
        determine deredden flux intensity relative to Hb=100,
        depending on the law used
    EXPLANATION:
   
    CALLING SEQUENCE:
        flux_deredden = deredden_flux(wavelengths, flux, m_ext)
   
    INPUTS:
        wavelength -     Wavelength in Angstrom
        flux -      absolute flux intensity
        m_ext -      logarithmic extinction
        ext_law -     extinction law
           ext_law='GAL' ; Howarth Galactic
           ext_law='GAL2' ; Savage and Mathis
           ext_law='CCM' ; CCM galactic
           ext_law='JBK' ; Whitford, Seaton, Kaler
           ext_law='FM' ; Fitxpatrick
           ext_law='SMC' ; Prevot SMC
           ext_law='LMC' ; Howarth LMC
    RETURN:  deredden relative intensity
   
    REVISION HISTORY:
        Python code by A. Danehkar, 31/08/2012
   """
   if (ext_law is not None):   
      extlaw = ext_law
   else:   
      extlaw = 'GAL'
   if (rv is not None):   
      r_v = rv
   else:   
      r_v = 3.1
   if (fmlaw is not None):   
      fm_law = fmlaw
   else:   
      fm_law = 'STANDARD'
   fl = redlaw(wavelength, ext_law=extlaw, rv=r_v, fmlaw=fm_law)
   int_dered = flux * 10.0 ** (m_ext * (1 + fl))
   return int_dered

def redlaw(wavelength, ext_law=None, rv=None, fmlaw=None):
   """
    NAME:
        redlaw
    PURPOSE:
        determine reddening law function for the line at the wavelength of lambda,
        depending on the law used
    EXPLANATION:
   
    CALLING SEQUENCE:
        fl = redlaw(wavelength, ext_law)
   
    INPUTS:
        wavelength -     Wavelength in Angstrom
        ext_law -     extinction law
           ext_law='GAL' ; Howarth Galactic
           ext_law='GAL2' ; Savage and Mathis
           ext_law='CCM' ; CCM galactic
           ext_law='JBK' ; Whitford, Seaton, Kaler
           ext_law='FM' ; Fitxpatrick
           ext_law='SMC' ; Prevot SMC
           ext_law='LMC' ; Howarth LMC
    RETURN:  f(lambda)
   
    REVISION HISTORY:
        Originally from IRAF STSDAS SYNPHOT redlaw.x, ebmvxfunc.x
        Converted to Python code by A. Danehkar, 31/08/2012
   """
   
   if (ext_law is not None):   
      extlaw = ext_law
   else:   
      extlaw = 'GAL'
   _expr = extlaw
   if _expr == 'GAL':   
      fl = redlaw_gal(wavelength, rv=rv)
   elif _expr == 'GAL2':   
      fl = redlaw_gal2(wavelength)
   elif _expr == 'CCM':   
      fl = redlaw_ccm(wavelength, rv=rv)
   elif _expr == 'JBK':  
      fl = redlaw_jbk(wavelegth)
   elif _expr == 'FM':   
      fl = redlaw_fm(wavelength, fmlaw=fmlaw, rv=rv)
   elif _expr == 'SMC':   
      fl = redlaw_smc(wavelength)
   elif _expr == 'LMC':   
      fl = redlaw_lmc(wavelength)
   else:   
      print('ext_law cannnot find')
   
   return fl

def redlaw_gal(wave, rv=None):
   """
   NAME:
       redlaw_smc
   PURPOSE:
      reddening law function for Galactic Seaton1979+Howarth1983+CCM1983
   
   EXPLANATION:
   
   CALLING SEQUENCE:
       fl = redlaw_gal(wave)
   
   INPUTS:
       wave[] -  wavelength of emission line, Angstroms
   RETURN: extl[] -  extinction evaluation array
   
   REVISION HISTORY:
       Based on the UV Formulae from Seaton 1979, MNRAS, 187, 73
       1979MNRAS.187P..73S, the opt/NIR from Howarth 1983, MNRAS, 203, 301
       the FIR from Cardelli, Clayton and Mathis 1989, ApJ, 345, 245
       1989ApJ...345..245C
       Originally from IRAF STSDAS SYNPHOT ebmvxfunc.x, pyneb.extinction
       Converted to Python code by A. Danehkar, 31/08/2012
   """
   # Tabulated inverse wavelengths in microns:
   xtable = numpy.array([0., 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7])
   etable = numpy.array([0., 1.36, 1.64, 1.84, 2.04, 2.24, 2.44, 2.66, 2.88, 3.14, 3.36, 3.56, 3.77, 3.96, 4.15, 4.26, 4.40, 4.52, 4.64])
   
   if hasattr(wave, "__len__"):
     npts = len(wave)
     extl = numpy.zeros(npts)
   else:
     npts = 1
     extl = numpy.int32(0)
   if (rv is not None):   
      r_v = rv
   else:   
      r_v = 3.1
   for pix in range(0, (npts - 1)+(1)):
   # Convert wavelength in angstroms to 1/microns
      if hasattr(wave, "__len__"):
         wavel=wave[pix]
      else:
         wavel=wave
      x = 10000.e+0 / wavel
      
      if (x <= 1.1):   
         # Infrared: extend optical results linearly to 0 at 1/lam = 0
         # extl[pix] = etable[2] * x^2
         # Cardelli, Clayton and Mathis 1989, ApJ, 345, 245
         val = r_v * (0.574 * x ** 1.61) - 0.527 * x ** 1.61
      else:   
         if (x < 1.83):   
            # Howarth 1983, Galactic
            val = (r_v - 3.1) + (((1.86 * x ** 2.0) - (0.48 * x ** 3.0) - (0.1 * x)))
         else:   
            if (x < 2.75):   
               # Optical region interpolates in Seaton's table 3
               #extl[pix] = lin_interp(etable, xtable,  x)
               # Howarth 1983, MNRAS, 203, 301
               val = r_v + 2.56 * (x - 1.83) - 0.993 * (x - 1.83) ** 2.0 #
            else:   
               if (x < 3.65):   
                  # Ultraviolet uses analytic formulae from Seaton's table 2
                  val = (r_v - 3.2) + 1.56 + 1.048 * x + 1.01 / ((x - 4.6) ** 2 + 0.280)
               else:   
                  if (x < 7.14):   
                     # More ultraviolet
                     val = (r_v - 3.2) + 2.29 + 0.848 * x + 1.01 / ((x - 4.6) ** 2 + 0.280)
                     # And more ultraviolet still
                  else:   
                     x = min([x, 50.0])
                     val = (r_v - 3.2) + 16.17 - 3.20 * x + 0.2975 * x ** 2
      if hasattr(extl, "__len__"):
         extl[pix] = val
      else:
         extl=val
   return (extl / 3.63) - 1.0

def redlaw_gal2(wave):
   """
   NAME:
       redlaw_smc
   PURPOSE:
      reddening law function for Galactic Savage & Mathis 1979
   
   EXPLANATION:
   
   CALLING SEQUENCE:
       fl = redlaw_gal2(wave)
   
   INPUTS:
       wave[] -  wavelength of emission line, Angstroms
   RETURN: extl[] -  extinction evaluation array
   
   REVISION HISTORY:
       Based on Savage & Mathis 1979, ARA&A, vol. 17, 73-111
       Originally from IRAF STSDAS SYNPHOT ebmvxfunc.x
       Initial IRAF implementation, R. A. Shaw, 20/09/1994
       Return A(lambda)/A(V) instead, R. A. Shaw, 04/03/95
       Converted to Python code by A. Danehkar, 31/08/2012
   """
   # Tabulated inverse wavelengths in microns:
   xtable = numpy.array([0.00, 0.29, 0.45, 0.80, 1.11, 1.43, 1.82, 2.27, 2.50, 2.91, 3.65, 4.00, 4.17, 4.35, 4.57, 4.76, 5.00, 5.26, 5.56, 5.88, 6.25, 6.71, 7.18, 8.00, 8.50, 9.00, 9.50, 10.00])
   
   #  Tabulated extinction function, A(lambda)/E(B-V):
   etable = numpy.array([0.00, 0.16, 0.38, 0.87, 1.50, 2.32, 3.10, 4.10, 4.40, 4.90, 6.20, 7.29, 8.00, 8.87, 9.67, 9.33, 8.62, 8.00, 7.75, 7.87, 8.12, 8.15, 8.49, 9.65, 10.55, 11.55, 12.90, 14.40])
   
   if hasattr(wave, "__len__"):
     npts = len(wave)
     extl = numpy.zeros(npts)
   else:
     npts = 1
     extl = numpy.int32(0)
   for pix in range(0, (npts - 1)+(1)):
      if hasattr(wave, "__len__"):
         wavel=wave[pix]
      else:
         wavel=wave
      if (wavel < 1000.0):
         print('redlaw_smc: Invalid wavelength')
      # Convert wavelength in angstroms to 1/microns
      x = 10000.e+0 / wavel
      x = min([x, 10.0])
      
      # Linearly interpolate extinction law in 1/lam
      val = lin_interp(etable, xtable, x)
      #deriv1 = spl_init(xtab, extab)
      #val=spl_interp(xtab, extab, deriv1, x)
      if hasattr(extl, "__len__"):
         extl[pix] = val
      else:
         extl=val
   return (extl / 3.63) - 1.0

def redlaw_ccm(wave, rv=None):
   """
    NAME:
        redlaw_ccm
    PURPOSE:
       reddening law function of Cardelli, Clayton & Mathis
   
    EXPLANATION:
   
    CALLING SEQUENCE:
        fl = redlaw_ccm(wave, rv)
   
    INPUTS:
        wave[] -  wavelength of emission line, Angstroms
        rv - The ratio of extinction to reddening defined as R_V = A_V/E(B-V)
    RETURN: extl[] -  extinction evaluation array
   
    REVISION HISTORY:
        Based on Formulae by Cardelli, Clayton & Mathis 1989, ApJ 345, 245-256.
        1989ApJ...345..245C
        Originally from IRAF STSDAS SYNPHOT redlaw.x
        Initial IRAF implementation, based upon CCM module
            in onedspec.deredden, R. A. Shaw, 18/05/1993
        Converted to Python code by A. Danehkar, 31/08/2012
   """
   if hasattr(wave, "__len__"):
     npts = len(wave)
     extl = numpy.zeros(npts)
   else:
     npts = 1
     extl = numpy.int32(0)
   if (rv is not None):   
      r_v = rv
   else:   
      r_v = 3.1
   for pix in range(0, (npts - 1)+(1)):
      if hasattr(wave, "__len__"):
         wavel=wave[pix]
      else:
         wavel=wave
      if (wavel < 1000.0):   
         print('redlaw_ccm: Invalid wavelength')
      
      # Convert input wavelength to inverse microns
      x = 10000.e+0 / wavel
      
      # For wavelengths longward of the L passband, linearly interpolate
      # to 0. at 1/lambda = 0.  (a=0.08, b=-0.0734 @ x=0.29)
      if (x < 0.29e+0):   
         a = 0.2759 * x
         b = -0.2531 * x
      else:   
         if (x < 1.1e+0):   
            y = x ** 1.61
            a = 0.574 * y
            b = -0.527 * y
         else:   
            if (x < 3.3e+0):   
               y = x - 1.82
               a = 1 + y * (0.17699 + y * (-0.50447 + y * (-0.02427 + y * (0.72085 + y * (0.01979 + y * (-0.77530 + y * 0.32999))))))
               b = y * (1.41338 + y * (2.28305 + y * (1.07233 + y * (-5.38434 + y * (-0.62251 + y * (5.30260 - y * 2.09002))))))
            else:   
               if (x < 5.9e+0):   
                  y = (x - 4.67) ** 2
                  a = 1.752 - 0.316 * x - 0.104 / (y + 0.341)
                  b = -3.090 + 1.825 * x + 1.206 / (y + 0.263)
               else:   
                  if (x < 8.0e+0):   
                     y = (x - 4.67) ** 2
                     a = 1.752 - 0.316 * x - 0.104 / (y + 0.341)
                     b = -3.090 + 1.825 * x + 1.206 / (y + 0.263)
                     
                     y = x - 5.9
                     a = a - 0.04473 * y ** 2 - 0.009779 * y ** 3
                     b = b + 0.2130 * y ** 2 + 0.1207 * y ** 3
                     # Truncate range of ISEF to that for 1000 Ang.
                  else:   
                     if (x <= 10.0e+0):   
                        x = min([x, 10.0e+0])
                        y = x - 8.e+0
                        a = -1.072 - 0.628 * y + 0.137 * y ** 2 - 0.070 * y ** 3
                        b = 13.670 + 4.257 * y - 0.420 * y ** 2 + 0.374 * y ** 3
      # Compute A(lambda)/A(V)
      y = a * r_v + b
      if hasattr(extl, "__len__"):
         extl[pix] = y
      else:
         extl= y
   return (extl / ((1.015452 * r_v) + 0.461000)) - 1.0


def redlaw_jbk(wave):
   """
   NAME:
       redlaw_smc
   PURPOSE:
      reddening law function for Galactic Whitford1958 + Seaton1977 + Kaler1976
   
   EXPLANATION:
   
   CALLING SEQUENCE:
        fl = redlaw_jbk(wave)
   
   INPUTS:
       wave[] -  wavelength of emission line, Angstroms
   RETURN: extl[] -  extinction evaluation array
   
   REVISION HISTORY:
       Based on Whitford (1958), extended to the UV by Seaton (1977),
       adapted by Kaler (1976).
       Originally from IRAF STSDAS SYNPHOT redlaw.x
       Initial IRAF implementation, R. A. Shaw, 13/05/1993
       Converted to Python code by A. Danehkar, 31/08/2012
   """
   # Tabulated wavelengths, Angstroms:
   refw = numpy.array([1150., 1200., 1250., 1300., 1350., 1400., 1450., 1500., 1550., 1600., 1650., 1700., 1750., 1800., 1850., 1900., 1950., 2000., 2050., 2100., 2150., 2200., 2250., 2300., 2350., 2400., 2450., 2500., 2550., 2600., 2650., 2700., 2750., 2800., 2850., 2900., 2950., 3000., 3050., 3100., 3333., 3500., 3600., 3700., 3800., 3900., 4000., 4100., 4200., 4300., 4400., 4500., 4600., 4700., 4800., 4861.3, 5000., 5100., 5200., 5300., 5400., 5500., 5600., 5700., 5800., 5900., 6000., 6100., 6200., 6300., 6400., 6500., 6600., 6700., 6800., 6900., 7000., 7200., 7400., 7600., 7800., 8000., 8200., 8400., 8600., 8800., 9000., 9500., 10000., 11000., 12000., 14000., 16000., 20000., 1.e+6])
   
   #  Tabulated extinction function:
   extab = numpy.array([1.96, 1.78, 1.61, 1.49, 1.37, 1.29, 1.24, 1.20, 1.20, 1.20, 1.17, 1.13, 1.11, 1.10, 1.12, 1.17, 1.25, 1.35, 1.45, 1.53, 1.60, 1.62, 1.52, 1.40, 1.28, 1.17, 1.06, 0.98, 0.9, 0.84, 0.77, 0.72, 0.68, 0.64, 0.60, 0.57, 0.53, 0.51, 0.48, 0.46, 0.385, 0.358, 0.33, 0.306, 0.278, 0.248, 0.220, 0.195, 0.168, 0.143, 0.118, 0.095, 0.065, 0.040, 0.015, 0.000, -0.030, -0.055, -0.078, -0.10, -0.121, -0.142, -0.164, -0.182, -0.201, -0.220, -0.238, -0.254, -0.273, -0.291, -0.306, -0.321, -0.337, -0.351, -0.365, -0.377, -0.391, -0.416, -0.441, -0.465, -0.490, -0.510, -0.529, -0.548, -0.566, -0.582, -0.597, -0.633, -0.663, -0.718, -0.763, -0.840, -0.890, -0.960, -1.000])
   
   xtable = 10000.e+0 / refw
   if hasattr(wave, "__len__"):
     npts = len(wave)
     extl = numpy.zeros(npts)
   else:
     npts = 1
     extl = numpy.int32(0)
   for pix in range(0, (npts - 1)+(1)):
      if hasattr(wave, "__len__"):
         wavel=wave[pix]
      else:
         wavel=wave
      if (wavel < 1000.0):   
         print('redlaw_smc: Invalid wavelength')
      # Convert wavelength in angstroms to 1/microns
      x = 10000.e+0 / wavel
      x = min([x, 10.0])
      
      # Linearly interpolate extinction law in 1/lam
      val = lin_interp(extab, xtable, x)
      #deriv1 = spl_init(xtab, extab)
      #val=spl_interp(xtab, extab, deriv1, x)
      if hasattr(extl, "__len__"):
         extl[pix] = val
      else:
         extl=val
   return extl

def redlaw_fm(wave, fmlaw=None, rv=None):
   """
   NAME:
       redlaw_fm
   PURPOSE:
      reddening law function for Small Magellanic Cloud
   
   EXPLANATION:
   
   CALLING SEQUENCE:
       fl = redlaw_smc(wave)
   
   INPUTS:
       wave[] -  wavelength of emission line, Angstroms
   RETURN: extl[] -  extinction evaluation array
   
   REVISION HISTORY:
       Based on Formulae by Fitzpatrick 1999, PASP, 11, 63
       1999PASP..111...63F, Fitzpatrick & Massa 1990,
       ApJS, 72, 163, 1990ApJS...72..163F
       Adopted from NASA IDL Library & PyAstronomy
       Revised in Python code by A. Danehkar, 30/12/2016
   """
   
   # Tabulated inverse wavelengths in microns:
   if hasattr(wave, "__len__"):
     npts = len(wave)
     extl = numpy.zeros(npts)
   else:
     npts = 1
     extl = numpy.int32(0)
   if (rv is not None):   
      r_v = rv
   else:   
      r_v = 3.1
   x0 = 4.596
   gamma1 = 0.99
   c3 = 3.23
   c4 = 0.41
   c2 = -0.824 + 4.717 / r_v
   c1 = 2.030 - 3.007 * c2
   if (fmlaw is not None):   
      if fmlaw == 'LMC2':   
         x0 = 4.626
         gamma1 = 1.05
         c4 = 0.42
         c3 = 1.92
         c2 = 1.31
         c1 = -2.16
      else:   
         if fmlaw == 'AVGLMC':   
            x0 = 4.596
            gamma1 = 0.91
            c4 = 0.64
            c3 = 2.73
            c2 = 1.11
            c1 = -1.28
         else:   
            x0 = 4.596
            gamma1 = 0.99
            c3 = 3.23
            c4 = 0.41
            c2 = -0.824 + 4.717 / r_v
            c1 = 2.030 - 3.007 * c2
   for pix in range(0, (npts - 1)+(1)):
   # Convert input wavelength to inverse microns
      if hasattr(wave, "__len__"):
         wavel=wave[pix]
      else:
         wavel=wave
      x = 10000.e+0 / wavel
      curve = x * 0.
      
      # Compute UV portion of A(lambda)/E(B-V) curve using FM fitting function and
      # R-dependent coefficients
      xcutuv = numpy.array([10000.0/2700.0])
      xspluv = 10000.0/numpy.array([2700.0,2600.0])
      
      iuv = numpy.where(x >= xcutuv)[0]
      n_uv = len(iuv)
      iopir = numpy.where(x <= xcutuv)[0]
      nopir = len(iopir)
      if (n_uv > 0):   
         xuv = numpy.concatenate((xspluv,x[iuv]))
      else:   
         xuv = xspluv
      
      yuv = c1 + c2 * xuv
      yuv = yuv + c3 * xuv ** 2 / ((xuv ** 2 - x0 ** 2) ** 2 + (xuv * gamma1) ** 2)
      yuv = yuv + c4 * (0.5392 * ((numpy.maximum(xuv, 5.9)) - 5.9) ** 2 + 0.05644 * ((numpy.maximum(xuv, 5.9)) - 5.9) ** 3)
      yuv = yuv + r_v
      yspluv = yuv[0:2]                  # save spline points
      
      if (n_uv > 0):   
         curve[iuv] = yuv[2:]      # remove spline points
      
      # Compute optical portion of A(lambda)/E(B-V) curve
      # using cubic spline anchored in UV, optical, and IR
      xsplopir = numpy.concatenate(([0],10000.0/numpy.array([26500.0,12200.0,6000.0,5470.0,4670.0,4110.0])))
      ysplir   = numpy.array([0.0,0.26469,0.82925])*r_v/3.1 
      ysplop   = numpy.array((numpy.polyval([-4.22809e-01, 1.00270, 2.13572e-04][::-1],r_v), 
            numpy.polyval([-5.13540e-02, 1.00216, -7.35778e-05][::-1],r_v), 
            numpy.polyval([ 7.00127e-01, 1.00184, -3.32598e-05][::-1],r_v), 
            numpy.polyval([ 1.19456, 1.01707, -5.46959e-03, 7.97809e-04, -4.45636e-05][::-1],r_v) ))
      
      ysplopir = numpy.concatenate((ysplir, ysplop))
      
      if (nopir > 0):   
         tck = interpolate.splrep(numpy.concatenate((xsplopir,xspluv)),numpy.concatenate((ysplopir,yspluv)),s=0)
         if hasattr(extl, "__len__"):
            curve[iopir] = interpolate.splev(x[iopir], tck)
         else:
            curve = interpolate.splev(x, tck)
      if hasattr(extl, "__len__"):
         extl[pix] = curve
      else:
         extl=curve
   return (extl / 3.63) - 1.0


def redlaw_smc(wave):
   """
   NAME:
       redlaw_smc
   PURPOSE:
      reddening law function for Small Magellanic Cloud
   
   EXPLANATION:
   
   CALLING SEQUENCE:
       fl = redlaw_smc(wave)
   
   INPUTS:
       wave[] -  wavelength of emission line, Angstroms
   RETURN: extl[] -  extinction evaluation array
   
   REVISION HISTORY:
       Based on Prevot et al. (1984), A&A, 132, 389-392
       1984A%26A...132..389P
       Originally from IRAF STSDAS SYNPHOT redlaw.x, ebmvxfunc.x
       Initial IRAF implementation, R. A. Shaw, 20/09/1994
       Return A(lambda)/A(V) instead, R. A. Shaw, 04/03/95
       Converted to Python code by A. Danehkar, 31/08/2012
   """
   # Tabulated inverse wavelengths in microns:
   xtab = numpy.array([0.00, 0.29, 0.45, 0.80, 1.11, 1.43, 1.82, 2.35, 2.70, 3.22, 3.34, 3.46, 3.60, 3.75, 3.92, 4.09, 4.28, 4.50, 4.73, 5.00, 5.24, 5.38, 5.52, 5.70, 5.88, 6.07, 6.27, 6.48, 6.72, 6.98, 7.23, 7.52, 7.84])
   
   # Tabulated extinction function, E(lambda-V)/E(B-V):
   extab = numpy.array([-3.10, -2.94, -2.72, -2.23, -1.60, -0.78, 0.00, 1.00, 1.67, 2.29, 2.65, 3.00, 3.15, 3.49, 3.91, 4.24, 4.53, 5.30, 5.85, 6.38, 6.76, 6.90, 7.17, 7.71, 8.01, 8.49, 9.06, 9.28, 9.84, 10.80, 11.51, 12.52, 13.54])
   
   if hasattr(wave, "__len__"):
     npts = len(wave)
     extl = numpy.zeros(npts)
   else:
     npts = 1
     extl = numpy.int32(0)
   for pix in range(0, (npts - 1)+(1)):
      if hasattr(wave, "__len__"):
         wavel=wave[pix]
      else:
         wavel=wave
      if (wavel < 1000.0):   
         print('redlaw_smc: Invalid wavelength')
      # Convert wavelength in angstroms to 1/microns
      x = 10000.e+0 / wavel
      x = min([x, 10.0])
      
      # Linearly interpolate extinction law in 1/lam
      val = lin_interp(extab, xtab, x)
      #deriv1 = spl_init(xtab, extab)
      #val=spl_interp(xtab, extab, deriv1, x)
      
      if hasattr(extl, "__len__"):
         extl[pix] = val + 3.1
      else:
         extl=val + 3.1
   return (extl / 3.242) - 1.0


def redlaw_lmc(wave):
   """
   NAME:
       redlaw_ccm
   PURPOSE:
      reddening law function for the Large Magellanic Cloud
   
   EXPLANATION:
   
   CALLING SEQUENCE:
       fl = redlaw_ccm(wave)
   
   INPUTS:
       wave[] -  wavelength of emission line, Angstroms
   RETURN: extl[] -  extinction evaluation array
   
   REVISION HISTORY:
       Based on Formulae by Howarth 1983, MNRAS, 203, 301
       1983MNRAS.203..301H
       Originally from IRAF STSDAS SYNPHOT ebmvlfunc.x, redlaw.x
       Initial IRAF implementation, R. A. Shaw, 18/10/1994
       Return A(lambda)/A(V) instead, R. A. Shaw, 14/03/95
       Converted to Python code by A. Danehkar, 31/08/2012
   """
   # Tabulated inverse wavelengths in microns:
   xtab = numpy.array([0.00, 0.29, 0.45, 0.80, 1.11, 1.43, 1.82])
   
   # Tabulated extinction function, A(lambda)/E(B-V), from Savage & Mathis:
   extab = numpy.array([0.00, 0.16, 0.38, 0.87, 1.50, 2.32, 3.10])
   
   if hasattr(wave, "__len__"):
     npts = len(wave)
     extl = numpy.zeros(npts)
   else:
     npts = 1
     extl = numpy.int32(0)
   for pix in range(0, (npts - 1)+(1)):
      if hasattr(wave, "__len__"):
         wavel=wave[pix]
      else:
         wavel=wave
      if (wavel < 1000.0):   
         print('redlaw_lmc: Invalid wavelength')
      
      # Convert input wavelength to inverse microns
      x = 10000.e+0 / wavel
      
      # Infrared - optical
      if (x <= 1.82):   
         # linear interpolation of Savage & Mathis 1979
         #  val = lin_interp(extab, xtab,  x)
         #  extl[pix] = val ;+ 3.1
         # Infrared - extend optical results linearly to 0 at 1/lam = 0
         val = ((1.86 - 0.48 * x) * x - 0.1) * x
      else:   
         # The following polynomial evaluations assume R = 3.1
         # Renormalize extinction function to A(lambda)/A(V)
         if (x <= 2.75):   
            #  Violet
            val = 3.1 + (2.04 + 0.094 * (x - 1.82)) * (x - 1.82)
         else:   
            #  Ultraviolet out to lambda = 1000 A
            x = min([x, 10.0])
            val = 3.1 - 0.236 + 0.462 * x + 0.105 * x ** 2 + 0.454 / ((x - 4.557) ** 2 + 0.293)
      if hasattr(extl, "__len__"):
         extl[pix] = val
      else:
         extl=val
   return (extl / 3.57) - 1.0

def lin_interp(vv, xx, xout):
   """
      linear interpolation/extrapolaton
   """
   # Make a copy so we don't overwrite the input arguments.
   v = vv
   x = xx
   interpfunc = interpolate.interp1d(xx,vv, kind='linear')
   vout=interpfunc(xout)
   
   return vout

