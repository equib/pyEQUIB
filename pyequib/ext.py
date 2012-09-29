"""
This module contains functions for extinction corrections using different laws
"""

import numpy, os
import array, math
from scipy import interpolate

"""class _Extinction(object):
    
    def __init__(self):
        pass
""" 

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


