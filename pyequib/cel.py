"""
This module contains functions for plasma diagnostics and abundance analysis 
from collisionally excited lines (CELs)
"""

import numpy, os
import array, math

def equib_cfd(x, xx, npt, ndim, hmh, d):
   """
   NAME:
       equib_cfd
   PURPOSE:
  
   EXPLANATION:
  
   CALLING SEQUENCE:
       (x, xx, npt, ndim, hmh, d)=equib_cfd(x, xx, npt, ndim, hmh, d)
  
   INPUTS:
       X -     X parameter
       XX -    XX parameter
       NPT -   NPT parameter
       NDIM -  NDIM parameter
       HMH -   HMH parameter
       D -     D parameter
   REVISION HISTORY:
       Converted from FORTRAN EQUIB to Python, 15/09/2013
   """
   #NPT= long(0)
   #NDIM= long(0)
   nptm = numpy.int32(0)
   i = numpy.int32(0)
   j = numpy.int32(0)
   #X= double(0)
   #XX=dblarr(NDIM+1)
   #HMH=dblarr(NDIM+1,NDIM+1)
   #D=dblarr(NDIM+1)
   x1 = numpy.float64(0)
   x2 = numpy.float64(0)
   a1 = numpy.float64(0)
   a2 = numpy.float64(0)
   hi = numpy.float64(0)
   if (x < xx[1]):   
      #print, XX[1]
      return (x, xx, npt, ndim, hmh, d)
   if (x > xx[npt]):   
      #print, XX[NPT]
      return (x, xx, npt, ndim, hmh, d)
   nptm = npt - 1
   for i in range(1, (nptm)+(1)):
      if (x < xx[i + 1]):   
         x1 = xx[i + 1] - x
         x2 = x - xx[i]
         hi = xx[i + 1] - xx[i]
         a1 = x1 * (x1 * x1 / (6 * hi) - hi / 6)
         a2 = x2 * (x2 * x2 / (6 * hi) - hi / 6)
         for j in range(1, (npt)+(1)):
            d[j] = a1 * hmh[i,j] + a2 * hmh[i + 1,j]
         d[i] = d[i] + x1 / hi
         d[i + 1] = d[i + 1] + x2 / hi
         return (x, xx, npt, ndim, hmh, d)
   
   return (x, xx, npt, ndim, hmh, d)

