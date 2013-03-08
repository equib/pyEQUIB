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

def equib_cfy(x, y, xx, yy, npt, ndim, hmh, d):
   """
   NAME:
       equib_cfy
   PURPOSE:
  
   EXPLANATION:
  
   CALLING SEQUENCE:
       (x, y, xx, yy, npt, ndim, hmh, d) = equib_cfy(x, y, xx, yy, npt, ndim, hmh, d)
  
   INPUTS:
       X -     XX parameter
       Y -     GH parameter
       XX -    Y parameter
       YY -    NPT parameter
       NPT -   IOPT parameter
       NDIM -  NDIM parameter
       HMH -   NDIMT3 parameter
       D -     HMH parameter
   REVISION HISTORY:
       Converted from FORTRAN EQUIB to Python, 15/09/2013
   """
   #NPT= long(0)
   #NDIM= long(0)
   
   j = numpy.int32(0)
   #;XX=dblarr(NDIM+1)
   #YY=dblarr(NDIM+1)
   #HMH=dblarr(NDIM+1,NDIM+1)
   #D=dblarr(NDIM+1)
   #X= double(0)
   #Y= double(0)
   tt = numpy.float64(0)
   if (x < xx[1]):   
      y = yy[1]
   if (x > xx[npt]):   
      y = yy[npt]
   tt = 0.0
   for j in range(1, (npt)+(1)):
      tt = tt + d[j] * yy[j]
   y = tt
   
   return (x, y, xx, yy, npt, ndim, hmh, d)

def equib_deriv(xy, d, x, n, ndim):
   """
   NAME:
       equib_deriv
   PURPOSE:
       Calculate the first derivative of the lagrangian interpolator
       of a function F, tabulated at the N points XY(I), I=1 to N.
       The derivative is given as the coefficients of F(I), I=1 to N,
       in the array D(I), I=1 to N.
   EXPLANATION:
  
   CALLING SEQUENCE:
       equib_deriv, XY, D, X, N, NDIM
  
   INPUTS:
       XY -     XX parameter
       D -      D parameter
       X -      X parameter
       N -      N parameter
       NDIM -   NDIM parameter
   REVISION HISTORY:
       Converted from FORTRAN EQUIB to Python, 15/09/2013
   """
   #N= long(0)
   #NDIM= long(0)
   
   i = numpy.int32(0)
   j = numpy.int32(0)
   k = numpy.int32(0)
   #XY=dblarr(NDIM+1)
   #D=dblarr(NDIM+1)
   #X=double(0)
   p1 = numpy.float64(0)
   p2 = numpy.float64(0)
   s = numpy.float64(0)
   
   for i in range(1, (n)+(1)):
      p1 = 1.
      s = 0.
      for j in range(1, (n)+(1)):
         if (j != i):   
            p1 = p1 * (xy[i] - xy[j])
            p2 = 1.
            for k in range(1, (n)+(1)):
               if ((k != i) and (k != j)):   
                  p2 = p2 * (x - xy[k])
            s = s + p2
      d[i] = s / p1
   
   return (xy, d, x, n, ndim)

def equib_elu(gh, n, ndim):
   """
   NAME:
       equib_elu
   PURPOSE:
  
   EXPLANATION:
  
   CALLING SEQUENCE:
       (gh, n, ndim) = equib_elu(gh, n, ndim)
  
   INPUTS:
       GH -     GH parameter
       N -      N parameter
       NDIM -   NDIM parameter
   REVISION HISTORY:
       Converted from FORTRAN EQUIB to Python, 15/09/2013
   """
   #N= long(0)
   #NDIM= long(0)
   
   indx = numpy.int32(0)
   i = numpy.int32(0)
   j = numpy.int32(0)
   jp = numpy.int32(0)
   #GH=dblarr(NDIM+1)
   indx = 0
   for i in range(1, (n)+(1)):
      for j in range(1, 4):
         jp = i + j - 2
         if ((jp >= 1) and (jp <= n)):   
            indx = indx + 1
            if (i > 1):   
               if (j == 1):   
                  gh[indx] = gh[indx] / gh[indx - 2]
               if (j == 2):   
                  gh[indx] = gh[indx] - gh[indx - 1] * gh[indx - 2]
   
   return (gh, n, ndim)

