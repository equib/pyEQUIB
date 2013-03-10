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

def equib_hgen(xx, gh, y, npt, iopt, ndim, ndimt3, hmh):
   """
   NAME:
       equib_hgen
   PURPOSE:
       Cubic spline interpolation
       The equation for the second derivatives at internal points
       is of the form G*YPP=B, where G has been evaluated and LU
       decomposed.
       this routine writes B=HMH*Y and then solves YPP=G**(-1)*HMH*Y,
       =HMH*Y.
       Three options are provided for boundary conditions-
       IOPT = 0  YPP=0 at end points
       IOPT = 1  YP=0  at end points
       IOPT = 2  YP at end points from lagarnge interpolant of a set of
       internal points.
  
   EXPLANATION:
  
   CALLING SEQUENCE:
       (xx, gh, y, npt, iopt, ndim, ndimt3, hmh) = equib_hgen(xx, gh, y, npt, iopt, ndim, ndimt3, hmh)
  
   INPUTS:
       XX -     XX parameter
       GH -     GH parameter
       Y -      Y parameter
       NPT -    NPT parameter
       IOPT -   IOPT parameter
       NDIM -   NDIM parameter
       NDIMT3 - NDIMT3 parameter
       HMH -    HMH parameter
   REVISION HISTORY:
       Converted from FORTRAN EQUIB to Python, 15/09/2013
   """
   #NPT= long(0)
   #IOPT= long(0)
   #NDIM= long(0)
   #NDIMT3= long(0)
   
   ndim3 = numpy.int32(0)
   nip = numpy.int32(0)
   i = numpy.int32(0)
   j = numpy.int32(0)
   k = numpy.int32(0)
   npm = numpy.int32(0)
   indx = numpy.int32(0)
   #XX=dblarr(NDIM+1)
   #GH=dblarr(NDIMT3+1)
   #Y=dblarr(NDIM+1)
   #HMH=dblarr(NDIM+1,NDIM+1)
   xy = numpy.zeros(5 + 1)
   d = numpy.zeros(5 + 1)
   c = numpy.zeros((2 + 1,5 + 1))
   a0 = numpy.float64(0)
   an1 = numpy.float64(0)
   h1 = numpy.float64(0)
   h2 = numpy.float64(0)
   # Case of derivative boundary condition, with
   if (iopt == 2):   
      # derivatives from NIP-point Lagrange at
      ndim3 = 5
      # internal points
      nip = 3
      for j in range(1, 3):
         for i in range(1, (nip)+(1)):
            k = (npt - nip) * (j - 1)
            xy[i] = xx[k + i]
         k = 1 + (npt - 1) * (j - 1)
         (xy, d, xx[k], nip, ndim3)=equib_deriv(xy, d, xx[k], nip, ndim3)
         for i in range(1, (nip)+(1)):
            c[j,i] = d[i]
   # Set up matrix equation G*YPP=HMH*Y
   a0 = xx[2] - xx[1]
   an1 = xx[npt] - xx[npt - 1]
   npm = npt - 2
   for i in range(1, (npm)+(1)):
      h1 = 6. / (xx[i + 1] - xx[i])
      h2 = 6. / (xx[i + 2] - xx[i + 1])
      for j in range(1, (npt)+(1)):
         hmh[i,j] = 0.
         if (j == i):   
            hmh[i,j] = h1
         if (j == i + 2):   
            hmh[i,j] = h2
         if (j == i + 1):   
            hmh[i,j] = -h1 - h2
   #Correct matrix for case of
   if ((iopt == 1) or (iopt == 2)):   
      # derivative boundary conditions
      hmh[1,1] = hmh[1,1] + 3 / a0
      hmh[1,2] = hmh[1,2] - 3 / a0
      hmh[npm,npt - 1] = hmh[npm,npt - 1] - 3 / an1
      hmh[npm,npt] = hmh[npm,npt] + 3 / an1
   if (iopt == 2):   
      for j in range(1, (nip)+(1)):
         hmh[1,j] = hmh[1,j] + 3 * c[1,j]
         k = npt + j - nip
         hmh[npm,k] = hmh[npm,k] - 3 * c[2,j]
   #for I=1,NPM do begin
   #endfor
   # Solve matrix equation with results in the form
   for i in range(1, (npt)+(1)):
   # YPP=HMH*Y. matrix g has been LU decomposed
      y[1] = hmh[1,i]
      indx = 0
      for j in range(2, (npm)+(1)):
         indx = indx + 3
         y[j] = hmh[j,i] - gh[indx] * y[j - 1]
      indx = indx + 1
      y[npm] = y[npm] / gh[indx]
      for j in range(2, (npm)+(1)):
         k = npm - j + 1
         indx = indx - 3
         y[k] = (y[k] - gh[indx + 1] * y[k + 1]) / gh[indx]
      for j in range(1, (npm)+(1)):
         hmh[j + 1,i] = y[j]
      #Insert values for second derivative at end
      hmh[1,i] = 0.
      # points: first and last rows of the matrix
      hmh[npt,i] = 0.
   # Case of derivative boundary conditions
   if (iopt > 0):   
      for j in range(1, (npt)+(1)):
         hmh[1,j] = -0.5 * hmh[2,j]
         hmh[npt,j] = -0.5 * hmh[npt - 1,j]
      hmh[1,1] = hmh[1,1] - 3 / (a0 * a0)
      hmh[1,2] = hmh[1,2] + 3 / (a0 * a0)
      hmh[npt,npt - 1] = hmh[npt,npt - 1] + 3 / (an1 * an1)
      hmh[npt,npt] = hmh[npt,npt] - 3 / (an1 * an1)
   if (iopt == 2):   
      for j in range(1, (nip)+(1)):
         hmh[1,j] = hmh[1,j] - 3 * c[1,j] / a0
         k = npt + j - nip
         hmh[npt,k] = hmh[npt,k] + 3 * c[2,j] / an1
   #for I=1,NPT do begin
   #endfor
   
   return (xx, gh, y, npt, iopt, ndim, ndimt3, hmh)

def equib_ghgen(gh, xx, npt, iopt, ndim, ndimt3):
   """
   NAME:
       equib_ghgen
   PURPOSE:
  
   EXPLANATION:
  
   CALLING SEQUENCE:
       (gh, xx, npt, iopt, ndim, ndimt3) = equib_ghgen(gh, xx, npt, iopt, ndim, ndimt3)
  
   INPUTS:
       GH -     GH parameter
       XX -     XX parameter
       NPT -    NPT parameter
       IOPT -   IOPT parameter
       NDIM -   NDIM parameter
       NDIMT3 - NDIMT3 parameter
   REVISION HISTORY:
       Converted from FORTRAN EQUIB to Python, 15/09/2013
   """
   #NPT= long(0)
   #IOPT= long(0)
   #NDIM= long(0)
   #NDIMT3= long(0)
   
   indx = numpy.int32(0)
   nptm = numpy.int32(0)
   i = numpy.int32(0)
   j = numpy.int32(0)
   ip = numpy.int32(0)
   jp = numpy.int32(0)
   ik = numpy.int32(0)
   #XX=dblarr(NDIM+1)
   #GH=dblarr(NDIMT3+1)
   indx = 0
   nptm = npt - 1
   for i in range(2, (nptm)+(1)):
      ip = i - 1
      for j in range(1, 4):
         jp = ip + j - 2
         if ((jp >= 1) and (jp <= nptm - 1)):   
            indx = indx + 1
            if (j == 2):   
               gh[indx] = 2 * (xx[i + 1] - xx[i - 1])
            else:   
               ik = i + (j - 1) / 2
               gh[indx] = xx[ik] - xx[ik - 1]
   if (iopt >= 1):   
      gh[1] = gh[1] - (xx[2] - xx[1]) / 2.
      gh[indx] = gh[indx] - (xx[npt] - xx[npt - 1]) / 2.
   
   return (gh, xx, npt, iopt, ndim, ndimt3)

def equib_luslv(a, b, n, m):
   """
   NAME:
       equib_luslv
   PURPOSE:
       Solving linear equations
   EXPLANATION:
  
   CALLING SEQUENCE:
       equib_luslv, A, B, N, M
  
   INPUTS:
       A -     A parameter
       B -     B parameter
       N -     N parameter
       M -     M parameter
   REVISION HISTORY:
       Converted from FORTRAN EQUIB to Python, 15/09/2013
   """
   #M= long(0)
   #N= long(0)
   #A=dblarr(M+1,M+1)
   #B=dblarr(M+1)
   
   (a, n, m)=equib_lured(a, n, m)
   (a, b, n, m)=equib_reslv(a, b, n, m)
   
   return (a, b, n, m)

def equib_lured(a, n, nr):
   """
   NAME:
       equib_lured
   PURPOSE:
  
   EXPLANATION:
  
   CALLING SEQUENCE:
       equib_lured, A, N, NR
  
   INPUTS:
       A -     A parameter
       N -     N parameter
       NR -     NR parameter
   REVISION HISTORY:
       Converted from FORTRAN EQUIB to Python, 15/09/2013
   """
   # N= long(0)
   # NR= long(0)
   
   nm1 = numpy.int32(0)
   i = numpy.int32(0)
   j = numpy.int32(0)
   k = numpy.int32(0)
   ip1 = numpy.int32(0)
   #A=dblarr(NR+1,NR+1)
   fact = numpy.float64(0)
   if (n == 1):   
      return (a, n, nr)
   nm1 = n - 1
   for i in range(1, (nm1)+(1)):
      ip1 = i + 1
      for k in range(ip1, (n)+(1)):
         fact = a[k,i] / a[i,i]
         for j in range(ip1, (n)+(1)):
            a[k,j] = a[k,j] - a[i,j] * fact
   
   return (a, n, nr)

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

