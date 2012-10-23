"""
This module contains functions for plasma diagnostics and abundance analysis 
from optical recombination lines (ORLs)
"""


import numpy, os
import array, math
from scipy import interpolate

def h_balmer_line_ratios(temp, dens, line):
   """
    NAME:
        h_balmer_line_ratios
    PURPOSE:
        return Hydrogen Case B Balmer line ratio for
        given electron temperature and density
        Table 4.4 in D. E. Osterbrock & G. J. Ferland,
        Astrophysics of Gaseius Nebulae and
        Active Galactic Nuclei, 2nd Ed., 2006
    EXPLANATION:
   
    CALLING SEQUENCE:
        H_balmer_theory = h_balmer_line_ratios(temp, dens, line)
   
    INPUTS:
        temp -     electron temperature in K
        dens -     electron density in cm-3
        line -     Balmer line number
           1 : Ha/Hb
           2 : Hg/Hb
           3 : Hd/Hb
           4 : H10/Hb
           5 : H15/Hb
           6 : H20/Hb
    RETURN:  Hydrogen Case B Balmer line ratio
   
    REVISION HISTORY:
        from Table 4.4 in D. E. Osterbrock & 
        G. J. Ferland, Astrophysics of Gaseius Nebulae
        and Active Galactic Nuclei, 2nd Ed., 2006
        Python code by A. Danehkar, 31/08/2012
   """
   dens_grid = numpy.array([1.0e2, 1.0e4, 1.0e6])
   temp_grid = numpy.array([5000.0, 10000.0, 20000.0])
   _expr = line
   # D. E. Osterbrock & G. J. Ferland,
   # Astrophysics of Gaseius Nebulae and Active Galactic Nuclei, 2nd Ed., 2006
   # Table 4.4, Hydrogen Case B Balmer line ratio
   if _expr == 1:   
      hr_grid = numpy.array([[3.041, 3.001, 2.918], [2.863, 2.847, 2.806], [2.747, 2.739, 2.725]]) # Ha/Hb
   elif _expr == 2:   
      hr_grid = numpy.array([[0.458, 0.460, 0.465], [0.468, 0.469, 0.471], [0.475, 0.476, 0.476]]) # Hg/Hb
   elif _expr == 3:   
      hr_grid = numpy.array([[0.251, 0.253, 0.258], [0.259, 0.260, 0.262], [0.264, 0.264, 0.266]]) # Hd/Hb
   elif _expr == 4:   
      hr_grid = numpy.array([[0.0515, 0.0520, 0.0616], [0.0530, 0.0533, 0.0591], [0.0540, 0.0541, 0.0575]]) # H10/Hb
   elif _expr == 5:   
      hr_grid = numpy.array([[0.01534, 0.01628, 0.02602], [0.01561, 0.01620, 0.02147], [0.01576, 0.01612, 0.01834]]) # H15/Hb
   elif _expr == 6:   
      hr_grid = numpy.array([[0.00657, 0.00819, 0.01394], [0.00662, 0.00755, 0.01058], [0.00664, 0.00717, 0.00832]]) # H20/Hb
   else:
      raise RuntimeError('no match found for expression')
   
   # Linearly interpolate extinction law in 1/lam
   if (temp < temp_grid[0] or temp > temp_grid[2]):   
      print 'ouside temperature range!'
      return 0
   if (dens < dens_grid[0] or dens > dens_grid[2]):   
      print 'ouside density range!'
      return 0
   
   # Linearly interpolate density
   hr_dns0 = lin_interp(hr_grid[0,:], dens_grid, dens)
   hr_dns1 = lin_interp(hr_grid[1,:], dens_grid, dens)
   hr_dns2 = lin_interp(hr_grid[2,:], dens_grid, dens)
   
   # Linearly interpolate temperature
   hr_dns_grid = numpy.array([hr_dns0, hr_dns1, hr_dns2])
   hr_tmp = lin_interp(hr_dns_grid, temp_grid, temp)
   
   return hr_tmp

def hb_emissivity(temp, dens):
   """
    NAME:
        hb_emissivity
    PURPOSE:
        return Hydrogen Beta Balmer emissivity for
        given electron temperature and density
        Table 4.4 in D. E. Osterbrock & G. J. Ferland,
        Astrophysics of Gaseius Nebulae and
        Active Galactic Nuclei, 2nd Ed., 2006
    EXPLANATION:
   
    CALLING SEQUENCE:
        h_i = hb_emissivity(temp, dens)
   
    INPUTS:
        temp -     electron temperature in K
        dens -     electron density in cm-3
    RETURN:  Hydrogen Beta Balmer emissivity (Case B)
   
    REVISION HISTORY:
        from Table 4.4 in D. E. Osterbrock & 
        G. J. Ferland, Astrophysics of Gaseius Nebulae
        and Active Galactic Nuclei, 2nd Ed., 2006
        Python code by A. Danehkar, 31/08/2012
   """
   dens_grid = numpy.array([1.0e2, 1.0e4, 1.0e6])
   temp_grid = numpy.array([5000.0, 10000.0, 20000.0])
   hr_grid = numpy.array([numpy.array([2.20, 2.22, 2.29]), numpy.array([1.23, 1.24, 1.25]), numpy.array([0.658, 0.659, 0.661])]) # 4*pi*j_Hb/(ne*np)
   
   hr_grid = 1.0e-25 * hr_grid # erg cm3 s-1
   # Linearly interpolate extinction law in 1/lam
   if (temp < temp_grid[0] or temp > temp_grid[2]):   
      print 'ouside temperature range!'
      return 0
   if (dens < dens_grid[0] or dens > dens_grid[2]):   
      print 'ouside density range!'
      return 0
   
   # Linearly interpolate density
   hr_dns0 = lin_interp(hr_grid[0,:], dens_grid, dens)
   hr_dns1 = lin_interp(hr_grid[1,:], dens_grid, dens)
   hr_dns2 = lin_interp(hr_grid[2,:], dens_grid, dens)
   
   # Linearly interpolate temperature
   hr_dns_grid = numpy.array([hr_dns0, hr_dns1, hr_dns2])
   hr_tmp = lin_interp(hr_dns_grid, temp_grid, temp)
   
   return hr_tmp

def hb_emissivity(temp, dens):
   """
    NAME:
        hb_emissivity
    PURPOSE:
        return Hydrogen Beta Balmer emissivity for
        given electron temperature and density
        Table 4.4 in D. E. Osterbrock & G. J. Ferland,
        Astrophysics of Gaseius Nebulae and
        Active Galactic Nuclei, 2nd Ed., 2006
    EXPLANATION:
   
    CALLING SEQUENCE:
        h_i = hb_emissivity(temp, dens)
   
    INPUTS:
        temp -     electron temperature in K
        dens -     electron density in cm-3
    RETURN:  Hydrogen Beta Balmer emissivity (Case B)
   
    REVISION HISTORY:
        from Table 4.4 in D. E. Osterbrock & 
        G. J. Ferland, Astrophysics of Gaseius Nebulae
        and Active Galactic Nuclei, 2nd Ed., 2006
        Python code by A. Danehkar, 31/08/2012
   """
   dens_grid = numpy.array([1.0e2, 1.0e4, 1.0e6])
   temp_grid = numpy.array([5000.0, 10000.0, 20000.0])
   hr_grid = numpy.array([numpy.array([2.20, 2.22, 2.29]), numpy.array([1.23, 1.24, 1.25]), numpy.array([0.658, 0.659, 0.661])]) # 4*pi*j_Hb/(ne*np)
   
   hr_grid = 1.0e-25 * hr_grid # erg cm3 s-1
   # Linearly interpolate extinction law in 1/lam
   if (temp < temp_grid[0] or temp > temp_grid[2]):   
      print 'ouside temperature range!'
      return 0
   if (dens < dens_grid[0] or dens > dens_grid[2]):   
      print 'ouside density range!'
      return 0
   
   # Linearly interpolate density
   hr_dns0 = lin_interp(hr_grid[0,:], dens_grid, dens)
   hr_dns1 = lin_interp(hr_grid[1,:], dens_grid, dens)
   hr_dns2 = lin_interp(hr_grid[2,:], dens_grid, dens)
   
   # Linearly interpolate temperature
   hr_dns_grid = numpy.array([hr_dns0, hr_dns1, hr_dns2])
   hr_tmp = lin_interp(hr_dns_grid, temp_grid, temp)
   
   return hr_tmp

def hb_eff_rec_coef(temp, dens):
   """
    NAME:
        hb_eff_rec_coef
    PURPOSE:
        return effective recombination coefficient
        of Hydrogen Beta Balmer for
        given electron temperature and density
        Table 4.4 in D. E. Osterbrock & G. J. Ferland,
        Astrophysics of Gaseius Nebulae and
        Active Galactic Nuclei, 2nd Ed., 2006
    EXPLANATION:
   
    CALLING SEQUENCE:
        hb_eff = hb_rec_coef(temp, dens)
   
    INPUTS:
        temp -     electron temperature in K
        dens -     electron density in cm-3
    RETURN:  Hydrogen Beta Balmer emissivity (Case B)
   
    REVISION HISTORY:
        from Table 4.4 in D. E. Osterbrock & 
        G. J. Ferland, Astrophysics of Gaseius Nebulae
        and Active Galactic Nuclei, 2nd Ed., 2006
        Python code by A. Danehkar, 31/08/2012
   """
   dens_grid = numpy.array([1.0e2, 1.0e4, 1.0e6])
   temp_grid = numpy.array([5000.0, 10000.0, 20000.0])
   hr_grid = numpy.array([numpy.array([5.37, 5.43, 5.59]), numpy.array([3.02, 3.03, 3.07]), numpy.array([1.61, 1.61, 1.62])]) # alpha_hb_eff
   
   hr_grid = 1.0e-14 * hr_grid # cm3 s-1
   # Linearly interpolate extinction law in 1/lam
   if (temp < temp_grid[0] or temp > temp_grid[2]):   
      print 'ouside temperature range!'
      return 0
   if (dens < dens_grid[0] or dens > dens_grid[2]):   
      print 'ouside density range!'
      return 0
   
   # Linearly interpolate density
   hr_dns0 = lin_interp(hr_grid[0,:], dens_grid, dens)
   hr_dns1 = lin_interp(hr_grid[1,:], dens_grid, dens)
   hr_dns2 = lin_interp(hr_grid[2,:], dens_grid, dens)
   
   # Linearly interpolate temperature
   hr_dns_grid = numpy.array([hr_dns0, hr_dns1, hr_dns2])
   hr_tmp = lin_interp(hr_dns_grid, temp_grid, temp)
   
   return hr_tmp

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



