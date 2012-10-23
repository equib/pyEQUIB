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

def he_i_emissivity_smits(temp, dens, line):
   """
    NAME:
        he_i_emissivity_smits
    PURPOSE:
        return helium emissivity for
        given electron temperature and density
        Smits D. P., 1996, MNRAS, 278, 683
        1996MNRAS.278..683S
    EXPLANATION:
   
    CALLING SEQUENCE:
        helium_emissivity = he_i_emissivity_smits(temp, dens, line)
   
    INPUTS:
        temp -     electron temperature in K
        dens -     electron density in cm-3
        line -     Balmer line number
           1 : He I 3889
           2 : He I 4026
           3 : He I 4387
           4 : He I 4471
           5 : He I 4922
           6 : He I 5876
           7 : He I 6678
           8 : He I 7065
           9 : He I 7281
           10 : He I 10830
   
    RETURN:  Hydrogen Case B Balmer line ratio
   
    REVISION HISTORY:
        Tables 1 and 2 in Smits 1996, MNRAS, 278, 683
        1996MNRAS.278..683S
        Python code by A. Danehkar, 31/08/2012
   """
   dens_grid = numpy.array([1.0e2, 1.0e4, 1.0e6])
   temp_grid = numpy.array([312.5, 625.0, 1250.0, 2500.0, 5000.0, 10000.0, 20000.0])
   
   # Table 1  in Smits 1996MNRAS.278..683S
   emiss4471_grid = numpy.array([numpy.array([9.364e-25, 5.854e-25, 3.578e-25, 2.104e-25, 1.174e-25, 6.146e-26, 3.001e-26]), numpy.array([1.049e-24, 6.195e-25, 3.674e-25, 2.129e-25, 1.179e-25, 6.155e-26, 3.001e-26]), numpy.array([1.512e-24, 7.501e-25, 4.048e-25, 2.228e-25, 1.202e-25, 6.196e-26, 3.010e-26])])
   
   # Table 2  in Smits 1996MNRAS.278..683S
   # Relative fluxes for He I 3889
   emiss3889_rel_grid = numpy.array([numpy.array([1.262, 1.340, 1.452, 1.617, 1.860, 2.215, 2.722]), numpy.array([1.261, 1.342, 1.455, 1.621, 1.865, 2.219, 2.727]), numpy.array([1.307, 1.382, 1.486, 1.646, 1.886, 2.238, 2.739])])
   
   # Relative fluxes for He I 3889
   emiss4026_rel_grid = numpy.array([numpy.array([0.396, 0.406, 0.419, 0.434, 0.450, 0.466, 0.479]), numpy.array([0.401, 0.408, 0.420, 0.435, 0.451, 0.466, 0.479]), numpy.array([0.413, 0.418, 0.427, 0.439, 0.453, 0.467, 0.480])])
   
   # Relative fluxes for He I 4387
   emiss4387_rel_grid = numpy.array([numpy.array([0.107, 0.110, 0.113, 0.116, 0.120, 0.123, 0.125]), numpy.array([0.109, 0.110, 0.113, 0.117, 0.120, 0.123, 0.125]), numpy.array([0.112, 0.113, 0.115, 0.118, 0.121, 0.124, 0.125])])
   
   # Relative fluxes for He I 4922
   emiss4922_rel_grid = numpy.array([numpy.array([0.273, 0.273, 0.272, 0.271, 0.269, 0.266, 0.263]), numpy.array([0.274, 0.273, 0.272, 0.270, 0.269, 0.266, 0.263]), numpy.array([0.274, 0.273, 0.272, 0.271, 0.269, 0.266, 0.263])])
   
   # Relative fluxes for He I 5876
   emiss5876_rel_grid = numpy.array([numpy.array([4.382, 4.035, 3.655, 3.295, 2.984, 2.752, 2.542]), numpy.array([4.047, 3.846, 3.548, 3.235, 2.952, 2.715, 2.534]), numpy.array([3.656, 3.532, 3.364, 3.131, 2.894, 2.686, 2.519])])
   
   # Relative fluxes for He I 6678
   emiss6678_rel_grid = numpy.array([numpy.array([1.274, 1.171, 1.058, 0.950, 0.856, 0.777, 0.714]), numpy.array([1.181, 1.116, 1.026, 0.932, 0.847, 0.772, 0.712]), numpy.array([1.068, 1.024, 0.974, 0.902, 0.830, 0.764, 0.708])])
   
   # Relative fluxes for He I 7065
   emiss7065_rel_grid = numpy.array([numpy.array([0.215, 0.230, 0.254, 0.293, 0.356, 0.461, 0.639]), numpy.array([0.216, 0.230, 0.254, 0.294, 0.356, 0.461, 0.639]), numpy.array([0.220, 0.234, 0.257, 0.295, 0.357, 0.461, 0.637])])
   
   # Relative fluxes for He I 7281
   emiss7281_rel_grid = numpy.array([numpy.array([0.028, 0.031, 0.035, 0.042, 0.054, 0.075, 0.112]), numpy.array([0.028, 0.031, 0.035, 0.042, 0.054, 0.075, 0.112]), numpy.array([0.027, 0.030, 0.035, 0.042, 0.054, 0.075, 0.112])])
   
   # Relative fluxes for He I 10830
   emiss10830_rel_grid = numpy.array([numpy.array([4.223, 4.072, 3.980, 3.990, 4.313, 5.584, 8.557]), numpy.array([4.033, 3.969, 3.928, 4.702, 13.10, 38.08, 90.87]), numpy.array([3.884, 3.853, 4.034, 7.482, 20.41, 51.69, 117.1])])
   
   _expr = line
   if _expr == 1:   
      emiss_grid = emiss3889_rel_grid * emiss4471_grid # He I 3889
   elif _expr == 2:   
      emiss_grid = emiss4026_rel_grid * emiss4471_grid # He I 4026
   elif _expr == 3:   
      emiss_grid = emiss4387_rel_grid * emiss4471_grid # He I 4387
   elif _expr == 4:   
      emiss_grid = emiss4471_grid # He I 4471
   elif _expr == 5:   
      emiss_grid = emiss4922_rel_grid * emiss4471_grid # He I 4922
   elif _expr == 6:   
      emiss_grid = emiss5876_rel_grid * emiss4471_grid # He I 5876
   elif _expr == 7:   
      emiss_grid = emiss6678_rel_grid * emiss4471_grid # He I 6678
   elif _expr == 8:   
      emiss_grid = emiss7065_rel_grid * emiss4471_grid # He I 7065
   elif _expr == 9:   
      emiss_grid = emiss7281_rel_grid * emiss4471_grid # He I 7281
   elif _expr == 10:   
      emiss_grid = emiss10830_rel_grid * emiss4471_grid # He I 10830
   else:
      raise RuntimeError('no match found for expression')
   
   # Linearly interpolate extinction law in 1/lam
   if (temp < temp_grid[0] or temp > temp_grid[6]):   
      print 'ouside temperature range!'
      return 0
   if (dens < dens_grid[0] or dens > dens_grid[2]):   
      print 'ouside density range!'
      return 0
   
   # Linearly interpolate density
   hr_tmp0 = lin_interp(emiss_grid[0,:], temp_grid, temp)
   hr_tmp1 = lin_interp(emiss_grid[1,:], temp_grid, temp)
   hr_tmp2 = lin_interp(emiss_grid[2,:], temp_grid, temp)
   
   # Linearly interpolate temperature
   hr_tmp_grid = numpy.array([hr_tmp0, hr_tmp1, hr_tmp2])
   hr_tmp = lin_interp(hr_tmp_grid, dens_grid, dens)
   
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



