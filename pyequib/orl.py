"""
This module contains functions for plasma diagnostics and abundance analysis 
from optical recombination lines (ORLs)
"""

# A. Danehkar
#
# Version 0.1, 15/08/2016
# First Release
#

import numpy, os
import array, math
from scipy import interpolate

def gamma4861(temp, dens):
   """
    NAME:
        gamm4861
    PURPOSE:
        determine the value of Log10 (gamm(H Beta))
        = Log10( 4*Pai*j(HBeta)/NpNe) at temperature Te and density Ne
        Storey P. J., Hummer D. G., 1995, MNRAS, 272, 41
    EXPLANATION:
   
    CALLING SEQUENCE:
        gamm4861_theory = gamma4861(temp,dens)
   
    INPUTS:
        temp -     electron temperature in K
        dens -     electron density in cm-3
    RETURN:  Log10 (gamm(H Beta))
   
    REVISION HISTORY:
        from Table 4.4 in D. E. Osterbrock & 
        G. J. Ferland, Astrophysics of Gaseius Nebulae
        and Active Galactic Nuclei, 2nd Ed., 2006
        Python code by A. Danehkar, 31/08/2012
   """
   hb_ems = math.log10(hb_emissivity(temp, dens))
   return hb_ems
   
def gamma4471(temp, dens):
   """
    NAME:
        gamma4471
    PURPOSE:
        determine the value of Log10 (gamm(HeI4471))
        = Log10( 4*Pai*j(HeI 4471)/N(He+)Ne) at temperature Te and density Ne
        Smits D. P., 1996, MNRAS, 278, 6
    EXPLANATION:
   
    CALLING SEQUENCE:
        gamm4471_theory = gamma4471(temp,dens)
   
    INPUTS:
        temp -     electron temperature in K
        dens -     electron density in cm-3
    RETURN:  Log10 (gamm(HeI4471))
   
    REVISION HISTORY:
        Python code by A. Danehkar, 31/08/2012
   """
   ems4471 = he_i_emissivity_smits(temp, dens, 4)
   return math.log10(ems4471)
   
def gamma5876(temp, dens):
   """
    NAME:
        gamma5876
    PURPOSE:
        determine the value of Log10 (gamm(HeI5876))
        = Log10( 4*Pai*j(HeI 5876)/N(He+)Ne) at temperature Te and density Ne
        Smits D. P., 1996, MNRAS, 278, 68
    EXPLANATION:
   
    CALLING SEQUENCE:
        gamm5876_theory = gamm5876(temp,dens)
   
    INPUTS:
        temp -     electron temperature in K
        dens -     electron density in cm-3
    RETURN:  Log10 (gamm(HeI5876))
   
    REVISION HISTORY:
        Python code by A. Danehkar, 31/08/2012
   """
   ems5876 = he_i_emissivity_smits(temp, dens, 6)
   return math.log10(ems5876)
   
def gamma6678(temp, dens):
   """
    NAME:
        gamma6678
    PURPOSE:
        determine the value of Log10 (gamm(HeI6678))
        = Log10( 4*Pai*j(HeI 6678)/N(He+)Ne) at temperature Te and density Ne
        Smits D. P., 1996, MNRAS, 278, 683
    EXPLANATION:
   
    CALLING SEQUENCE:
        gamm6678_theory = gamm6678(temp,dens)
   
    INPUTS:
        temp -     electron temperature in K
        dens -     electron density in cm-3
    RETURN:  Log10 (gamm(HeI6678))
   
    REVISION HISTORY:
        Python code by A. Danehkar, 31/08/2012
   """
   ems6678 = he_i_emissivity_smits(temp, dens, 7)
   return math.log10(ems6678)
   
def hb_ems_aeff(temp, dens):
   """
    NAME:
        hb_ems_aeff
    PURPOSE:
        determine the value of Aeff and emissivity of H_beta
        Table 4.4 in D. E. Osterbrock & G. J. Ferland, 
        Astrophysics of Gaseius Nebulae and
        Active Galactic Nuclei, 2nd Ed., 2006
    EXPLANATION:
   
    CALLING SEQUENCE:
        [ems, aeff]=hb_ems_aeff(temp, dens)
   
    INPUTS:
        temp  - electron temperature in K
        dens  - electron density in cm-3
    OUTPUTS:
        {aeff:double(0.0), ems:double(0.0)}
        hbeta.aeff - effective recombination coefficient of H_beta
        hbeta.ems   - emissivity of H_beta
   
    REVISION HISTORY:
        Python code by A. Danehkar, 10/05/2013
   """
   aeff = hb_eff_rec_coef(temp, dens)
   logems = math.log10(hbeta.aeff/double(4861.33/1.98648e-8))
   ems = 10.0 ** logems
   return (ems, aeff)


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
      print('outside temperature range!')
      return 0
   if (dens < dens_grid[0] or dens > dens_grid[2]):   
      print('outside density range!')
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
      print('outside temperature range!')
      return 0
   if (dens < dens_grid[0] or dens > dens_grid[2]):   
      print('outside density range!')
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
      print('outside temperature range!')
      return 0
   if (dens < dens_grid[0] or dens > dens_grid[2]):   
      print('outside density range!')
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
      print('outside temperature range!')
      return 0
   if (dens < dens_grid[0] or dens > dens_grid[2]):   
      print('outside density range!')
      return 0
   
   # Linearly interpolate density
   hr_tmp0 = lin_interp(emiss_grid[0,:], temp_grid, temp)
   hr_tmp1 = lin_interp(emiss_grid[1,:], temp_grid, temp)
   hr_tmp2 = lin_interp(emiss_grid[2,:], temp_grid, temp)
   
   # Linearly interpolate temperature
   hr_tmp_grid = numpy.array([hr_tmp0, hr_tmp1, hr_tmp2])
   hr_tmp = lin_interp(hr_tmp_grid, dens_grid, dens)
   
   return hr_tmp

def h_i_tot_rec_coef_sh(temp, dens, case_name=None):
   """
    NAME:
        h_i_tot_rec_coef_sh
    PURPOSE:
        return total hydrogen recombination coefficient
        given electron temperature and density
        Storey & Hummer, MNRAS, 272, 41S, 1995
    EXPLANATION:
   
    CALLING SEQUENCE:
        hb_eff = h_i_tot_rec_coef_sh(temp, dens)
   
    INPUTS:
        temp -     electron temperature in K
        dens -     electron density in cm-3
    RETURN:  Hydrogen Beta Balmer emissivity (Case B)
   
    REVISION HISTORY:
        Storey & Hummer, MNRAS, 272, 41S, 1995
        1995MNRAS.272...41S
        Python code by A. Danehkar, 31/08/2012
   """
   case_num = 0
   if (case_name is not None):   
      _expr = case_name
      if _expr == 'A':   
         case_num = 1
      elif _expr == 'B':   
         case_num = 0
      else:
         raise RuntimeError('no match found for expression')
   else:   
      case_num = 0
   
   if case_num == 1:   
      # Case A
      dens_grid = numpy.array([1.e2, 1.e3, 1.e4, 1.e5, 1.e6, 1.e7, 1.e8, 1.e9, 1.e10])
      temp_grid = numpy.array([5.e2, 1.e3, 3.e3, 5.e3, 7.5e3, 1.e4, 1.25e4, 1.5e4, 2.e4, 3.e4])
      # alpha_hb_eff ; cm3 s-1
      hr_grid = numpy.array([numpy.array([3.251e-12, 3.351e-12, 3.536e-12, 3.880e-12, 4.552e-12, 5.943e-12, 9.129e-12, 1.769e-11, 4.694e-11]), numpy.array([2.038e-12, 2.069e-12, 2.125e-12, 2.229e-12, 2.424e-12, 2.802e-12, 3.575e-12, 5.308e-12, 9.817e-12]), numpy.array([9.690e-13, 9.735e-13, 9.819e-13, 9.973e-13, 1.026e-12, 1.079e-12, 1.179e-12, 1.374e-12, 1.778e-12]), numpy.array([6.809e-13, 6.827e-13, 6.861e-13, 6.923e-13, 7.038e-13, 7.252e-13, 7.651e-13, 8.406e-13, 9.890e-13]), numpy.array([5.120e-13, 5.128e-13, 5.145e-13, 5.174e-13, 5.230e-13, 5.333e-13, 5.524e-13, 5.883e-13, 6.569e-13]), numpy.array([4.169e-13, 4.174e-13, 4.183e-13, 4.201e-13, 4.234e-13, 4.294e-13, 4.408e-13, 4.619e-13, 5.018e-13]), numpy.array([3.547e-13, 3.550e-13, 3.557e-13, 3.568e-13, 3.590e-13, 3.630e-13, 3.705e-13, 3.844e-13, 4.107e-13]), numpy.array([3.104e-13, 3.106e-13, 3.111e-13, 3.119e-13, 3.134e-13, 3.163e-13, 3.216e-13, 3.315e-13, 3.501e-13]), numpy.array([2.507e-13, 2.509e-13, 2.511e-13, 2.516e-13, 2.525e-13, 2.541e-13, 2.572e-13, 2.630e-13, 2.737e-13]), numpy.array([1.843e-13, 1.844e-13, 1.845e-13, 1.847e-13, 1.851e-13, 1.858e-13, 1.872e-13, 1.898e-13, 1.947e-13])])
      
      # Linearly interpolate extinction law in 1/lam
      if (temp < temp_grid[0] or temp > temp_grid[9]):   
         print('outside temperature range!')
         return 0
      if (dens < dens_grid[0] or dens > dens_grid[8]):   
         print('outside density range!')
         return 0
      
      # Linearly interpolate density
      hr_dns0 = lin_interp(hr_grid[0,:], dens_grid, dens)
      hr_dns1 = lin_interp(hr_grid[1,:], dens_grid, dens)
      hr_dns2 = lin_interp(hr_grid[2,:], dens_grid, dens)
      hr_dns3 = lin_interp(hr_grid[3,:], dens_grid, dens)
      hr_dns4 = lin_interp(hr_grid[4,:], dens_grid, dens)
      hr_dns5 = lin_interp(hr_grid[5,:], dens_grid, dens)
      hr_dns6 = lin_interp(hr_grid[6,:], dens_grid, dens)
      hr_dns7 = lin_interp(hr_grid[7,:], dens_grid, dens)
      hr_dns8 = lin_interp(hr_grid[8,:], dens_grid, dens)
      hr_dns9 = lin_interp(hr_grid[9,:], dens_grid, dens)
      
      
      # Linearly interpolate temperature
      hr_dns_grid = numpy.array([hr_dns0, hr_dns1, hr_dns2, hr_dns3, hr_dns4, hr_dns5, hr_dns6, hr_dns7, hr_dns8, hr_dns9])
      hr_tmp = lin_interp(hr_dns_grid, temp_grid, temp)
   else:   
      # Case B
      dens_grid = numpy.array([1.e2, 1.e3, 1.e4, 1.e5, 1.e6, 1.e7, 1.e8, 1.e9, 1.e10, 1.e11, 1.e12, 1.e13, 1.e14])
      temp_grid = numpy.array([5.e2, 1.e3, 3.e3, 5.e3, 7.5e3, 1.e4, 1.25e4, 1.5e4, 2.e4, 3.e4])
      # alpha_hb_eff ; cm3 s-1
      hr_grid = numpy.array([numpy.array([2.493e-12, 2.573e-12, 2.720e-12, 2.998e-12, 3.542e-12, 4.681e-12, 7.330e-12, 1.462e-11, 4.044e-11, 1.700e-10, 1.119e-09, 9.959e-09, 9.756e-08]), numpy.array([1.512e-12, 1.535e-12, 1.579e-12, 1.658e-12, 1.810e-12, 2.106e-12, 2.717e-12, 4.111e-12, 7.823e-12, 2.037e-11, 7.981e-11, 4.937e-10, 4.266e-09]), numpy.array([6.708e-13, 6.740e-13, 6.798e-13, 6.907e-13, 7.109e-13, 7.486e-13, 8.204e-13, 9.615e-13, 1.257e-12, 1.941e-12, 3.811e-12, 1.040e-11, 4.325e-11]), numpy.array([4.522e-13, 4.534e-13, 4.556e-13, 4.597e-13, 4.674e-13, 4.816e-13, 5.083e-13, 5.592e-13, 6.601e-13, 8.726e-13, 1.371e-12, 2.755e-12, 7.792e-12]), numpy.array([3.273e-13, 3.278e-13, 3.288e-13, 3.306e-13, 3.341e-13, 3.404e-13, 3.524e-13, 3.749e-13, 4.181e-13, 5.047e-13, 6.908e-13, 1.138e-12, 2.453e-12]), numpy.array([2.585e-13, 2.588e-13, 2.594e-13, 2.604e-13, 2.623e-13, 2.658e-13, 2.724e-13, 2.848e-13, 3.083e-13, 3.541e-13, 4.477e-13, 6.551e-13, 1.200e-12]), numpy.array([2.144e-13, 2.147e-13, 2.149e-13, 2.156e-13, 2.167e-13, 2.190e-13, 2.230e-13, 2.307e-13, 2.452e-13, 2.728e-13, 3.276e-13, 4.426e-13, 7.303e-13]), numpy.array([1.836e-13, 1.837e-13, 1.839e-13, 1.843e-13, 1.851e-13, 1.866e-13, 1.893e-13, 1.944e-13, 2.040e-13, 2.221e-13, 2.571e-13, 3.278e-13, 5.036e-13]), numpy.array([1.428e-13, 1.429e-13, 1.430e-13, 1.432e-13, 1.436e-13, 1.444e-13, 1.458e-13, 1.484e-13, 1.532e-13, 1.621e-13, 1.788e-13, 2.106e-13, 2.959e-13]), numpy.array([9.911e-14, 9.913e-14, 9.917e-14, 9.924e-14, 9.937e-14, 9.962e-14, 1.001e-13, 1.009e-13, 1.025e-13, 1.054e-13, 1.103e-13, 1.190e-13, 1.528e-13])])
      
      # Linearly interpolate extinction law in 1/lam
      if (temp < temp_grid[0] or temp > temp_grid[9]):   
         print('outside temperature range!')
         return 0
      if (dens < dens_grid[0] or dens > dens_grid[12]):   
         print('outside density range!')
         return 0
      
      # Linearly interpolate density
      hr_dns0 = lin_interp(hr_grid[0,:], dens_grid, dens)
      hr_dns1 = lin_interp(hr_grid[1,:], dens_grid, dens)
      hr_dns2 = lin_interp(hr_grid[2,:], dens_grid, dens)
      hr_dns3 = lin_interp(hr_grid[3,:], dens_grid, dens)
      hr_dns4 = lin_interp(hr_grid[4,:], dens_grid, dens)
      hr_dns5 = lin_interp(hr_grid[5,:], dens_grid, dens)
      hr_dns6 = lin_interp(hr_grid[6,:], dens_grid, dens)
      hr_dns7 = lin_interp(hr_grid[7,:], dens_grid, dens)
      hr_dns8 = lin_interp(hr_grid[8,:], dens_grid, dens)
      hr_dns9 = lin_interp(hr_grid[9,:], dens_grid, dens)
      
      # Linearly interpolate temperature
      hr_dns_grid = numpy.array([hr_dns0, hr_dns1, hr_dns2, hr_dns3, hr_dns4, hr_dns5, hr_dns6, hr_dns7, hr_dns8, hr_dns9])
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

