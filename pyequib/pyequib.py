# -*- coding: utf-8 -*-

"""
This module contains functions for Plasma Diagnostics and Abundance Analysis
"""

# A. Danehkar
#
# Version 0.2.0, 01/10/2020
#
#

import numpy as np
#import pandas as pd
import atomneb
from scipy import interpolate

__all__ =  ["calc_temperature", "calc_density",
            "calc_populations", "calc_crit_density",
            "calc_emissivity", "calc_abundance",
            "print_ionic", "get_omij_temp",
            "calc_emiss_h_beta",
            "calc_emiss_he_i_rl", "calc_emiss_he_ii_rl",
            "calc_emiss_c_ii_rl", "calc_emiss_c_iii_rl",
            "calc_emiss_n_ii_rl", "calc_emiss_n_iii_rl",
            "calc_emiss_o_ii_rl", "calc_emiss_ne_ii_rl",
            "calc_abund_he_i_rl", "calc_abund_he_ii_rl",
            "calc_abund_c_ii_rl", "calc_abund_c_iii_rl",
            "calc_abund_n_ii_rl", "calc_abund_n_iii_rl",
            "calc_abund_o_ii_rl", "calc_abund_ne_ii_rl",
            "redlaw", "redlaw_gal", "redlaw_gal2",
            "redlaw_ccm", "redlaw_jbk", "redlaw_fm",
            "redlaw_smc", "redlaw_lmc",
            "deredden_flux", "deredden_relflux"]

def calc_temperature(line_flux_ratio=None, density=None,
                     upper_levels=None, lower_levels=None,
                     elj_data=None, omij_data=None, aij_data=None,
                     low_temperature=None, high_temperature=None,
                     num_temperature=None, min_density=None):
   """
        This function determines electron temperature from given
        flux intensity ratio for specified ion with upper level(s)
        lower level(s) by solving atomic level populations and
        line emissivities in statistical equilibrium
        for given electron density.

    :Returns:
       type=double. This function returns the electron temperature.

    :Keywords:
        line_flux_ratio  :     in, required, type=float
                               flux intensity ratio
        density          :     in, required, type=float
                               electron density
        upper_levels     :     in, required, type=string,
                               upper atomic level(s) e.g '1,2/', '1,2,1,3/'
        lower_levels     :     in, required, type=string
                               lower atomic level(s) e.g '1,2/', '1,2,1,3/'
        elj_data         :     in, required, type=array/object
                               energy levels (Ej) data
        omij_data        :     in, required, type=array/object
                               collision strengths (omega_ij) data
        aij_data         :     in, required, type=array/object
                               transition probabilities (Aij) data
        low_temperature  :     in, optional, type=float
                               lower temperature range
        high_temperature  :     in, optional, type=float
                               upper temperature range
        num_temperature  :     in, optional, type=integer
                               number of the iteration step
        min_density      :     in, optional, type=float
                               lower density range

    :Examples:
       For example::

        >>> import pyequib
        >>> import atomneb
        >>> import os
        >>> base_dir = '../externals/atomneb/'
        >>> data_dir = os.path.join('atomic-data', 'chianti70')
        >>> atom_elj_file = os.path.join(base_dir,data_dir, 'AtomElj.fits')
        >>> atom_omij_file = os.path.join(base_dir,data_dir, 'AtomOmij.fits')
        >>> atom_aij_file = os.path.join(base_dir,data_dir, 'AtomAij.fits')
        >>> atom='s'
        >>> ion='ii'
        >>> s_ii_elj=atomneb.read_elj(atom_elj_file, atom, ion, level_num=5) # read Energy Levels (Ej)
        >>> s_ii_omij=atomneb.read_omij(atom_omij_file, atom, ion) # read Collision Strengths (Omegaij)
        >>> s_ii_aij=atomneb.read_aij(atom_aij_file, atom, ion) # read Transition Probabilities (Aij)
        >>> upper_levels='1,2,1,3/'
        >>> lower_levels='1,5/'
        >>> density = np.float64(2550)
        >>> line_flux_ratio=np.float64(10.753)
        >>> temperature=pyequib.calc_temperature(line_flux_ratio=line_flux_ratio, density=density,
        >>>                              upper_levels=upper_levels, lower_levels=lower_levels,
        >>>                              elj_data=s_ii_elj, omij_data=s_ii_omij,
        >>>                              aij_data=s_ii_aij)
        >>> print("Electron Temperature:", temperature)
           Electron Temperature:       7920.2865

    :Categories:
      Plasma Diagnostics, Collisionally Excited Lines

    :Dirs:
     ./
         Main routines

    :Author:
      Ashkbiz Danehkar

    :Copyright:
      This library is released under a GNU General Public License.

    :Version:
      0.3.0

    :History:
        15/09/2013, A. Danehkar, Translated from FORTRAN to IDL code.

        20/10/2016, A. Danehkar, Replaced str2int with strnumber.

        20/10/2016, A. Danehkar, Replaced CFY, SPLMAT, and CFD with
             IDL function INTERPOL( /SPLINE).

        20/10/2016, A. Danehkar, Replaced LUSLV with IDL LAPACK function
                          LA_LINEAR_EQUATION.

        15/11/2016, A. Danehkar, Replaced LA_LINEAR_EQUATION (not work in GDL)
              with IDL function LUDC & LUSOL.

        19/11/2016, A. Danehkar, Replaced INTERPOL (not accurate) with
                       SPL_INIT & SPL_INTERP.

        20/11/2016, A. Danehkar, Made a new function calc_populations()
          for solving atomic level populations and separated it from
          calc_abundance(), calc_density() and calc_temperature().

        10/03/2017, A. Danehkar, Integration with AtomNeb, now uses atomic data
                         input elj_data, omij_data, aij_data.

        12/06/2017, A. Danehkar, Cleaning the function, and remove unused varibales
                           from calc_temperature().

        27/02/2019, A. Danehkar, Fix a bug in the atomic level assumption, and
                           use the simplified calc_populations() routine.

        04/03/2019, A. Danehkar, Use the get_omij_temp() routine.

        24/05/2019, A. Danehkar, Add the optional temperature range.

        03/10/2020, A. Danehkar, Transferred from IDL to Python.

    FORTRAN HISTORY:

        03/05/1981, I.D.Howarth,  Version 1.

        05/05/1981, I.D.Howarth,  Minibug fixed!

        07/05/1981, I.D.Howarth,  Now takes collision rates or strengths.

        03/08/1981, S.Adams,      Interpolates collision strengths.

        07/08/1981, S.Adams,      Input method changed.

        19/11/1984, R.E.S.Clegg,  SA files entombed in scratch disk. Logical
                                  filenames given to SA's data files.

        08/1995, D.P.Ruffle, Changed input file format. Increased matrices.

        02/1996, X.W.Liu,   Tidy up. SUBROUTINES SPLMAT, HGEN, CFY and CFD
                            modified such that matrix sizes (i.e. maximum
                            of Te and maximum no of levels) can now be cha
                            by modifying the parameters NDIM1, NDIM2 and N
                            in the Main program. EASY!
                            Now takes collision rates as well.
                            All variables are declared explicitly
                            Generate two extra files (ionpop.lis and ionra
                            of plain stream format for plotting.

        06/1996, C.J.Pritchet, Changed input data format for cases IBIG=1,2.
                            Fixed readin bug for IBIG=2 case.
                            Now reads reformatted upsilons (easier to see
                            and the 0 0 0 data end is excluded for these c
                            The A values have a different format for IBIG=.

        2006, B.Ercolano,   Converted to F90.
   """
   #common share1, Atomic_Data_Path
   h_planck = 6.62606957e-27 # erg s
   c_speed = 2.99792458e10 # cm/s
   
   if (line_flux_ratio is not None) == 0:   
      print('flux intensity ratio is not given')
      return 0
   if (density is not None) == 0:   
      print('Density is not set')
      return 0
   if (elj_data is not None) == 0:   
      print('Energy Levels data (elj_data) are not set')
      return 0
   if (omij_data is not None) == 0:   
      print('Collision Strengths (omij_data) are not set')
      return 0
   if (aij_data is not None) == 0:   
      print('Transition Probabilities (aij_data) are not set')
      return 0
   if (upper_levels is not None) == 0:   
      print('Upper levels (upper_levels) are not given')
      return 0
   if (lower_levels is not None) == 0:   
      print('Lower levels (lower_levels) are not given')
      return 0
   if (density <= 0.e0):   
      print('density = ', density)
      return 0
   if (low_temperature is not None):   
      temp_min = low_temperature
   else:   
      temp_min = 5000.0
   if (high_temperature is not None):   
      temp_max = high_temperature
   else:   
      temp_max = 20000.0
   if (num_temperature is not None):   
      temp_num = num_temperature
   else:   
      temp_num = 4
   if (min_density is not None):   
      dens_min = min_density
   else:   
      dens_min = 1.0
   
   iteration = np.int32(0)
   
   level_num = np.int32(0)
   int1 = np.int32(0)
   ind = np.int32(0)
   it = np.int32(0)
   
   tempi = np.float64(0)
   tinc = np.float64(0)
   densi = np.float64(0)
   dinc = np.float64(0)
   temperature = np.float64(0)
   eji = np.float64(0)
   wav = np.float64(0)
   emis_sum_a = np.float64(0)
   emis_sum_b = np.float64(0)
   qx = np.float64(0)
   ax = np.float64(0)
   ex = np.float64(0)
   frat = np.float64(0)
   dee = np.float64(0)
   ltext = ''#
   
   result1 = np.float64(0)
   
   level_num = len(elj_data)
   t_num = len(omij_data['strength'][0])
   omij_num = len(omij_data)
   
   wava = np.zeros(level_num + 1)
   wavb = np.zeros(level_num + 1)
   omij = np.zeros((level_num, level_num, t_num))
   check_value = np.zeros(2)
   
   label1 = (level_num + 1)*['']
   
   upper_levels_str = do_strsplit(upper_levels, ',',escapech='/')
   lower_levels_str = do_strsplit(lower_levels, ',',escapech='/')
   
   upper_levels_num = np.int32(len(upper_levels_str)/2)
   lower_levels_num = np.int32(len(lower_levels_str)/2)
   
   itrana = np.zeros((2 + 1, upper_levels_num + 1))
   itranb = np.zeros((2 + 1, lower_levels_num + 1))

   itrana[:,:] = 0
   itranb[:,:] = 0

   upper_levels_i = np.int32(0)
   for i in range(0, upper_levels_num):
      itrana[0,i] = do_str2int(upper_levels_str[upper_levels_i])
      itrana[1,i] = do_str2int(upper_levels_str[upper_levels_i + 1])
      upper_levels_i = upper_levels_i + 2
      #if upper_levels_i >= upper_levels_num:
      #   break

   lower_levels_i = np.int32(0)
   for i in range(0, lower_levels_num):
      itranb[0,i] = do_str2int(lower_levels_str[lower_levels_i])
      itranb[1,i] = do_str2int(lower_levels_str[lower_levels_i + 1])
      lower_levels_i = lower_levels_i + 2
      #if lower_levels_i >= lower_levels_num:
      #   break

   irats = 0
   #level_max=max([max(ITRANA),max(ITRANB)]) ! mistake
   level_max = level_num
   aij = aij_data['aij'][0]
   aij=aij.T
   elj = elj_data['ej']
   # set temperature iterations
   # start of iterations
   # ****************************
   for iteration in range(1, 10):
      if (iteration == 1):   
         tempi = temp_min
      else:   
         tempi = check_value[0]
      int1 = temp_num
      tinc = (temp_max - temp_min) / ((int1 - 1) ** (iteration))
      #int1=8
      #TINC=(55000.0)/((int1-1)^(iteration))
      densi = density
      if (densi <= dens_min):   
         densi = dens_min
      results = np.zeros((2, int1))
      if (tempi < temp_min):   
         tempi = temp_min # add
      # Start of temperature iteration
      for jt in range(1, (int1)+(1)):
         temperature = tempi + (jt - 1) * tinc
         if ((temperature <= 0.e0) | (density <= 0.e0)):
            print('temperature = ', temperature, ', density = ', density)
            return 0
         if level_max > level_num:   
            print("error outside level range")
            return 0
         omij_t = get_omij_temp(temperature=temperature, omij_data=omij_data, level_num=level_max, irats=irats)
         nlj = calc_populations(temperature=temperature, density=densi, elj_data=elj_data, omij_data=omij_data, aij_data=aij_data, eff_omij=omij_t, level_num=level_max, irats=irats)
         # Search ITRANA, ITRANB for transitions & sum up
         emis_sum_a = np.float64(0)
         emis_sum_b = np.float64(0)
         for ikt in range(0, (upper_levels_num - 1)+(1)):
            i = np.int32(itrana[0,ikt])
            j = np.int32(itrana[1,ikt])
            emissivity_line = np.float64(0)
            if (aij[j - 1,i - 1] != 0.e0):
               eji = elj[j - 1] - elj[i - 1]
               wav = 1.e8 / eji
               emissivity_line = nlj[j - 1] * aij[j - 1,i - 1] * h_planck * c_speed * 1.e8 / wav
               emis_sum_a = emis_sum_a + emissivity_line
         for ikt in range(0, (lower_levels_num - 1)+(1)):
            i = np.int32(itranb[0,ikt])
            j = np.int32(itranb[1,ikt])
            emissivity_line = np.float64(0)
            if (aij[j - 1,i - 1] != 0.e0):
               eji = elj[j - 1] - elj[i - 1]
               wav = 1.e8 / eji
               emissivity_line = nlj[j - 1] * aij[j - 1,i - 1] * h_planck * c_speed * 1.e8 / wav
               emis_sum_b = emis_sum_b + emissivity_line
         frat = emis_sum_a / emis_sum_b
         results[0, jt - 1] = temperature
         results[1, jt - 1] = frat - line_flux_ratio
         
         for ia in range(0, upper_levels_num):
            i1 = np.int32(itrana[0,ia])
            i2 = np.int32(itrana[1,ia])
            dee = elj[i2 - 1] - elj[i1 - 1]
            wava[ia] = 1.e8 / dee
         for ib in range(0, lower_levels_num):
            i1 = np.int32(itranb[0,ib])
            i2 = np.int32(itranb[1,ib])
            dee = elj[i2 - 1] - elj[i1 - 1]
            wavb[ib] = 1.e8 / dee
         # End of the temperature iteration
      # iteration and detect the sign change.
      for i in range(2, (int1)+(1)):
         check = 0
         if (check_sign(results[1, i - 1], results[1, 0]) != results[1, i - 1]):
            #if this condition, the values have a different sign
            check_value[:] = results[:, i - 2] # the value before the sign change returned
            check = 1
            break
      if ((check == 0) & (iteration < 9)):    # check if there is any change of sign,
         #and checks if it should be upper or lower limit
         if (abs(results[1, 0])) < (abs(results[1, int1 - 1])):
            check_value[:] = results[0,:]
         else:   
            if (abs(results[1, int1 - 1]) < abs(results[1, 0])):
               check_value[:] = results[:, int1 - 2]
            else:   
               print('check_value is wrong')
               return 0
      else:   
         if ((check == 0) and (iteration == 9)):    #check if no change of sign,
            #and checks if it should be upper or lower limit
            if (abs(results[1, 0]) < abs(results[1, int1 - 1])):
               check_value[:] = results[:,0]
            else:   
               if (abs(results[1, int1 - 1]) < abs(results[1, 0])):
                  check_value[:] = results[:, int1 - 1]
               else:   
                  print('check_value is wrong')
                  return 0
   # end of iterations
   #****************************
   result1 = check_value[0]
   return result1

def calc_density(line_flux_ratio=None, temperature=None,
                 upper_levels=None, lower_levels=None,
                 elj_data=None, omij_data=None, aij_data=None,
                 low_density=None, high_density=None,
                 num_density=None, min_temperature=None):
   """
        This function determines electron density from given
        flux intensity ratio for specified ion with upper level(s)
        lower level(s) by solving atomic level populations and
        line emissivities in statistical equilibrium
        for given electron temperature.

    :Returns:
       type=double. This function returns the electron density.

    :Keywords:
        line_flux_ratio  :     in, required, type=float
                               flux intensity ratio
        temperature      :     in, required, type=float
                               electron temperature
        upper_levels     :     in, required, type=string
                               upper atomic level(s) e.g '1,2/', '1,2,1,3/'
        lower_levels     :     in, required, type=string
                               lower atomic level(s) e.g '1,2/', '1,2,1,3/'
        elj_data         :     in, required, type=array/object
                               energy levels (Ej) data
        omij_data        :     in, required, type=array/object
                               collision strengths (omega_ij) data
        aij_data         :     in, required, type=array/object
                               transition probabilities (Aij) data
        low_density      :     in, optional, type=float
                               lower density range
        high_density      :     in, optional, type=float
                               upper density range
        num_density      :     in, optional, type=integer
                               number of the iteration step
        min_temperature  :     in, optional, type=float
                               minimum temperature

    :Examples:
       For example::

        >>> import pyequib
        >>> import atomneb
        >>> import os
        >>> base_dir = '../externals/atomneb/'
        >>> data_dir = os.path.join('atomic-data', 'chianti70')
        >>> atom_elj_file = os.path.join(base_dir,data_dir, 'AtomElj.fits')
        >>> atom_omij_file = os.path.join(base_dir,data_dir, 'AtomOmij.fits')
        >>> atom_aij_file = os.path.join(base_dir,data_dir, 'AtomAij.fits')
        >>> atom='s'
        >>> ion='ii'
        >>> s_ii_elj=atomneb.read_elj(atom_elj_file, atom, ion, level_num=5) # read Energy Levels (Ej)
        >>> s_ii_omij=atomneb.read_omij(atom_omij_file, atom, ion) # read Collision Strengths (Omegaij)
        >>> s_ii_aij=atomneb.read_aij(atom_aij_file, atom, ion) # read Transition Probabilities (Aij)\
        >>> upper_levels='1,2/'
        >>> lower_levels='1,3/'
        >>> temperature=np.float64(7000.0)#
        >>> line_flux_ratio=np.float64(1.506)#
        >>> density=pyequib.calc_density(line_flux_ratio=line_flux_ratio, temperature=temperature,
        >>>                      upper_levels=upper_levels, lower_levels=lower_levels,
        >>>                      elj_data=s_ii_elj, omij_data=s_ii_omij,
        >>>                      aij_data=s_ii_aij)
        >>> print("Electron Density:", density)
           Electron Density:       2312.6395

    :Categories:
      Plasma Diagnostics, Collisionally Excited Lines

    :Dirs:
     ./
         Main routines

    :Author:
      Ashkbiz Danehkar

    :Copyright:
      This library is released under a GNU General Public License.

    :Version:
      0.3.0

    :History:
        15/09/2013, A. Danehkar, Translated from FORTRAN to IDL code.

        20/10/2016, A. Danehkar, Replaced str2int with strnumber.

        20/10/2016, A. Danehkar, Replaced CFY, SPLMAT, and CFD with
             IDL function INTERPOL( /SPLINE).

        20/10/2016, A. Danehkar, Replaced LUSLV with IDL LAPACK function
                          LA_LINEAR_EQUATION.

        15/11/2016, A. Danehkar, Replaced LA_LINEAR_EQUATION (not work in GDL)
              with IDL function LUDC & LUSOL.

        19/11/2016, A. Danehkar, Replaced INTERPOL (not accurate) with
                       SPL_INIT & SPL_INTERP.

        20/11/2016, A. Danehkar, Made a new function calc_populations()
          for solving atomic level populations and separated it from
          calc_abundance(), calc_density() and calc_temperature().

        10/03/2017, A. Danehkar, Integration with AtomNeb, now uses atomic data
                         input elj_data, omij_data, aij_data.

        12/06/2017, A. Danehkar, Cleaning the function, and remove unused varibales
                           from calc_density().

        27/02/2019, A. Danehkar, Fix a bug in the atomic level assumption, and
                           use the simplified calc_populations() routine.

        04/03/2019, A. Danehkar, Use the get_omij_temp() routine.

        24/05/2019, A. Danehkar, Add the optional density range.

        03/10/2020, A. Danehkar, Transferred from IDL to Python.

    FORTRAN HISTORY:

        03/05/1981, I.D.Howarth,  Version 1.

        05/05/1981, I.D.Howarth,  Minibug fixed!

        07/05/1981, I.D.Howarth,  Now takes collision rates or strengths.

        03/08/1981, S.Adams,      Interpolates collision strengths.

        07/08/1981, S.Adams,      Input method changed.

        19/11/1984, R.E.S.Clegg,  SA files entombed in scratch disk. Logical
                                  filenames given to SA's data files.

        08/1995, D.P.Ruffle, Changed input file format. Increased matrices.

        02/1996, X.W.Liu,   Tidy up. SUBROUTINES SPLMAT, HGEN, CFY and CFD
                            modified such that matrix sizes (i.e. maximum
                            of Te and maximum no of levels) can now be cha
                            by modifying the parameters NDIM1, NDIM2 and N
                            in the Main program. EASY!
                            Now takes collision rates as well.
                            All variables are declared explicitly
                            Generate two extra files (ionpop.lis and ionra
                            of plain stream format for plotting.

        06/1996, C.J.Pritchet, Changed input data format for cases IBIG=1,2.
                            Fixed readin bug for IBIG=2 case.
                            Now reads reformatted upsilons (easier to see
                            and the 0 0 0 data end is excluded for these c
                            The A values have a different format for IBIG=.

        2006, B.Ercolano,   Converted to F90.
   """
   # common share1, Atomic_Data_Path
   h_planck = 6.62606957e-27 # erg s
   c_speed = 2.99792458e10 # cm/s

   if (line_flux_ratio is not None) == 0:
      print('flux intensity ratio is not given')
      return 0
   if (temperature is not None) == 0:
      print('Temperature is not set')
      return 0
   if (elj_data is not None) == 0:
      print('Energy Levels data (elj_data) are not set')
      return 0
   if (omij_data is not None) == 0:
      print('Collision Strengths (omij_data) are not set')
      return 0
   if (aij_data is not None) == 0:
      print('Transition Probabilities (aij_data) are not set')
      return 0
   if (upper_levels is not None) == 0:
      print('Upper levels (upper_levels) are not given')
      return 0
   if (lower_levels is not None) == 0:
      print('Lower levels (lower_levels) are not given')
      return 0
   if (temperature <= 0.e0):
      print('temperature = ', temperature)
      return 0
   if (low_density is not None):
      dens_min = low_density
   else:
      dens_min = 1.0
   if (high_density is not None):
      dens_max = high_density
   else:
      dens_max = 100000.0
   if (num_density is not None):
      dens_num = num_density
   else:
      dens_num = 4
   if (min_temperature is not None):
      temp_min = min_temperature
   else:
      temp_min = 5000.0

   iteration = np.int32(0)

   level_num = np.int32(0)
   int1 = np.int32(0)
   ind = np.int32(0)
   it = np.int32(0)

   tempi = np.float64(0)
   densi = np.float64(0)
   dinc = np.float64(0)
   density = np.float64(0)
   eji = np.float64(0)
   wav = np.float64(0)
   emis_sum_a = np.float64(0)
   emis_sum_b = np.float64(0)
   qx = np.float64(0)
   ax = np.float64(0)
   ex = np.float64(0)
   frat = np.float64(0)
   dee = np.float64(0)
   ltext = ''#

   result1 = np.float64(0)

   level_num = len(elj_data)
   t_num = len(omij_data['strength'][0])
   omij_num = len(omij_data)

   wava = np.zeros(level_num + 1)
   wavb = np.zeros(level_num + 1)
   omij = np.zeros((level_num, level_num, t_num))
   check_value = np.zeros(2)

   label1 = (level_num + 1)*['']

   upper_levels_str = do_strsplit(upper_levels, ',',escapech='/')
   lower_levels_str = do_strsplit(lower_levels, ',',escapech='/')

   upper_levels_num = np.int32(len(upper_levels_str)/2)
   lower_levels_num = np.int32(len(lower_levels_str)/2)

   itrana = np.zeros((2 + 1, upper_levels_num + 1))
   itranb = np.zeros((2 + 1, lower_levels_num + 1))

   itrana[:,:] = 0
   itranb[:,:] = 0

   upper_levels_i = np.int32(0)
   for i in range(0, upper_levels_num ):
      itrana[0,i] = do_str2int(upper_levels_str[upper_levels_i])
      itrana[1,i] = do_str2int(upper_levels_str[upper_levels_i + 1])
      upper_levels_i = upper_levels_i + 2
      #if upper_levels_i >= upper_levels_num:
      #   break

   lower_levels_i = np.int32(0)
   for i in range(0, lower_levels_num):
      itranb[0,i] = do_str2int(lower_levels_str[lower_levels_i])
      itranb[1,i] = do_str2int(lower_levels_str[lower_levels_i + 1])
      lower_levels_i = lower_levels_i + 2
      #if lower_levels_i >= lower_levels_num:
      #   break

   irats = 0
   #level_max=max([max(ITRANA),max(ITRANB)]) ! mistake
   level_max = level_num
   aij = aij_data['aij'][0]
   aij = aij.T
   elj = elj_data['ej']
   tempi = temperature
   if (tempi < temp_min):
      tempi = temp_min # add

   omij_t = get_omij_temp(temperature=tempi, omij_data=omij_data, level_num=level_num, irats=irats)
   # set density iterations
   # start of iterations
   # ****************************
   for iteration in range(1, 10):
      if (iteration == 1):
         densi = dens_min
      else:
         densi = check_value[0]
      ind = dens_num
      dinc = (dens_max - dens_min) / ((ind - 1) ** (iteration))
      #IND=8
      #DINC=(1000000.0)/((IND-1)^(iteration))
      results = np.zeros((2, ind))
      if (densi <= dens_min):
         densi = dens_min
      # Start of density iteration
      for jjd in range(1, (ind)+(1)):
         density = densi + (jjd - 1) * dinc
         if ((temperature <= 0.e0) | (density <= 0.e0)):
            print('temperature = ', temperature, ', density = ', density)
            return 0
         if level_max > level_num:
            print("error outside level range")
            return 0
         nlj = calc_populations(temperature=temperature, density=density, elj_data=elj_data, omij_data=omij_data, aij_data=aij_data, eff_omij=omij_t, level_num=level_max, irats=irats)
         # Search ITRANA, ITRANB for transitions & sum up
         emis_sum_a = np.float64(0)
         emis_sum_b = np.float64(0)
         for ikt in range(0, upper_levels_num):
            i = np.int32(itrana[0,ikt])
            j = np.int32(itrana[1,ikt])
            emissivity_line = np.float64(0)
            if (aij[j - 1,i - 1] != 0.e0):
               eji = elj[j - 1] - elj[i - 1]
               wav = 1.e8 / eji
               emissivity_line = nlj[j - 1] * aij[j - 1,i - 1] * h_planck * c_speed * 1.e8 / wav
               emis_sum_a = emis_sum_a + emissivity_line
         for ikt in range(0, lower_levels_num):
            i = np.int32(itranb[0,ikt])
            j = np.int32(itranb[1,ikt])
            emissivity_line = np.float64(0)
            if (aij[j - 1,i - 1] != 0.e0):
               eji = elj[j - 1] - elj[i - 1]
               wav = 1.e8 / eji
               emissivity_line = nlj[j - 1] * aij[j - 1,i - 1] * h_planck * c_speed * 1.e8 / wav
               emis_sum_b = emis_sum_b + emissivity_line
         frat = emis_sum_a / emis_sum_b
         results[0, jjd - 1] = density
         results[1, jjd - 1] = frat - line_flux_ratio
      # End of the denity iteration
      for ia in range(0, (upper_levels_num - 1)+(1)):
         i1 = np.int32(itrana[ia,0])
         i2 = np.int32(itrana[ia,1])
         dee = elj[i2 - 1] - elj[i1 - 1]
         wava[ia] = 1.e8 / dee
      for ib in range(0, (lower_levels_num - 1)+(1)):
         i1 = np.int32(itranb[ib,0])
         i2 = np.int32(itranb[ib,1])
         dee = elj[i2 - 1] - elj[i1 - 1]
         wavb[ib] = 1.e8 / dee
      int1 = ind
      # iteration and detect the sign change.
      for i in range(2, (int1)+(1)):
         check = 0
         if (check_sign(results[1, i - 1], results[1, 0]) != results[1, i - 1]):
            #if this condition, the values have a different sign
            check_value[:] = results[:,i - 2] # the value before the sign change returned
            check = 1
            break
      if ((check == 0) & (iteration < 9)):    # check if there is any change of sign,
         #and checks if it should be upper or lower limit
         if (abs(results[1, 0])) < (abs(results[1, int1 - 1])):
            check_value[:] = results[:,0]
         else:
            if (abs(results[1, int1 - 1]) < abs(results[1, 0])):
               check_value[:] = results[:, int1 - 2]
            else:
               print('check_value is wrong')
               return 0
      else:
         if ((check == 0) & (iteration == 9)):    #check if no change of sign,
            #and checks if it should be upper or lower limit
            if (abs(results[1, 0]) < abs(results[1, int1 - 1])):
               check_value[:] = results[:, 0]
            else:
               if (abs(results[1, int1 - 1]) < abs(results[1, 0])):
                  check_value[:] = results[:, int1 - 1]
               else:
                  print('check_value is wrong')
                  return 0
   # end of iterations
   #****************************
   result1 = check_value[0]
   return result1

def calc_populations(temperature=None, density=None,
                     elj_data=None, omij_data=None, aij_data=None,
                     eff_omij=None, level_num=None, irats=None):
   """
        This function solves atomic level populations in statistical equilibrium
        for given electron temperature and density.
   
    :Returns:
       type=array/object. This function returns the atomic level populations.
   
    :Keywords:
        temperature :   in, required, type=float
                        electron temperature
        density     :   in, required, type=float
                        electron density
        elj_data    :   in, required, type=array/object
                               energy levels (Ej) data
        omij_data   :   in, required, type=array/object
                               collision strengths (omega_ij) data
        aij_data    :   in, required, type=array/object
                               transition probabilities (Aij) data
        eff_Omij    :   in, type=array/object
                        effective collision strengths (Omij_T) at given temperature
        level_num   :   in, type=int
                        Number of levels
        irats       :   in, type=int
                        Else Coll. rates = tabulated values * 10 ** irats
   
    :Examples:
       For example::
   
        >>> import pyequib
        >>> import atomneb
        >>> import os
        >>> base_dir = '../externals/atomneb/'
        >>> data_dir = os.path.join('atomic-data', 'chianti70')
        >>> atom_elj_file = os.path.join(base_dir,data_dir, 'AtomElj.fits')
        >>> atom_omij_file = os.path.join(base_dir,data_dir, 'AtomOmij.fits')
        >>> atom_aij_file = os.path.join(base_dir,data_dir, 'AtomAij.fits')
        >>> atom='s'
        >>> ion='ii'
        >>> s_ii_elj=atomneb.read_elj(atom_elj_file, atom, ion, level_num=5) # read Energy Levels (Ej)
        >>> s_ii_omij=atomneb.read_omij(atom_omij_file, atom, ion) # read Collision Strengths (Omegaij)
        >>> s_ii_aij=atomneb.read_aij(atom_aij_file, atom, ion) # read Transition Probabilities (Aij)\
        >>> density = np.float64(1000)
        >>> temperature=np.float64(10000.0)#
        >>> nlj=pyequib.calc_populations(temperature=temperature, density=density,
        >>>                      elj_data=s_ii_elj, omij_data=s_ii_omij,
        >>>                      aij_data=s_ii_aij)
        >>> print('Atomic Level Populations:', nlj)
           Atomic Level Populations:    0.96992832    0.0070036315     0.023062261   2.6593671e-06   3.1277019e-06
   
    :Categories:
      Plasma Diagnostics, Abundance Analysis, Collisionally Excited Lines
   
    :Dirs:
     ./
         Subroutines
   
    :Author:
      Ashkbiz Danehkar
   
    :Copyright:
      This library is released under a GNU General Public License.
   
    :Version:
      0.3.0
   
    :History:
        15/09/2013, A. Danehkar, Translated from FORTRAN to IDL code.
   
        20/10/2016, A. Danehkar, Replaced str2int with strnumber.
   
        20/10/2016, A. Danehkar, Replaced CFY, SPLMAT, and CFD with
             IDL function INTERPOL( /SPLINE).
   
        20/10/2016, A. Danehkar, Replaced LUSLV with IDL LAPACK function
                          LA_LINEAR_EQUATION.
   
        15/11/2016, A. Danehkar, Replaced LA_LINEAR_EQUATION (not work in GDL)
              with IDL function LUDC & LUSOL.
   
        19/11/2016, A. Danehkar, Replaced INTERPOL (not accurate) with
                       SPL_INIT & SPL_INTERP.
   
        20/11/2016, A. Danehkar, Made a new function calc_populations()
          for solving atomic level populations and separated it from
          calc_abundance(), calc_density() and calc_temperature().
   
        10/03/2017, A. Danehkar, Integration with AtomNeb, now uses atomic data
                         input elj_data, omij_data, aij_data.
   
        12/06/2017, A. Danehkar, Cleaning the function, and remove unused varibales
                           from calc_populations().
   
        27/02/2019, A. Danehkar, Simplify the calc_populations() routine
                           for external usage.
   
        04/03/2019, A. Danehkar, Use the get_omij_temp() routine.

        03/10/2020, A. Danehkar, Transferred from IDL to Python.

    FORTRAN HISTORY:
   
        03/05/1981, I.D.Howarth,  Version 1.
   
        05/05/1981, I.D.Howarth,  Minibug fixed!
   
        07/05/1981, I.D.Howarth,  Now takes collision rates or strengths.
   
        03/08/1981, S.Adams,      Interpolates collision strengths.
   
        07/08/1981, S.Adams,      Input method changed.
   
        19/11/1984, R.E.S.Clegg,  SA files entombed in scratch disk. Logical
                                  filenames given to SA's data files.
   
        08/1995, D.P.Ruffle, Changed input file format. Increased matrices.
   
        02/1996, X.W.Liu,   Tidy up. SUBROUTINES SPLMAT, HGEN, CFY and CFD
                            modified such that matrix sizes (i.e. maximum
                            of Te and maximum no of levels) can now be cha
                            by modifying the parameters NDIM1, NDIM2 and N
                            in the Main program. EASY!
                            Now takes collision rates as well.
                            All variables are declared explicitly
                            Generate two extra files (ionpop.lis and ionra
                            of plain stream format for plotting.
   
        06/1996, C.J.Pritchet, Changed input data format for cases IBIG=1,2.
                            Fixed readin bug for IBIG=2 case.
                            Now reads reformatted upsilons (easier to see
                            and the 0 0 0 data end is excluded for these c
                            The A values have a different format for IBIG=.
   
        2006, B.Ercolano,   Converted to F90.
   """
   h_planck = 4.13566766225e-15 # eV.s #6.62606957e-27 # erg.s
   c_speed = 2.99792458e10 # cm/s
   k_b = 8.617330350e-5 # eV/K # 1.3806485279e-16 # erg/K
   
   pi = 3.1415926535897931e0
   h_bar_planck = 1.054571800e-27 #erg.s
   me = 9.10938356e-28 # gram
   k_b_erg = 1.38064852e-16 # erg/K
   beta1 = np.float64(h_bar_planck ** 2 * np.sqrt(2 * pi / (k_b_erg * me)) / me) # 8.629D-06
   
   if (temperature is not None) == 0:   
      print('Temperature is not set')
      return 0
   if (density is not None) == 0:   
      print('Density is not set')
      return 0
   if (elj_data is not None) == 0:   
      print('elj_data is not set')
      return 0
   if (omij_data is not None) == 0:   
      print('omij_data is not set')
      return 0
   if (aij_data is not None) == 0:   
      print('aij_data is not set')
      return 0
   if (level_num is not None) == 0:   
      level_num = len(elj_data)
   t_num = len(omij_data['strength'][0]) # Number of temperature intervals
   if (irats is not None) == 0:   
      irats = 0
   
   t_lin_list = omij_data['strength'][0]
   t_log_list = np.log10(t_lin_list) # temperature intervals (array)
   
   aij =aij_data['aij'][0] # Transition Probabilities (A_ij)
   aij=aij.T
   ej = elj_data['ej'] # Energy Levels (E_j) in cm-1
   gj = np.int32(elj_data['j_v'] * 2. + 1.) # Ground Levels (G_j)
   
   qj = np.zeros(t_num)
   qij = np.zeros((level_num, level_num))
   equilib_eqs = np.zeros((level_num, level_num))
   nlj = np.zeros(level_num)
   
   qij[:,:] = np.float64(0)
   equilib_eqs[:,:] = np.float64(0)
   nlj[:] = np.float64(0)
   
   t_log = np.log10(temperature)
   
   if (t_num == 1):   
      print('Coll. strengths available for 1 Te only - assuming const')
   else:   
      if (t_num == 2):   
         print('Coll. strengths available for 2 Te only - linear interp')
   # Derive the interpolated effective collision strengths (Omij_T) from collision strengths data (Omij)
   # Obtain collisional de-excitation and excitation rates (Qij) from the effective collision strengths (Omij_T)
   if (eff_omij is not None) == 0:   
      omij_t = get_omij_temp(temperature=temperature, omij_data=omij_data, level_num=level_num, irats=irats)
   else:   
      omij_t = eff_omij
   omij_t=np.float64(omij_t.T)
   for i in range(2, level_num+1):
      for j in range(i, level_num+1):
         d_e = np.float64(ej[j - 1] - ej[i - 2]) * h_planck * c_speed # delta Energy in eV# convert from cm-1 to eV
         # Calculate the Boltzmann factor
         exp_de_kt = np.float64(np.exp(-d_e / (k_b * temperature))) # Maxwell-Boltzmann distribution
         # Obtain collisional de-excitation and excitation rates from the effective collision strengths Omij_T
         if (irats == 0):
            qij[i - 2,j - 1] = beta1 * omij_t[i - 2,j - 1] * exp_de_kt / (np.float64(gj[i - 2]) * np.sqrt(temperature)) # collisional excitation rates
            qij[j - 1,i - 2] = beta1 * omij_t[i - 2,j - 1] / (np.float64(gj[j - 1]) * np.sqrt(temperature)) # collisional de-excitation rates
         else:   
            qij[i - 2,j - 1] = omij_t[i - 2,j - 1] * exp_de_kt * 10. ** irats # collisional excitation rates
            qij[j - 1,i - 2] = np.float64(gj[i - 2]) * qij[i - 2,j - 1] / (exp_de_kt * np.float64(gj[j - 1])) # collisional de-excitation rates
   # Calculate the critical densities
   # N_crit_i = Sum_{j} (Aij) / Sum_{j} (Qij)
   #A_i_sum = TOTAL(Aij, 2) # Sum each of the columns in Aij
   #Q_i_sum = TOTAL(Qij, 2) # Sum each of the columns in Qij
   #N_crit=A_i_tot/Q_i_tot # critical densities
   for i in range(2, level_num+1):
      for j in range(1, level_num+1):
         if (j != i):   
            # the equations for the equilibrium level populations:
            # collisional de-excitation eqs -  collisional excitation eqs = 0
            # Sum_{j ne i} (Ne * Nj * Qji) + Sum_{j > i} (Nj Aji)
            #    - (Sum_{j ne i} (Ne * Ni * Qij) + Sum_{j < i} (Ni Aij)) = 0
            equilib_eqs[i - 1,j - 1] = equilib_eqs[i - 1,j - 1] + density * qij[j - 1,i - 1] # collisional de-excitation
            equilib_eqs[i - 1,i - 1] = equilib_eqs[i - 1,i - 1] - density * qij[i - 1,j - 1] # collisional excitation
            if (j > i):   
               equilib_eqs[i - 1,j - 1] = equilib_eqs[i - 1,j - 1] + aij[j - 1,i - 1] # collisional de-excitation
            else:
               equilib_eqs[i - 1,i - 1] = equilib_eqs[i - 1,i - 1] - aij[i - 1,j - 1] # collisional excitation
   b0 = np.float64(-equilib_eqs[1:(level_num - 1)+1,0])
   equilib_eqs[0:(level_num - 2)+1,0:(level_num - 2)+1] = equilib_eqs[1:(level_num - 1)+1,1:(level_num - 1)+1]
   a0 = (equilib_eqs[0:(level_num - 2)+1,0:(level_num - 2)+1])
   # Solve the equations for the equilibrium level populations
   # A.X = B
   # A: Matrix for the equilibrium level populations equations (i,j), Equilib_Eqs[*,*]
   # B: Vector for the equilibrium level populations equations (i,j=0), Equilib_Eqs[*,0]
   # X: Ionic population density (j), Nj
   # X=la_linear_equation(A0, B0[0:level_num-2])# this function does not work in GDL!
   # ludc(a0, index0)  # Decompose A0# supported by GDL
   # obtain the ionic population densities (j)
   # x0 = lusol(a0, index0, b0[0:(level_num - 2)+1]) # Compute the solution X# supported by GDL
   x0 = np.linalg.solve(a0, b0[0:(level_num - 2)+1])
   # Calculate the atomic level populations (Nlj) from the ionic population densities (Nj)
   nlj[1:(level_num - 1)+1] = x0[0:(level_num - 2)+1]
   nlj[0] = 1.e0
   n_tot = np.float64(sum(nlj[0:(level_num - 1)+1]))
   nlj[0:(level_num - 1)+1] = nlj[0:(level_num - 1)+1] / n_tot
   return nlj

def calc_crit_density(temperature=None,
                      elj_data=None, omij_data=None, aij_data=None,
                      level_num=None, irats=None):
   """
        This function calculates critical densities in statistical equilibrium
        for given electron temperature.
   
    :Returns:
       type=array/object. This function returns the critical densities.
   
    :Keywords:
        temperature :   in, required, type=float
                        electron temperature
        elj_data    :   in, required, type=array/object
                               energy levels (Ej) data
        omij_data   :   in, required, type=array/object
                               collision strengths (omega_ij) data
        aij_data    :   in, required, type=array/object
                               transition probabilities (Aij) data
        level_num   :   in, type=int
                        Number of levels
        irats       :   in, type=int
                        Else Coll. rates = tabulated values * 10 ** irats
   
    :Examples:
       For example::
   
        >>> import pyequib
        >>> import atomneb
        >>> import os
        >>> base_dir = '../externals/atomneb/'
        >>> data_dir = os.path.join('atomic-data', 'chianti70')
        >>> atom_elj_file = os.path.join(base_dir,data_dir, 'AtomElj.fits')
        >>> atom_omij_file = os.path.join(base_dir,data_dir, 'AtomOmij.fits')
        >>> atom_aij_file = os.path.join(base_dir,data_dir, 'AtomAij.fits')
        >>> atom='s'
        >>> ion='ii'
        >>> s_ii_elj=atomneb.read_elj(atom_elj_file, atom, ion, level_num=5) # read Energy Levels (Ej)
        >>> s_ii_omij=atomneb.read_omij(atom_omij_file, atom, ion) # read Collision Strengths (Omegaij)
        >>> s_ii_aij=atomneb.read_aij(atom_aij_file, atom, ion) # read Transition Probabilities (Aij)\
        >>> temperature=np.float64(10000.0)
        >>> n_crit=pyequib.calc_crit_density(temperature=temperature,
        >>>                          elj_data=s_ii_elj, omij_data=s_ii_omij,
        >>>                          aij_data=s_ii_aij)
        >>> print('Critical Densities:', n_crit)
           Critical Densities:       0.0000000       5007.8396       1732.8414       1072685.0       2220758.1
   
    :Categories:
      Plasma Diagnostics, Abundance Analysis, Collisionally Excited Lines
   
    :Dirs:
     ./
         Main routines
   
    :Author:
      Ashkbiz Danehkar
   
    :Copyright:
      This library is released under a GNU General Public License.
   
    :Version:
      0.3.0
   
    :History:
        15/09/2013, A. Danehkar, Translated from FORTRAN to IDL code.
   
        20/10/2016, A. Danehkar, Replaced str2int with strnumber.
   
        20/10/2016, A. Danehkar, Replaced CFY, SPLMAT, and CFD with
             IDL function INTERPOL( /SPLINE).
   
        20/10/2016, A. Danehkar, Replaced LUSLV with IDL LAPACK function
                          LA_LINEAR_EQUATION.
   
        15/11/2016, A. Danehkar, Replaced LA_LINEAR_EQUATION (not work in GDL)
              with IDL function LUDC & LUSOL.
   
        19/11/2016, A. Danehkar, Replaced INTERPOL (not accurate) with
                       SPL_INIT & SPL_INTERP.
   
        20/11/2016, A. Danehkar, Made a new function calc_populations()
          for solving atomic level populations and separated it from
          calc_abundance(), calc_density() and calc_temperature().
   
        10/03/2017, A. Danehkar, Integration with AtomNeb, now uses atomic data
                         input elj_data, omij_data, aij_data.
   
        12/06/2017, A. Danehkar, Cleaning the function, and remove unused varibales
                           from calc_populations().
   
        27/02/2019, A. Danehkar, Simplify the calc_populations() routine
                           for external usage.
   
        01/03/2019, A. Danehkar, Create the calc_crit_density() routine
                           from the calc_populations() routine.
   
        04/03/2019, A. Danehkar, Use the get_omij_temp() routine.

        03/10/2020, A. Danehkar, Transferred from IDL to Python.

    FORTRAN HISTORY:
   
        03/05/1981, I.D.Howarth,  Version 1.
   
        05/05/1981, I.D.Howarth,  Minibug fixed!
   
        07/05/1981, I.D.Howarth,  Now takes collision rates or strengths.
   
        03/08/1981, S.Adams,      Interpolates collision strengths.
   
        07/08/1981, S.Adams,      Input method changed.
   
        19/11/1984, R.E.S.Clegg,  SA files entombed in scratch disk. Logical
                                  filenames given to SA's data files.
   
        08/1995, D.P.Ruffle, Changed input file format. Increased matrices.
   
        02/1996, X.W.Liu,   Tidy up. SUBROUTINES SPLMAT, HGEN, CFY and CFD
                            modified such that matrix sizes (i.e. maximum
                            of Te and maximum no of levels) can now be cha
                            by modifying the parameters NDIM1, NDIM2 and N
                            in the Main program. EASY!
                            Now takes collision rates as well.
                            All variables are declared explicitly
                            Generate two extra files (ionpop.lis and ionra
                            of plain stream format for plotting.
   
        06/1996, C.J.Pritchet, Changed input data format for cases IBIG=1,2.
                            Fixed readin bug for IBIG=2 case.
                            Now reads reformatted upsilons (easier to see
                            and the 0 0 0 data end is excluded for these c
                            The A values have a different format for IBIG=.
   
        2006, B.Ercolano,   Converted to F90.
   """
   h_planck = 4.13566766225e-15 # eV.s #6.62606957e-27 # erg.s
   c_speed = 2.99792458e10 # cm/s
   k_b = 8.617330350e-5 # eV/K # 1.3806485279e-16 # erg/K
   
   pi = 3.1415926535897931e0
   h_bar_planck = 1.054571800e-27 #erg.s
   me = 9.10938356e-28 # gram
   k_b_erg = 1.38064852e-16 # erg/K
   beta1 = h_bar_planck ** 2 * np.sqrt(2 * pi / (k_b_erg * me)) / me # 8.629D-06
   
   if (temperature is not None) == 0:   
      print('Temperature is not set')
      return 0
   if (elj_data is not None) == 0:   
      print('elj_data is not set')
      return 0
   if (omij_data is not None) == 0:   
      print('omij_data is not set')
      return 0
   if (aij_data is not None) == 0:   
      print('aij_data is not set')
      return 0
   if (level_num is not None) == 0:   
      level_num = len(elj_data)
   t_num = len(omij_data['strength'][0]) # Number of temperature intervals
   omij_num = len(omij_data)
   omij = np.zeros((t_num, level_num, level_num))
   for k in range(1, (omij_num - 1)+(1)):
      i = omij_data['level1'][k]
      j = omij_data['level2'][k]
      if (i <= level_num)&( j <= level_num):
         omij[0:(t_num - 1)+1,i - 1,j - 1] = omij_data['strength'][k]
   if (irats is not None) == 0:   
      irats = 0
   
   t_lin_list = omij_data['strength'][0]
   t_log_list = np.log10(t_lin_list) # temperature intervals (array)

   aij = aij_data['aij'][0] # Transition Probabilities (A_ij)
   aij=aij.T
   ej = elj_data['ej'] # Energy Levels (E_j) in cm-1
   gj = np.int32(elj_data['j_v'] * 2. + 1.)  # Ground Levels (G_j)
   
   qj = np.zeros(t_num)
   qij = np.zeros((level_num, level_num))
   equilib_eqs = np.zeros((level_num, level_num))
   nlj = np.zeros(level_num)
   
   qij[:,:] = np.float64(0)
   equilib_eqs[:,:] = np.float64(0)
   nlj[:] = np.float64(0)
   
   t_log = np.log10(temperature)
   
   if (t_num == 1):   
      print('Coll. strengths available for 1 Te only - assuming const')
   else:   
      if (t_num == 2):   
         print('Coll. strengths available for 2 Te only - linear interp')
   # Derive the interpolated effective collision strengths (Omij_T) from collision strengths data (Omij)
   # Obtain collisional de-excitation and excitation rates (Qij) from the effective collision strengths (Omij_T)
   omij_t = get_omij_temp(temperature=temperature, omij_data=omij_data, level_num=level_num, irats=irats)
   omij_t = np.float64(omij_t.T)
   for i in range(2, (level_num)+(1)):
      for j in range(i, (level_num)+(1)):
         d_e = np.float64(ej[j - 1] - ej[i - 2]) * h_planck * c_speed # delta Energy in eV# convert from cm-1 to eV
         # Calculate the Boltzmann factor
         exp_de_kt = np.exp(-d_e / (k_b * temperature)) # Maxwell-Boltzmann distribution
         # Obtain collisional de-excitation and excitation rates from the effective collision strengths Omij_T
         if (irats == 0):   
            qij[i - 2,j - 1] = beta1 * omij_t[i - 2,j - 1] * exp_de_kt / (np.float64(gj[i - 2]) * np.sqrt(temperature)) # collisional excitation rates
            qij[j - 1,i - 2] = beta1 * omij_t[i - 2,j - 1] / (np.float64(gj[j - 1]) * np.sqrt(temperature)) # collisional de-excitation rates
         else:   
            qij[i - 2,j - 1] = omij_t[i - 2,j - 1] * exp_de_kt * 10. ** irats # collisional excitation rates
            qij[j - 1,i - 2] = np.float64(gj[i - 2]) * qij[i - 2,j - 1] / (exp_de_kt * np.float64(gj[j - 1])) # collisional de-excitation rates
   # Calculate the critical densities
   # N_crit_i = Sum_{j} (Aij) / Sum_{j} (Qij)
   a_i_sum = np.sum(aij.T, axis=0) # Sum each of the columns in Aij
   q_i_sum = np.sum(qij.T, axis=0) # Sum each of the columns in Qij
   len1=np.amin([len(a_i_sum), len(q_i_sum)])
   n_crit = a_i_sum[0:len1] / q_i_sum[0:len1] # critical densities
   return n_crit

def calc_emissivity(temperature=None, density=None,
                    atomic_levels=None,
                    elj_data=None, omij_data=None, aij_data=None):
   """
        This function calculates line emissivities for specified ion with level(s) by
        solving atomic level populations and in statistical equilibrium
        for given electron density and temperature.
   
    :Returns:
       type=double. This function returns the line emissivity.
   
    :Keywords:
        temperature   :     in, required, type=float
                            electron temperature
        density       :     in, required, type=float
                            electron density
        atomic_levels :     In, required, type=string
                            level(s) e.g '1,2/', '1,2,1,3/'
        elj_data      :     in, required, type=array/object
                            energy levels (Ej) data
        omij_data     :     in, required, type=array/object
                            collision strengths (omega_ij) data
        aij_data      :     in, required, type=array/object
                            transition probabilities (Aij) data
   
    :Examples:
       For example::
   
        >>> import pyequib
        >>> import atomneb
        >>> import os
        >>> base_dir = '../externals/atomneb/'
        >>> data_dir = os.path.join('atomic-data', 'chianti70')
        >>> atom_elj_file = os.path.join(base_dir,data_dir, 'AtomElj.fits')
        >>> atom_omij_file = os.path.join(base_dir,data_dir, 'AtomOmij.fits')
        >>> atom_aij_file = os.path.join(base_dir,data_dir, 'AtomAij.fits')
        >>> atom='o'
        >>> ion='iii'
        >>> o_iii_elj=atomneb.read_elj(atom_elj_file, atom, ion, level_num=5) # read Energy Levels (Ej)
        >>> o_iii_omij=atomneb.read_omij(atom_omij_file, atom, ion) # read Collision Strengths (Omegaij)
        >>> o_iii_aij=atomneb.read_aij(atom_aij_file, atom, ion) # read Transition Probabilities (Aij)
        >>> temperature=np.float64(10000.0)
        >>> density=np.float64(5000.0)
        >>> atomic_levels='3,4/'
        >>> emiss5007=np.float64(0.0)
        >>> emiss5007=pyequib.calc_emissivity(temperature=temperature, density=density,
        >>>                           atomic_levels=atomic_levels,
        >>>                           elj_data=o_iii_elj, omij_data=o_iii_omij,
        >>>                           aij_data=o_iii_aij
        >>> print('Emissivity(O III 5007):', emiss5007)
           Emissivity(O III 5007):   3.6041012e-21
   
    :Categories:
      Abundance Analysis, Collisionally Excited Lines, Emissivity
   
    :Dirs:
     ./
         Main routines
   
    :Author:
      Ashkbiz Danehkar
   
    :Copyright:
      This library is released under a GNU General Public License.
   
    :Version:
      0.3.0
   
    :History:
        15/09/2013, A. Danehkar, Translated from FORTRAN to IDL code.
   
        20/10/2016, A. Danehkar, Replaced str2int with strnumber.
   
        20/10/2016, A. Danehkar, Replaced CFY, SPLMAT, and CFD with
             IDL function INTERPOL( /SPLINE).
   
        20/10/2016, A. Danehkar, Replaced LUSLV with IDL LAPACK function
                          LA_LINEAR_EQUATION.
   
        15/11/2016, A. Danehkar, Replaced LA_LINEAR_EQUATION (not work in GDL)
              with IDL function LUDC & LUSOL.
   
        19/11/2016, A. Danehkar, Replaced INTERPOL (not accurate) with
                       SPL_INIT & SPL_INTERP.
   
        20/11/2016, A. Danehkar, Made a new function calc_populations()
          for solving atomic level populations and separated it from
          calc_abundance(), calc_density() and calc_temperature().
   
        21/11/2016, A. Danehkar, Made a new function calc_emissivity()
                         for calculating line emissivities and separated it
                         from calc_abundance().
   
        10/03/2017, A. Danehkar, Integration with AtomNeb, now uses atomic data
                         input elj_data, omij_data, aij_data.
   
        12/06/2017, A. Danehkar, Cleaning the function, and remove unused varibales
                           from calc_emissivity().
   
        27/06/2019, A. Danehkar, Use the simplified calc_populations() routine.

        03/10/2020, A. Danehkar, Transferred from IDL to Python.

    FORTRAN HISTORY:
   
        03/05/1981, I.D.Howarth,  Version 1.
   
        05/05/1981, I.D.Howarth,  Minibug fixed!
   
        07/05/1981, I.D.Howarth,  Now takes collision rates or strengths.
   
        03/08/1981, S.Adams,      Interpolates collision strengths.
   
        07/08/1981, S.Adams,      Input method changed.
   
        19/11/1984, R.E.S.Clegg,  SA files entombed in scratch disk. Logical
                                  filenames given to SA's data files.
   
        08/1995, D.P.Ruffle, Changed input file format. Increased matrices.
   
        02/1996, X.W.Liu,   Tidy up. SUBROUTINES SPLMAT, HGEN, CFY and CFD
                            modified such that matrix sizes (i.e. maximum
                            of Te and maximum no of levels) can now be cha
                            by modifying the parameters NDIM1, NDIM2 and N
                            in the Main program. EASY!
                            Now takes collision rates as well.
                            All variables are declared explicitly
                            Generate two extra files (ionpop.lis and ionra
                            of plain stream format for plotting.
   
        06/1996, C.J.Pritchet, Changed input data format for cases IBIG=1,2.
                            Fixed readin bug for IBIG=2 case.
                            Now reads reformatted upsilons (easier to see
                            and the 0 0 0 data end is excluded for these c
                            The A values have a different format for IBIG=.
   
        2006, B.Ercolano,   Converted to F90.
   """
   #global atomic_data_path
   
   h_planck = 6.62606957e-27 # erg s
   c_speed = 2.99792458e10 # cm/s                    
   if (temperature is not None) == 0:   
      print('Temperature is not set')
      return 0
   if (density is not None) == 0:   
      print('Density is not set')
      return 0
   if (elj_data is not None) == 0:   
      print('Energy Levels data (elj_data) are not set')
      return 0
   if (omij_data is not None) == 0:   
      print('Collision Strengths (omij_data) are not set')
      return 0
   if (aij_data is not None) == 0:   
      print('Transition Probabilities (aij_data) are not set')
      return 0
   if (atomic_levels is not None) == 0:   
      print('Atomic levels (atomic_levels) are not given')
      return 0
   if ((temperature <= 0.e0) | (density <= 0.e0)):
      print('temperature = ', temperature, ', density = ', density)
      return 0
   
   irats = np.int32(0)
   itemp = np.int32(0)
   ikt = np.int32(0)
   
   eji = np.float64(0)
   wav = np.float64(0)

   levels_str = do_strsplit(atomic_levels, ',', escapech='/')
   levels_num = np.int32(len(levels_str) / 2)
   itranc = np.zeros((2 + 1, levels_num + 1))
   itranc[:,:] = 0
   levels_i = 0

   levels_i = np.int32(0)
   for i in range(0, levels_num):
       itranc[0, i] = do_str2int(levels_str[levels_i])
       itranc[1, i] = do_str2int(levels_str[levels_i + 1])
       levels_i = levels_i + 2
       # if levels_i >= levels_num:
       #   break
   irats = 0
   aij = aij_data['aij'][0]
   aij=aij.T
   elj = elj_data['ej']
   
   if ((temperature <= 0.e0) | (density <= 0.e0)):
      print('temperature = ', temperature, ', density = ', density)
      return 0
   
   nlj = calc_populations(temperature=temperature, density=density, elj_data=elj_data, omij_data=omij_data, aij_data=aij_data, irats=irats)
   
   emissivity_all = np.float64(0)
   for ikt in range(0, levels_num):
      i = np.int32(itranc[0,ikt])
      j = np.int32(itranc[1,ikt])
      emissivity_line = np.float64(0)
      if (aij[j - 1,i - 1] != 0.e0):
         eji = elj[j - 1] - elj[i - 1]
         wav = 1.e8 / eji
         emissivity_line = nlj[j - 1] * aij[j - 1,i - 1] * h_planck * c_speed * 1.e8 / (wav * density)
         emissivity_all = emissivity_all + emissivity_line
   return emissivity_all

def calc_abundance(temperature=None, density=None,
                   line_flux=None, atomic_levels=None,
                   elj_data=None, omij_data=None, aij_data=None,
                   h_i_aeff_data=None):
   """
        This function determines the ionic abundance from the observed
        flux intensity for specified ion with level(s)
        by solving atomic level populations and
        line emissivities in statistical equilibrium
        for given electron density and temperature.
   
    :Returns:
       type=double. This function returns the ionic abundanc.
   
    :Keywords:
        temperature   :     in, required, type=float
                            electron temperature
        density       :     in, required, type=float
                            electron density
        line_flux     :     in, required, type=float
                            line flux intensity
        atomic_levels :     in, required, type=string
                            level(s) e.g '1,2/', '1,2,1,3/'
        elj_data      :     in, required, type=array/object
                            energy levels (Ej) data
        omij_data     :     in, required, type=array/object
                            collision strengths (omega_ij) data
        aij_data      :     in, required, type=array/object
                            transition probabilities (Aij) data
        h_i_aeff_data :     in, required, type=array/object
                            H I recombination coefficients
   
    :Examples:
       For example::
   
        >>> import pyequib
        >>> import atomneb
        >>> import os
        >>> base_dir = '../externals/atomneb/'
        >>> data_dir = os.path.join('atomic-data', 'chianti70')
        >>> atom_elj_file = os.path.join(base_dir,data_dir, 'AtomElj.fits')
        >>> atom_omij_file = os.path.join(base_dir,data_dir, 'AtomOmij.fits')
        >>> atom_aij_file = os.path.join(base_dir,data_dir, 'AtomAij.fits')
        >>> data_rc_dir = os.path.join('atomic-data-rc')
        >>> atom_rc_sh95_file= os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
        >>> atom='o'
        >>> ion='iii'
        >>> o_iii_elj=atomneb.read_elj(atom_elj_file, atom, ion, level_num=5) # read Energy Levels (Ej)
        >>> o_iii_omij=atomneb.read_omij(atom_omij_file, atom, ion) # read Collision Strengths (Omegaij)
        >>> o_iii_aij=atomneb.read_aij(atom_aij_file, atom, ion) # read Transition Probabilities (Aij)
        >>> atom='h'
        >>> ion='ii' # H I
        >>> hi_rc_data=atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)
        >>> h_i_aeff_data=hi_rc_data['aeff'][0]
        >>> temperature=np.float64(10000.0)
        >>> density=np.float64(5000.0)
        >>> atomic_levels='3,4/'
        >>> iobs5007=np.float64(1200.0)
        >>> abb5007=np.float64(0.0)
        >>> abb5007=pyequib.calc_abundance(temperature=temperature, density=density,
        >>>                        line_flux=iobs5007, atomic_levels=atomic_levels,
        >>>                        elj_data=o_iii_elj, omij_data=o_iii_omij,
        >>>                        aij_data=o_iii_aij, h_i_aeff_data=hi_rc_data['aeff'][0])
        >>> print('N(O^2+)/N(H+):', abb5007)
           N(O^2+)/N(H+):   0.00041256231
   
    :Categories:
      Abundance Analysis, Collisionally Excited Lines
   
    :Dirs:
     ./
         Main routines
   
    :Author:
      Ashkbiz Danehkar
   
    :Copyright:
      This library is released under a GNU General Public License.
   
    :Version:
      0.3.0
   
    :History:
        15/09/2013, A. Danehkar, Translated from FORTRAN to IDL code.
   
        20/10/2016, A. Danehkar, Replaced str2int with strnumber.
   
        20/10/2016, A. Danehkar, Replaced CFY, SPLMAT, and CFD with
             IDL function INTERPOL( /SPLINE).
   
        20/10/2016, A. Danehkar, Replaced LUSLV with IDL LAPACK function
                          LA_LINEAR_EQUATION.
   
        15/11/2016, A. Danehkar, Replaced LA_LINEAR_EQUATION (not work in GDL)
              with IDL function LUDC & LUSOL.
   
        19/11/2016, A. Danehkar, Replaced INTERPOL (not accurate) with
                       SPL_INIT & SPL_INTERP.
   
        20/11/2016, A. Danehkar, Made a new function calc_populations()
          for solving atomic level populations and separated it from
          calc_abundance(), calc_density() and calc_temperature().
   
        21/11/2016, A. Danehkar, Made a new function calc_emissivity()
                         for calculating line emissivities and separated it
                         from calc_abundance().
   
        10/03/2017, A. Danehkar, Integration with AtomNeb, now uses atomic data
                         input elj_data, omij_data, aij_data.
   
        12/06/2017, A. Danehkar, Cleaning the function, and remove unused varibales
                           from calc_abundance().

        03/10/2020, A. Danehkar, Transferred from IDL to Python.

    FORTRAN HISTORY:
   
        03/05/1981, I.D.Howarth,  Version 1.
   
        05/05/1981, I.D.Howarth,  Minibug fixed!
   
        07/05/1981, I.D.Howarth,  Now takes collision rates or strengths.
   
        03/08/1981, S.Adams,      Interpolates collision strengths.
   
        07/08/1981, S.Adams,      Input method changed.
   
        19/11/1984, R.E.S.Clegg,  SA files entombed in scratch disk. Logical
                                  filenames given to SA's data files.
   
        08/1995, D.P.Ruffle, Changed input file format. Increased matrices.
   
        02/1996, X.W.Liu,   Tidy up. SUBROUTINES SPLMAT, HGEN, CFY and CFD
                            modified such that matrix sizes (i.e. maximum
                            of Te and maximum no of levels) can now be cha
                            by modifying the parameters NDIM1, NDIM2 and N
                            in the Main program. EASY!
                            Now takes collision rates as well.
                            All variables are declared explicitly
                            Generate two extra files (ionpop.lis and ionra
                            of plain stream format for plotting.
   
        06/1996, C.J.Pritchet, Changed input data format for cases IBIG=1,2.
                            Fixed readin bug for IBIG=2 case.
                            Now reads reformatted upsilons (easier to see
                            and the 0 0 0 data end is excluded for these c
                            The A values have a different format for IBIG=.
   
        2006, B.Ercolano,   Converted to F90.
   """
   ahb = np.float64(0)
   
   h_planck = 6.62606957e-27 # erg.s
   c_speed = 2.99792458e10 # cm/s 
   
   if (temperature is not None) == 0:   
      print('Temperature is not set')
      return 0
   if (density is not None) == 0:   
      print('Density is not set')
      return 0
   if (elj_data is not None) == 0:   
      print('Energy Levels data (elj_data) are not set')
      return 0
   if (omij_data is not None) == 0:   
      print('Collision Strengths (omij_data) are not set')
      return 0
   if (aij_data is not None) == 0:   
      print('Transition Probabilities (aij_data) are not set')
      return 0
   if (h_i_aeff_data is not None) == 0:   
      print('H I recombination coefficients (h_i_aeff_data) are not set')
      return 0
   if (atomic_levels is not None) == 0:   
      print('Atomic levels (atomic_levels) are not given')
      return 0
   if (line_flux is not None) == 0:   
      print('Line flux intensity (line_flux) is not given')
      return 0
   if ((temperature <= 0.e0)|(density <= 0.e0)):
      print('temperature = ', temperature, ', density = ', density)
      return 0
   
   #  T4=temperature*1.0D-4
   #  AHB=3.036D-14*T4^(-0.87D0) ; Brocklehurt 1971; Aller (1984), Physics of Thermal Gaseous Nebulae, p. 76
   #  WAVHB=4861.33D ;4861.D0
   #  emissivity_Hbeta=AHB*h_Planck*c_Speed*1.e8/WAVHB ; N(H+) * N(e-) (erg/s)
   # emissivity_Hbeta=1.387D-25*T4^(-0.983D0)* 10.D0^(-0.0424D0/T4) ;  Brocklehurst (1971); Aller (1984)
   emissivity_hbeta = calc_emiss_h_beta(temperature=temperature, density=density, h_i_aeff_data=h_i_aeff_data)
   emissivity_all = np.float64(0.0)
   emissivity_all = calc_emissivity(temperature=temperature, density=density, atomic_levels=atomic_levels, elj_data=elj_data, omij_data=omij_data, aij_data=aij_data)
   
   if emissivity_all == 0:   
      print("cannot calculate emissivity")
      return 0
   abund = (emissivity_hbeta / emissivity_all) * (line_flux / 100.0)
   return abund

def print_ionic(temperature=None, density=None,
                elj_data=None, omij_data=None, aij_data=None,
                h_i_aeff_data=None,
                printemissivity=None, printpopulations=None, printcritdensity=None):
   """
       This function prints the atom's transitions information,
       atomic level populations, critical densities, and emissivities
       for given temperature and density.
   
    :Keywords:
        temperature   :   in, required, type=float
                          electron temperature
        density       :   in, required, type=float
                          electron density
        elj_data      :   in, required, type=array/object
                          energy levels (Ej) data
        omij_data     :   in, required, type=array/object
                          collision strengths (omega_ij) data
        aij_data      :   in, required, type=array/object
                          transition probabilities (Aij) data
        h_i_aeff_data :   in, type=array/object
                          H I recombination coefficients
        printEmissivity  :   in, type=boolean
                             Set for printing Emissivities
        printPopulations :   in, type=boolean
                             Set for printing Populations
        printCritDensity  :  in, type=boolean
                             Set for printing Critical Densities
   
    :Examples:
       For example::
   
        >>> import pyequib
        >>> import atomneb
        >>> import os
        >>> base_dir = '../externals/atomneb/'
        >>> data_dir = os.path.join('atomic-data', 'chianti70')
        >>> atom_elj_file = os.path.join(base_dir,data_dir, 'AtomElj.fits')
        >>> atom_omij_file = os.path.join(base_dir,data_dir, 'AtomOmij.fits')
        >>> atom_aij_file = os.path.join(base_dir,data_dir, 'AtomAij.fits')
        >>> data_rc_dir = os.path.join('atomic-data-rc')
        >>> atom_rc_sh95_file= os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
        >>> atom='o'
        >>> ion='iii'
        >>> o_iii_elj=atomneb.read_elj(atom_elj_file, atom, ion, level_num=5) # read Energy Levels (Ej)
        >>> o_iii_omij=atomneb.read_omij(atom_omij_file, atom, ion) # read Collision Strengths (Omegaij)
        >>> o_iii_aij=atomneb.read_aij(atom_aij_file, atom, ion) # read Transition Probabilities (Aij)
        >>> atom='h'
        >>> ion='ii' # H I
        >>> hi_rc_data=atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)
        >>> temperature=np.float64(10000.0)#
        >>> density = np.float64(1000.)
        >>> pyequib.print_ionic, temperature=temperature, density=density,
        >>>              elj_data=o_iii_elj, omij_data=o_iii_omij,
        >>>              aij_data=o_iii_aij, h_i_aeff_data=hi_rc_data['aeff'][0]
           Temperature =   10000.0 K
           Density =    1000.0 cm-3
   
           Level    Populations   Critical Densities
           Level 1:   3.063E-01   0.000E+00
           Level 2:   4.896E-01   4.908E+02
           Level 3:   2.041E-01   3.419E+03
           Level 4:   4.427E-05   6.853E+05
           Level 5:   2.985E-09   2.547E+07
   
            2.597E-05
                88.34um
               (2-->1)
            2.859E-22
   
            0.000E+00   9.632E-05
                32.66um      51.81um
               (3-->1)     (3-->2)
            0.000E+00   7.536E-22
   
            2.322E-06   6.791E-03   2.046E-02
              4932.60A    4960.29A    5008.24A
               (4-->1)     (4-->2)     (4-->3)
            4.140E-25   1.204E-21   3.593E-21
   
            0.000E+00   2.255E-01   6.998E-04   1.685E+00
              2315.58A    2321.67A    2332.12A    4364.45A
               (5-->1)     (5-->2)     (5-->3)     (5-->4)
            0.000E+00   5.759E-24   1.779E-26   2.289E-23
   
           H-beta emissivity: 1.237E-25 N(H+) Ne  [erg/s]
   
    :Categories:
      Plasma Diagnostics, Abundance Analysis, Collisionally Excited Lines
   
    :Dirs:
     ./
         Main routines
   
    :Author:
      Ashkbiz Danehkar
   
    :Copyright:
      This library is released under a GNU General Public License.
   
    :Version:
      0.3.0
   
    :History:
        04/03/2019, A. Danehkar, create the print_ionic() routine.

        03/10/2020, A. Danehkar, Transferred from IDL to Python.
   """
   h_planck = 6.62606957e-27 # erg s
   c_speed = 2.99792458e10 # cm/s
   
   if (temperature is not None) == 1:   
      print('Temperature = {0:9.1f} K'.format(temperature))
   if (density is not None) == 1:   
      print('Density = {0:9.1f} cm-3'.format(density))
   if (elj_data is not None) == 0:   
      print('elj_data is not set')
      return
   if (omij_data is not None) == 0:   
      print('omij_data is not set')
      return
   if (aij_data is not None) == 0:   
      print('aij_data is not set')
      return
   if (printemissivity is not None) == 0:   
      printemissivity = 1 # default
   if (printpopulations is not None) == 0:   
      printpopulations = 1 # default
   if (printcritdensity is not None) == 0:   
      printcritdensity = 1 # default
   level_num = len(elj_data)
   if (printpopulations is not None) == 1:   
      if ((temperature is not None) == 1 ) & ((density is not None) == 1):
         nlj = calc_populations(temperature=temperature, density=density, elj_data=elj_data, omij_data=omij_data, aij_data=aij_data, level_num=level_num)
      else:   
         print('calc_populations needs temperature and density.')
   if (printcritdensity is not None) == 1:   
      if (temperature is not None) == 1:   
         n_crit = calc_crit_density(temperature=temperature, elj_data=elj_data, omij_data=omij_data, aij_data=aij_data, level_num=level_num)
      else:   
         print('calc_crit_density needs temperature.')
   print('')
   if ((temperature is not None) == 1) & ((((printcritdensity is not None) == 1 ) | ((printpopulations is not None) == 1))):
      if (density is not None) == 1:   
         print('Level    Populations   Critical Densities')
      else:   
         print('Level    Critical Densities')
      for i in range(1, (level_num)+(1)):
         s = ''
         level_str = 'Level  {0:1d}:'.format(i)
         s = s + level_str + '   '
         if (density is not None) == 1:   
            nlj_str = '{0:9.3E}'.format(nlj[i - 1])
            s = s + nlj_str + '   '
         n_crit_str = '{0:9.3E}'.format(n_crit[i - 1])
         s = s + n_crit_str
         print(s)
      print('')
   if (printemissivity is not None) == 1:
      aij = aij_data['aij'][0]
      aij = aij.T
      elj = elj_data['ej']
      for i in range(2, (level_num)+(1)):
         aij_str = ''
         transition_str = ''
         wavelength_str = ''
         emissivity_str = ''
         for j in range(1, (i - 1)+(1)):
            eji = elj[i - 1] - elj[j - 1]
            wav = 1.e8 / eji
            aij_str = aij_str + '{0:10.3E}  '.format(aij[i - 1,j - 1])
            if wav < 10000:   
               wavelength_str = wavelength_str + '{0:10.2f}A '.format(wav)
            else:   
               wavelength_str = wavelength_str + '{0:10.2f}um '.format(wav * 1.e-4)
            transition_str = transition_str + '    ({0:1d}-->'.format(i)+ '{0:1d}) '.format(j)
            if ((temperature is not None) == 1) & ( (density is not None) == 1):
               emissivity_line = nlj[i - 1] * aij[i - 1,j - 1] * h_planck * c_speed * 1.e8 / (wav * density)
               emissivity_str = emissivity_str + '{0:10.3E}  '.format(emissivity_line)
            else:   
               emissivity_str = emissivity_str + ' '
         print(aij_str)
         print(wavelength_str)
         print(transition_str)
         print(emissivity_str)
         print('')
      if (h_i_aeff_data is not None) == 1:   
         if ((temperature is not None) == 1 ) & ( (density is not None) == 1):
            emissivity_hbeta = calc_emiss_h_beta(temperature=temperature, density=density, h_i_aeff_data=h_i_aeff_data)
            emissivity_str = '{0:10.3E}'.format(emissivity_hbeta)
            print("H-beta emissivity:", emissivity_str, " N(H+) Ne  [erg/s]")
         else:   
            print('gamma_hb_4861 needs temperature and density.')
   
   return

def get_omij_temp(temperature=None, omij_data=None, elj_data=None,
                  level_num=None, irats=None):
   """
        This function derives the effective collision strengths (Omij_T) from
        the collision strengths (omega_ij) data for the given temperature.
   
    :Returns:
       type=array/object. This function returns the effective collision strengths (Omij_T).
   
    :Keywords:
        temperature :   in, required, type=float
                        electron temperature
        omij_data   :   in, required, type=array/object
                        collision strengths (omega_ij) data
        level_num   :   in, type=int
                        Number of levels
        irats       :   in, type=int
                        Else Coll. rates = tabulated values * 10 ** irats
   
    :Examples:
       For example::
   
        >>> import pyequib
        >>> import atomneb
        >>> import os
        >>> base_dir = '../externals/atomneb/'
        >>> data_dir = os.path.join('atomic-data', 'chianti70')
        >>> atom_elj_file = os.path.join(base_dir,data_dir, 'AtomElj.fits')
        >>> atom_omij_file = os.path.join(base_dir,data_dir, 'AtomOmij.fits')
        >>> atom_aij_file = os.path.join(base_dir,data_dir, 'AtomAij.fits')
        >>> atom='s'
        >>> ion='ii'
        >>> s_ii_elj=atomneb.read_elj(atom_elj_file, atom, ion, level_num=5) # read Energy Levels (Ej)
        >>> s_ii_omij=atomneb.read_omij(atom_omij_file, atom, ion) # read Collision Strengths (Omegaij)
        >>> s_ii_aij=atomneb.read_aij(atom_aij_file, atom, ion) # read Transition Probabilities (Aij)\
        >>> temperature=np.float64(10000.0)#
        >>> omij_t=pyequib.get_omij_temp(temperature=temperature, omij_data=s_ii_omij)
        >>> print('Effective Collision Strengths: ')
        >>> print(omij_t)
           Effective Collision Strengths:
           0.0000000       0.0000000       0.0000000       0.0000000       0.0000000
           2.7800000       0.0000000       0.0000000       0.0000000       0.0000000
           4.1600000       7.4600000       0.0000000       0.0000000       0.0000000
           1.1700000       1.8000000       2.2000000       0.0000000       0.0000000
           2.3500000       3.0000000       4.9900000       2.7100000       0.0000000
   
    :Categories:
      Plasma Diagnostics, Abundance Analysis, Collisionally Excited Lines
   
    :Dirs:
     ./
         Subroutines
   
    :Author:
      Ashkbiz Danehkar
   
    :Copyright:
      This library is released under a GNU General Public License.
   
    :Version:
      0.3.0
   
    :History:
        04/03/2019, A. Danehkar, create the get_omij_temp() routine.

        03/10/2020, A. Danehkar, Transferred from IDL to Python.
   """
   h_planck = 4.13566766225e-15 # eV.s #6.62606957e-27 # erg.s
   c_speed = 2.99792458e10 # cm/s
   k_b = 8.617330350e-5 # eV/K # 1.3806485279e-16 # erg/K
   
   if (temperature is not None) == 0:   
      print('Temperature is not set')
      return 0
   if (omij_data is not None) == 0:   
      print('omij_data is not set')
      return 0
   if (level_num is not None) == 0:   
      if (elj_data is not None) == 0:   
         level_num = max([max(omij_data['level1'][:]), max(omij_data['level2'][:])])
      else:   
         level_num = len(elj_data)
   if (irats is not None) == 0:   
      irats = 0
   else:   
      if (irats != 0):
         if (elj_data is not None) == 0:
            print('elj_data is not set. It is required for irats')
            return 0
         ej = elj_data['ej'] # Energy Levels (E_j) in cm-1
   t_log = np.float64(np.log10(temperature))
   t_num = len(omij_data['strength'][0]) # Number of temperature intervals
   t_lin_list = np.float64(omij_data['strength'][0])
   t_log_list = np.log10(t_lin_list) # temperature intervals (array)
   omij_num = len(omij_data)
   omij_t = np.zeros((level_num, level_num))
   for k in range(1, (omij_num - 1)+(1)):
      i = omij_data['level1'][k]
      j = omij_data['level2'][k]
      if (i <= level_num) & (j <= level_num):
         qj = np.float64(omij_data['strength'][k])
         if (irats != 0):   
            d_e =  np.float64(ej[j - 1] - ej[i - 2]) * h_planck * c_speed # delta Energy in eV; convert from cm-1 to eV
            # Calculate the Boltzmann factor
            exp_de_kt = np.exp(-d_e / (k_b * temperature)) # Maxwell-Boltzmann distribution
            qj = qj / exp_de_kt #Take out the exp. before interpolation
         if (t_num == 1):   
            omij_t[j - 1,i - 1] = qj
         else:   
            if (t_num == 2):   
               omij_t[j - 1,i - 1] = qj[0] + (qj[1] - qj[0]) / (t_log_list[1] - t_log_list[0]) * (t_log - t_log_list[0])
            else:   
               #Qj_T=interpol(Qj, T_log_list, T_log, /SPLINE)
               # Calculate interpolating cubic spline
               #qj_2 = spl_init(t_log_list, qj)
               # Calculate the interpolated Omij_T values corresponding to T_log
               # Obtain the effective collision strengths Omij_T
               #qj_t = spl_interp(t_log_list, qj, qj_2, t_log, double=True)
               # qj_t=np.interp(t_log, t_log_list, qj)
               interpfunc = interpolate.interp1d(t_log_list, qj, kind='cubic')
               qj_t = interpfunc(t_log)
               omij_t[j - 1,i - 1] = qj_t
   return omij_t

def calc_emiss_h_beta(temperature=None, density=None, h_i_aeff_data=None):
    """
         This function calculates the emissivity for H_beta 4861A
         Emis(Hbeta)= 4pi j(HBeta 4861 A)/Np Ne) for the given temperature and density
         by using the helium emissivities from
         Storey & Hummer, 1995MNRAS.272...41S.

     :Private:

     :Returns:
        type=double. This function returns the H beta emissivity 4pi j(HBeta 4861)/Np Ne).

     :Keywords:
         temperature     :   in, required, type=float
                             electron temperature
         density         :   in, required, type=float
                             electron density
         h_i_aeff_data   :   in, required, type=array/object
                             H I recombination coefficients

     :Categories:
       Abundance Analysis, Recombination Lines, Emissivity

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.3.0

     :History:
         Based on H I emissivities
         from Storey & Hummer, 1995MNRAS.272...41S.

         25/08/2012, A. Danehkar, IDL code written.

         11/03/2017, A. Danehkar, Integration with AtomNeb.

         10/07/2019, A. Danehkar, Change from logarithmic to linear

         03/10/2020, A. Danehkar, Transferred from IDL to Python.
    """
    if (temperature is not None) == 0:
        print('Temperature is not set')
        return 0
    if (density is not None) == 0:
        print('Density is not set')
        return 0
    if (h_i_aeff_data is not None) == 0:
        print('H I recombination coefficients (h_i_aeff_data) are not set')
        return 0

    # h_a_col= find_aeff_sh95_column(3, 2)
    linenum = find_aeff_sh95_column(4, 2, 25)

    teh2 = np.float64(temperature)
    neh2 = np.float64(density)
    line1 = np.int32(linenum - 1)
    emissivity = np.float64(0.0)

    h_i_ems = np.zeros((10, 13))
    temp1 = np.zeros(302)
    temp_grid = np.array([500., 1000., 3000., 5000., 7500., 10000., 12500., 15000., 20000., 30000.])

    nlines = 130
    h_i_aeff=h_i_aeff_data
    h_i_aeff=h_i_aeff.T
    for i in range(0, nlines):
        temp1 = h_i_aeff[:, i]
        loc1=(np.where(temp_grid == temp1[1]))
        loc1 = np.asarray(loc1[0])
        tpos = np.int32(loc1)
        npos = np.int32(round(np.log10(temp1[0]) - 2))
        h_i_ems[tpos, npos] = temp1[line1]  # temp[2:45]

    # restrict to the density & temperature ranges to 1995MNRAS.272...41S
    if (neh2 < 1.1e2):
        neh2 = 1.1e2
    if (neh2 > 1.e14):
        neh2 = 1.e14
    if (teh2 < 550.):
        teh2 = 550.
    if (teh2 > 30000.):
        teh2 = 30000.

    # get logarithmic density
    dens_log = np.log10(neh2)

    dens_grid = np.float64(np.arange(13) + 2)
    hb_emissivity = np.float64(0.0)
    # Bilinearly interpolate density & temperature
    # emiss_log =_interp2d(hi_ems1, temp_grid, dens_grid, TEh2, dens_log, [101,101], /cubic, /quintic)#, /trigrid) not work on GDL
    hb_emissivity = interp_2d(h_i_ems[:, :], teh2, dens_log, temp_grid, dens_grid)

    # ems_log = np.log10(hr_tmp/np.float64(4861.33/1.98648E-08))
    # h_beta_emissivity = 10.0^ems_log
    # h_beta_emissivity_log=np.log10(h_beta_emissivity(temp, density))

    return hb_emissivity

def find_aeff_sh95_column(lo_lev, hi_lev, lev_num):
    """
        This function locates and returns the data location
        of the given low energy level, high energy level,
        and the level number
        within the database of H I emissivities given by
        from Storey & Hummer, 1995MNRAS.272...41S.

    :Private:

    :Returns:
       type=double. This function returns the data location .

    :Params:
        lo_lev  :  in, required, type=float
                      low energy level

        hi_lev  :  in, required, type=float
                      high energy level

        lev_num :  in, required, type=float
                      level number

    :Author:
      Ashkbiz Danehkar

    :Copyright:
      This library is released under a GNU General Public License.

    :Version:
      0.3.0

    :History:
        Based on H I emissivities
        from Storey & Hummer, 1995MNRAS.272...41S.

        25/08/2012, A. Danehkar, IDL code written.

        11/03/2017, A. Danehkar, Integration with AtomNeb.
 
        03/10/2020, A. Danehkar, Transferred from IDL to Python.
    """
    # lev_num=25
    count = 2
    # lo_lev=3
    # hi_lev=2
    for k in range(lev_num, 0, -1):
        for l in range(1, k):
            count = count + 1
            if ((k == lo_lev) & (l == hi_lev)):
                return count
    return 0

def interp_2d(tab_in, x, y, x_tab, y_tab, xlog=None, ylog=None):
    """
       function interp_2d,tab_in,x,y,x_tab,y_tab,[xlog=xlog],[ylog=ylog]
       return an interpolated 2 Dim less than tab_in (no extrapolation)
       tab_in: fltarr([*,*,...],n_x,n_y)
       x,y: values where we want to interpolate
       x_tab,y_tab: fltarr(n_x),fltarr(n_y) : x and y vectors, must be
          sorted, no necessarily in increasing order.
       x_log and y_log: allows log interpolation
       C. Morisset (LAS, Marseille, 2000)

       03/10/2020, A. Danehkar, Transferred from IDL to Python.
    """
    # ON_ERROR, 2

    size_tab_in = [len(tab_in), len(tab_in[0])]
    n_dim_tab = 1 #size_tab_in[0]
    n_x = len(x_tab)
    n_y = len(y_tab)
    if (n_x != size_tab_in[n_dim_tab - 1]) | (n_y != size_tab_in[n_dim_tab]):
        print('Dimension of x_tab or y_tab incompatible with tab_in')
    tab_tmp = tab_in #reform(tab_in, size_tab_in[n_dim_tab + 2] / n_x / n_y, n_x, n_y)

    # No extrapolation:
    if ((x > max(x_tab))| (x < min(x_tab)) | ( y > max(y_tab)) | (y < min(y_tab))):
        print('No extrapollation')

    i_x = np.amax(np.where(x_tab < x))
    i_y = np.amax(np.where(y_tab < y))

    if x_tab[0] > x_tab[1]:
        x_incr = -1
    else:
        x_incr = 1
    if y_tab[0] > y_tab[1]:
        y_incr = -1
    else:
        y_incr = 1

    if x_incr == 1:
        i_x = np.amax(np.where(x_tab < x))
    else:
        i_x = np.amax(np.where(x_tab > x))
    if y_incr == 1:
        i_y = np.amax(np.where(y_tab < y))
    else:
        i_y = np.amax(np.where(y_tab > y))

    if (xlog is not None):
        f_x = 1. - (np.log10(x) - np.log10(x_tab[i_x])) / (np.log10(x_tab[i_x + 1]) - np.log10(x_tab[i_x]))
    else:
        f_x = 1. - (x - x_tab[i_x]) / (x_tab[i_x + 1] - x_tab[i_x])
    if (ylog is not None):
        f_y = 1. - (np.log10(y) - np.log10(y_tab[i_y])) / (np.log10(y_tab[i_y + 1])
                                                           - np.log10(y_tab[i_y]))
    else:
        f_y = 1. - (y - y_tab[i_y]) / (y_tab[i_y + 1] - y_tab[i_y])

    tab_out = tab_tmp[i_x, i_y] * f_x * f_y + tab_tmp[i_x + 1, i_y] * (1. - f_x) * f_y \
              + tab_tmp[i_x, i_y + 1] * f_x * (1. - f_y) \
              + tab_tmp[i_x + 1, i_y + 1] * (1. - f_x) * (1. - f_y)

    if n_dim_tab > 2:
        tab_out = tab_out #reform(tab_out, size_tab_in[1:(n_dim_tab - 2) + 1], overwrite=True)
    return tab_out

def calc_emiss_he_i_rl(temperature=None, density=None,
                       linenum=None, he_i_aeff_data=None):
   """
        This function calculates the emissivity
        for the given wavelength of He I recombination line
        by using the recombination coefficients from Porter et al.
        2012MNRAS.425L..28P.
   
    :Returns:
       type=double. This function returns the line emissivity.
   
    :Keywords:
        temperature    :    in, required, type=float
                            electron temperature
        density        :    in, required, type=float
                            electron density
        linenum        :    in, required, type=int
                            Line Number for Wavelength
   
                            Wavelength=4120.84:linenum=7,
   
                            Wavelength=4387.93: linenum=8,
   
                            Wavelength=4437.55: linenum=9,
   
                            Wavelength=4471.50: linenum=10,
   
                            Wavelength=4921.93: linenum=12,
   
                            Wavelength=5015.68: linenum=13,
   
                            Wavelength=5047.74: linenum=14,
   
                            Wavelength=5875.66: linenum=15,
   
                            Wavelength=6678.16: linenum=16,
   
                            Wavelength=7065.25: linenum=17,
   
                            Wavelength=7281.35: linenum=18.
   
        line_flux      :    in, required, type=float
                            line flux intensity
        he_i_aeff_data :    in, required, type=array/object
                            He I recombination coefficients
   
    :Examples:
       For example::
   
        >>> import pyequib
        >>> import atomneb
        >>> import os
        >>> base_dir = '../externals/atomneb/'
        >>> data_rc_dir = os.path.join('atomic-data-rc')
        >>> atom_rc_he_i_file= filepath('rc_he_ii_PFSD12.fits', root_dir=base_dir, subdir=data_rc_dir )
        >>> atom_rc_sh95_file= os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
        >>>
        >>> atom='he'
        >>> ion='ii' # He I
        >>> he_i_rc_data=atomneb.read_aeff_he_i_pfsd12(atom_rc_he_i_file, atom, ion)
        >>> he_i_aeff_data=he_i_rc_data['aeff'][0]
        >>> temperature=np.float64(10000.0)
        >>> density=np.float64(5000.0)
        >>> linenum=10# 4471.50
        >>> emiss_he_i=pyequib.calc_emiss_he_i_rl(temperature=temperature, density=density,
        >>>                                linenum=linenum,
        >>>                                he_i_aeff_data=he_i_aeff_data)
        >>> print('Emissivity:', emiss_he_i)
           Emissivity:   6.3822830e-26
   
    :Categories:
      Abundance Analysis, Recombination Lines, Emissivity
   
    :Dirs:
     ./
         Main routines
   
    :Author:
      Ashkbiz Danehkar
   
    :Copyright:
      This library is released under a GNU General Public License.
   
    :Version:
      0.3.0
   
    :History:
        Based on improved He I emissivities in the case B
        from Porter et al. 2012MNRAS.425L..28P
   
        15/12/2013, A. Danehkar, IDL code written.
   
        20/03/2017, A. Danehkar, Integration with AtomNeb.
   
        10/07/2019, A. Danehkar, Made a new function calc_emiss_he_i_rl()
                         for calculating line emissivities and separated it
                         from calc_abund_he_i_rl().

        03/10/2020, A. Danehkar, Transferred from IDL to Python.
   """
   if (temperature is not None) == 0:   
      print('Temperature is not set')
      return 0
   if (density is not None) == 0:   
      print('Density is not set')
      return 0
   if (he_i_aeff_data is not None) == 0:   
      print('He I recombination coefficients (he_i_aeff_data) are not set')
      return 0
   if (linenum is not None) == 0:   
      print('Line Number for Wavelength is not given')
      return 0
   if ((temperature <= 0.e0) | (density <= 0.e0)):
      print('temperature = ', temperature, ', density = ', density)
      return 0
   
   teh2 = np.float64(temperature)
   neh2 = np.float64(density)
   line1 = np.int32(linenum - 1)
   emissivity =  np.float64(0.0)
   
   # hei_ems=read_porter()
   hei_ems = np.zeros((21, 14)) #(21,14,44)
   temp1 = np.zeros(46)
   
   nlines = 294
   he_i_aeff = he_i_aeff_data
   he_i_aeff = he_i_aeff.T
   for i in range(0, nlines):
      temp1 = he_i_aeff[:,i]
      tpos = np.int32(round((temp1[0] / 1000) - 5))
      npos = np.int32(round(np.log10(temp1[1]) - 1))
      hei_ems[tpos,npos] = temp1[line1 + 2]#temp[2:45]
   
   # restrict to the density & temperature ranges to 2012MNRAS.425L..28P
   if (neh2 < 1.e1):   
      neh2 = 1.e1
   if (neh2 > 1.e14):   
      neh2 = 1.e14
   if (teh2 < 5000):   
      teh2 = 5000.
   if (teh2 > 25000):   
      teh2 = 25000.
   
   # get logarithmic density
   dens_log = np.log10(neh2)
   
   dens_grid = np.float64(np.arange(14) + 1)
   temp_grid = np.float64(1000 * (np.arange(21) + 5))
   
   hei_ems1 = hei_ems[:,:]
   # Bilinearly interpolate density & temperature
   # emiss_log =_interp2d(hei_ems1, temp_grid, dens_grid, TEh2, dens_log, [101,101], /cubic, /quintic)#, /trigrid) not work on GDL
   emiss_log = interp_2d(hei_ems1, teh2, dens_log, temp_grid, dens_grid)
   
   # wavl=he_i_aeff_data_Wavelength[line1]
   emissivity = 10.e0 ** (emiss_log)
   
   return emissivity

def calc_emiss_he_ii_rl(temperature=None, density=None,
                        he_ii_aeff_data=None):
   """
        This functioncalculates the emissivity
        for the He II recombination line 4686 A
        by using the helium emissivities from
        Storey & Hummer, 1995MNRAS.272...41S.
   
    :Returns:
       type=double. This function returns the line emissivity.
   
    :Keywords:
        temperature     :   in, required, type=float
                            electron temperature
        density         :   in, required, type=float
                            electron density
        he_ii_aeff_data :   in, required, type=array/object
                            He II recombination coefficients
   
    :Examples:
       For example::
   
        >>> import pyequib
        >>> import atomneb
        >>> import os
        >>> base_dir = '../externals/atomneb/'
        >>> data_rc_dir = os.path.join('atomic-data-rc')
        >>> atom_rc_he_i_file= os.path.join(base_dir,data_rc_dir, 'rc_he_ii_PFSD12.fits')
        >>> atom_rc_sh95_file= os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
        >>>
        >>> atom='he'
        >>> ion='iii' # He II
        >>> he_ii_rc_data=atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)
        >>> he_ii_aeff_data=he_ii_rc_data['aeff'][0]
        >>> temperature=np.float64(10000.0)
        >>> density=np.float64(5000.0)
        >>> he_ii_4686_flux = 135.833
        >>> emiss_he_ii=pyequib.calc_emiss_he_ii_rl(temperature=temperature, density=density,
        >>>                                 he_ii_aeff_data=he_ii_aeff_data)
        >>> print('Emissivity:', emiss_he_ii)
           Emissivity:   1.4989134e-24
   
    :Categories:
      Abundance Analysis, Recombination Lines, Emissivity
   
    :Dirs:
     ./
         Main routines
   
    :Author:
      Ashkbiz Danehkar
   
    :Copyright:
      This library is released under a GNU General Public License.
   
    :Version:
      0.3.0
   
    :History:
        Based on He II emissivities
        from Storey & Hummer, 1995MNRAS.272...41S.
   
        15/12/2013, A. Danehkar, IDL code written.
   
        02/04/2017, A. Danehkar, Integration with AtomNeb.
   
        10/07/2019, A. Danehkar, Change from logarithmic to linear
   
        10/07/2019, A. Danehkar, Made a new function calc_emiss_he_ii_rl()
                         for calculating line emissivities and separated it
                         from calc_abund_he_ii_rl().

        03/10/2020, A. Danehkar, Transferred from IDL to Python.
   """
   if (temperature is not None) == 0:   
      print('Temperature is not set')
      return 0
   if (density is not None) == 0:   
      print('Density is not set')
      return 0
   if (he_ii_aeff_data is not None) == 0:   
      print('He II recombination coefficients (he_ii_aeff_data) are not set')
      return 0
   if ((temperature <= 0.e0) | (density <= 0.e0)):
      print('temperature = ', temperature, ', density = ', density)
      return 0
   
   #h_a_col= find_aeff_sh95_column(3, 2)
   linenum = find_aeff_sh95_column(4, 3, 25)
   
   teh2 = np.float64(temperature)
   neh2 = np.float64(density)
   line1 = np.int32(linenum - 1)
   emissivity = np.float64(0.0)
   
   heii_ems = np.zeros((12, 13))
   temp1 = np.zeros(302)
   temp_grid = np.array([500., 1000., 3000., 5000., 7500., 10000., 12500., 15000., 20000., 30000., 50000., 100000.])
   
   nlines = 156

   he_ii_aeff = he_ii_aeff_data
   he_ii_aeff = he_ii_aeff.T
   for i in range(0, nlines):
       temp1 = he_ii_aeff[:, i]
       loc1 = (np.where(temp_grid == temp1[1]))
       loc1 = np.asarray(loc1[0])
       tpos = np.int32(loc1)
       npos = np.int32(round(np.log10(temp1[0]) - 2))
       heii_ems[tpos, npos] = temp1[line1]  # temp[2:45]
   
   # restrict to the density & temperature ranges to 2012MNRAS.425L..28P
   if (neh2 < 1.e2):   
      neh2 = 1.e2
   if (neh2 > 1.e14):   
      neh2 = 1.e14
   if (teh2 < 500.):   
      teh2 = 500.
   if (teh2 > 100000.):   
      teh2 = 100000.
   
   # get logarithmic density
   dens_log = np.log10(neh2)
   
   dens_grid = np.float64(np.arange(13) + 2)
   
   heii_ems1 = heii_ems[:,:]
   # Bilinearly interpolate density & temperature
   # emissivity =_interp2d(heii_ems1, temp_grid, dens_grid, TEh2, dens_log, [101,101], /cubic, /quintic)#, /trigrid) not work on GDL
   emissivity = interp_2d(heii_ems1, teh2, dens_log, temp_grid, dens_grid)
   # emissivity_log = np.log10(emissivity)
   
   return emissivity

def calc_emiss_c_ii_rl(temperature=None, density=None,
                       wavelength=None, c_ii_rc_data=None):
    """
        This function calculates the emissivity
        for the given wavelength of C II recombination line
        by using the recombination coefficients from
        from Davey et al. (2000) 2000A&AS..142...85D.

    :Returns:
       type=double. This function returns the line emissivity.

    :Keywords:
        temperature   :     in, required, type=float
                            electron temperature
        density       :     in, required, type=float
                            electron density
        wavelength    :     in, required, type=float
                            Line Wavelength in Angstrom
        c_ii_rc_data  :     in, required, type=array/object
                            C II recombination coefficients

    :Examples:
       For example::

        >>> import pyequib
        >>> import atomneb
        >>> import os
        >>> base_dir = '../externals/atomneb/'
        >>> data_rc_dir = os.path.join('atomic-data-rc')
        >>> atom_rc_all_file= os.path.join(base_dir,data_rc_dir, 'rc_collection.fits')
        >>> atom_rc_sh95_file= os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
        >>>
        >>> atom='c'
        >>> ion='iii' # C II
        >>> c_ii_rc_data=atomneb.read_aeff_collection(atom_rc_all_file, atom, ion)
        >>> temperature=np.float64(10000.0)
        >>> density=np.float64(5000.0)
        >>> wavelength=6151.43
        >>> emiss_c_ii=pyequib.calc_emiss_c_ii_rl(temperature=temperature, density=density,
        >>>                               wavelength=wavelength,
        >>>                               c_ii_rc_data=c_ii_rc_data)
        >>> print('Emissivity:', emiss_c_ii)
           Emissivity:   5.4719511e-26

    :Categories:
      Abundance Analysis, Recombination Lines, Emissivity

    :Dirs:
     ./
         Main routines

    :Author:
      Ashkbiz Danehkar

    :Copyright:
      This library is released under a GNU General Public License.

    :Version:
      0.3.0

    :History:
        Based on recombination coefficients for C II lines from
        Davey et al. 2000A&AS..142...85D.

        Adopted from MOCASSIN, Ercolano et al. 2005MNRAS.362.1038E.

        02/2003, Yong Zhang, added to MOCASSIN.

        10/05/2013, A. Danehkar, Translated to IDL code.

        15/04/2017, A. Danehkar, Integration with AtomNeb.

        10/07/2019, A. Danehkar, Made a new function calc_emiss_c_ii_rl()
                         for calculating line emissivities and separated it
                         from calc_abund_c_ii_rl().        
        
        03/10/2020, A. Danehkar, Transferred from IDL to Python.
    """
   # ciiRLstructure ={Wave:np.float64(0.0), Int:np.float64(0.0), Obs:np.float64(0.0), Abundance:np.float64(0.0)}
    h_planck = 6.62606957e-27 # erg s
    c_speed = 2.99792458e10 # cm/s
   
    if (temperature is not None) == 0:
        print('Temperature is not set')
        return 0
    if (density is not None) == 0:
        print('Density is not set')
        return 0
    if (c_ii_rc_data is not None) == 0:
        print('C II recombination coefficients (c_ii_rc_data) are not set')
        return 0
    if (wavelength is not None) == 0:
        print('Wavelength is not given')
        return 0
    if ((temperature <= 0.e0) | (density <= 0.e0)):
        print('temperature = ', temperature, ', density = ', density)
        return 0

    #c_ii_rc = np.asarray(c_ii_rc_data)

    lamb = np.float64(0.0)
    a = np.float64(0.0)
    b = np.float64(0.0)
    c = np.float64(0.0)
    d = np.float64(0.0)
    f = np.float64(0.0)
    aeff = np.float64(0.0)
    br = np.float64(1.0)
    temp4 = temperature / 10000.0
    loc1 = np.where(abs(c_ii_rc_data['wavelength'] - wavelength) <= 1.5)
    loc1=np.asarray(loc1[0])
    temp2 = len(loc1)
    if temp2 != 1:
        wavelength_min = np.amin(c_ii_rc_data['wavelength'][loc1])
        loc1 = np.where(c_ii_rc_data['wavelength'] == wavelength_min)
        loc1 = np.asarray(loc1[0])
    lamb = np.float64(c_ii_rc_data['wavelength'][loc1])
    a = np.float64(c_ii_rc_data['a'][loc1])
    b = np.float64(c_ii_rc_data['b'][loc1])
    c = np.float64(c_ii_rc_data['c'][loc1])
    d = np.float64(c_ii_rc_data['d'][loc1])
    f = np.float64(c_ii_rc_data['f'][loc1])
    aeff = 1.0e-14 * (a * (temp4 ** f))
    aeff = aeff * (1. + (b * (1. - temp4)) + (c * ((1. - temp4) ** 2)) + (d * ((1. - temp4) ** 3)))
    #ciiRLs_Int = 100.0*(aeff/hbeta_aeff)*br*(4861.33/lamb)*abund
    #abund=line_flux/ciiiRLs_Int
    emissivity = (np.float64(aeff * br) / np.float64(lamb)) * np.float64(h_planck * c_speed * 1.e8)
   
    return emissivity

def calc_emiss_c_iii_rl(temperature=None, density=None,
                        wavelength=None, c_iii_rc_data=None):
   """
         This function calculates the emissivity
         for the given wavelength of C III recombination line
         by using the recombination coefficients from
         Pequignot et al. 1991A&A...251..680P.

     :Returns:
        type=double. This function returns the line emissivity.

     :Keywords:
         temperature   :     in, required, type=float
                             electron temperature
         density       :     in, required, type=float
                             electron density
         wavelength    :     in, required, type=float
                             Line Wavelength in Angstrom
         c_iii_rc_data :     in, required, type=array/object
                             C III recombination coefficients

    :Examples:
       For example::

         >>> import pyequib
         >>> import atomneb
         >>> import os
         >>> base_dir = '../externals/atomneb/'
         >>> data_rc_dir = os.path.join('atomic-data-rc')
         >>> atom_rc_ppb91_file=os.path.join(base_dir,data_rc_dir, 'rc_PPB91.fits')
         >>> atom_rc_sh95_file= os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
         >>>
         >>> atom='c'
         >>> ion='iv' # C III
         >>> c_iii_rc_data=atomneb.read_aeff_ppb91(atom_rc_ppb91_file, atom, ion)
         >>> temperature=np.float64(10000.0)
         >>> density=np.float64(5000.0)
         >>> wavelength=4647.42
         >>> emiss_c_iii=pyequib.calc_emiss_c_iii_rl(temperature=temperature, density=density,
         >>>                                 wavelength=wavelength,
         >>>                                 c_iii_rc_data=c_iii_rc_data)
         >>> print('Emissivity:', emiss_c_iii)
            Emissivity:   7.5749632e-25

     :Categories:
       Abundance Analysis, Recombination Lines, Emissivity

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.3.0

     :History:
         Based on effective radiative recombination coefficients for C III lines from
         Pequignot, Petitjean, Boisson, C. 1991A&A...251..680P.

         18/05/2013, A. Danehkar, Translated to IDL code.

         06/04/2017, A. Danehkar, Integration with AtomNeb.

         10/07/2019, A. Danehkar, Made a new function calc_emiss_c_iii_rl()
                          for calculating line emissivities and separated it
                          from calc_abund_c_iii_rl().

         03/10/2020, A. Danehkar, Transferred from IDL to Python.
   """
   # ciiiRLstructure ={Wave:np.float64(0.0), Int:np.float64(0.0), Obs:np.float64(0.0), Abundance:np.float64(0.0)}
   h_planck = 6.62606957e-27 # erg s
   c_speed = 2.99792458e10 # cm/s 
   
   if (temperature is not None) == 0:   
      print('Temperature is not set')
      return 0
   if (density is not None) == 0:   
      print('Density is not set')
      return 0
   if (c_iii_rc_data is not None) == 0:   
      print('C III recombination coefficients (c_iii_rc_data) are not set')
      return 0
   if (wavelength is not None) == 0:   
      print('Wavelength is not given')
      return 0
   if ((temperature <= 0.e0) | (density <= 0.e0)):
      print('temperature = ', temperature, ', density = ', density)
      return 0

   lamb = np.float64(0.0)
   a = np.float64(0.0)
   b = np.float64(0.0)
   c = np.float64(0.0)
   d = np.float64(0.0)
   br = np.float64(0.0)
   aeff = np.float64(0.0)
   ion = ''
   
   z = 3.0 # ion level c^3+
   # equation (1) in 1991A&A...251..680P
   temp4 = 1.0e-4 * temperature / z ** 2
   loc1 = np.where(abs(c_iii_rc_data['wavelength'] - wavelength) <= 0.01)
   loc1=np.asarray(loc1[0])
   temp2 = len(loc1)
   if temp2 != 1:
      wavelength_min = np.amin(c_iii_rc_data['wavelength'][loc1])
      loc1 = np.where(c_iii_rc_data['wavelength'] == wavelength_min)
      loc1 = np.asarray(loc1[0])
   lamb = np.float64(c_iii_rc_data['wavelength'][loc1])
   a = np.float64(c_iii_rc_data['a'][loc1])
   b = np.float64(c_iii_rc_data['b'][loc1])
   c = np.float64(c_iii_rc_data['c'][loc1])
   d = np.float64(c_iii_rc_data['d'][loc1])
   br = np.float64(c_iii_rc_data['br'][loc1])
   # equation (1) in 1991A&A...251..680P
   aeff = 1.0e-13 * z * br
   aeff = aeff * (a * (temp4 ** b)) / (1. + c * (temp4 ** d))
   #ciiiRLs_Int = 100.0 * (aeff/hbeta_aeff) * (4861.33/lamb) * abund
   #abund=line_flux/ciiiRLs_Int
   emissivity = (np.float64(aeff) / np.float64(lamb)) * np.float64(h_planck * c_speed * 1.e8)
   
   return emissivity

def calc_emiss_n_ii_rl(temperature=None, density=None,
                       wavelength=None,
                       n_ii_rc_br=None, n_ii_rc_data=None):
   """
        This function calculates the emissivity
        for the given wavelength of N II recombination line
        by using the recombination coefficients from
        Escalante & Victor 1990ApJS...73..513E.

    :Returns:
       type=double. This function returns the line emissivity.

    :Keywords:
        temperature   :     in, required, type=float
                            electron temperature
        density       :     in, required, type=float
                            electron density
        wavelength    :     in, required, type=float
                            Line Wavelength in Angstrom
        n_ii_rc_br    :     in, required, type=array/object
                            N II branching ratios (Br)
        n_ii_rc_data  :     in, required, type=array/object
                            N II recombination coefficients

    :Examples:
       For example::

        >>> import pyequib
        >>> import atomneb
        >>> import os
        >>> base_dir = '../externals/atomneb/'
        >>> data_rc_dir = os.path.join('atomic-data-rc')
        >>> atom_rc_all_file= os.path.join(base_dir,data_rc_dir, 'rc_collection.fits')
        >>> atom_rc_sh95_file= os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
        >>> atom='h'
        >>> ion='ii' # H I
        >>> h_i_rc_data=atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)
        >>> h_i_aeff_data=h_i_rc_data['aeff'][0]
        >>> atom='n'
        >>> ion='iii' # N II
        >>> n_ii_rc_data=atomneb.read_aeff_collection(atom_rc_all_file, atom, ion)
        >>> n_ii_rc_data_br=atomneb.read_aeff_collection(atom_rc_all_file, atom, ion, br=True)
        >>> temperature=np.float64(10000.0)
        >>> density=np.float64(5000.0)
        >>> wavelength=4442.02
        >>> emiss_n_ii=pyequib.calc_emiss_n_ii_rl(temperature=temperature, density=density,
        >>>                               wavelength=wavelength,
        >>>                               n_ii_rc_br=n_ii_rc_data_br, n_ii_rc_data=n_ii_rc_data,
        >>>                               h_i_aeff_data=h_i_aeff_data)
        >>> print('Emissivity:', emiss_n_ii)
           Emissivity:   3.0397397e-26

    :Categories:
      Abundance Analysis, Recombination Lines, Emissivity

    :Dirs:
     ./
         Main routines

    :Author:
      Ashkbiz Danehkar

    :Copyright:
      This library is released under a GNU General Public License.

    :Version:
      0.3.0

    :History:
        Based on Effective recombination coefficients for N II lines from
        Escalante & Victor 1990ApJS...73..513E.

        Adopted from MIDAS Rnii script written by X.W.Liu.

        Revised based on scripts by Yong Zhang added to MOCASSIN, 02/2003
                          Ercolano et al. 2005MNRAS.362.1038E.

        10/05/2013, A. Danehkar, Translated to IDL code.

        25/04/2017, A. Danehkar, Integration with AtomNeb.

        10/07/2019, A. Danehkar, Made a new function calc_emiss_n_ii_rl()
                         for calculating line emissivities and separated it
                         from calc_abund_n_ii_rl().

        03/10/2020, A. Danehkar, Transferred from IDL to Python.
   """

#  niiRLstructure ={Wave:np.float64(0.0), #REAL*8
#              Int:np.float64(0.0),
#              Obs:np.float64(0.0),
#              abundance:np.float64(0.0),
#              g1:long(0), #INTEGER
#              g2:long(0), #INTEGER
#              Mult1:'', #CHARACTER*7
#              Term1:'', #CHARACTER*9
#              Term2:'' #CHARACTER*9
#              }
   h_planck = np.float64(6.62606957e-27) # erg s
   c_speed = np.float64(2.99792458e10) # cm/s
   
   if (temperature is not None) == 0:   
      print('Temperature is not set')
      return 0
   if (density is not None) == 0:   
      print('Density is not set')
      return 0
   if (n_ii_rc_data is not None) == 0:   
      print('N II recombination coefficients (n_ii_rc_data) are not set')
      return 0
   if (n_ii_rc_br is not None) == 0:   
      print('N II branching ratios (n_ii_rc_br) are not set')
      return 0
   if (wavelength is not None) == 0:   
      print('Wavelength is not given')
      return 0
   if ((temperature <= 0.e0) | (density <= 0.e0)):
      print('temperature = ', temperature, ', density = ', density)
      return 0
   
   wave = np.float64(0.0)
   rl_br = np.float64(0.0)
   g1 = np.float64(0.0)
   g2 = np.float64(0.0)
   temp4 = temperature / 10000.0
   loc1 = np.where(abs(n_ii_rc_br['wavelength'] - wavelength) <= 0.01)
   loc1=np.asarray(loc1[0])
   temp2 = len(loc1)
   if temp2 != 1:
      wavelength_min = np.amin(n_ii_rc_br['wavelength'][loc1])
      loc1 = np.where(n_ii_rc_br['wavelength'] == wavelength_min)
      loc1 = np.asarray(loc1[0])
   wave = np.float64(n_ii_rc_br['wavelength'][loc1])
   rl_br = np.float64(n_ii_rc_br['br'][loc1])

   loc1=loc1.item()
   #---------------------------------------
   if ((loc1 >= 0) & (loc1 <= 5)):
      # atomic transitions: 2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p E1 3P* - 3D  : 03  row
      #i = 2    #case A
      i = 3     #case B
      #---------------------------------------
   elif ((loc1 >= 6) & (loc1 <= 8)):
      # atomic transitions: 2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 3P* - 3S     : 04 row
      #i = 4    #Case: A
      i = 5     #Case: B
      #---------------------------------------
   elif ((loc1 >= 9) & (loc1 <= 14)):
      # atomic transitions: 2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 3P* - 3P     : 05 row
      #i = 6    #Case: A
      i = 7     #Case: B
      #---------------------------------------
   elif (loc1 == 15):
      # atomic transitions: 2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 1P* - 1P     : 08 row
      i = 8     #Case: A
      #i = 9    #Case: B
      #---------------------------------------
   elif (loc1 == 16):
      # atomic transitions: 2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 1P* - 1D     : 12 row
      i = 10     #Case: A
      #i = 11    #Case: B
      #---------------------------------------
   elif (loc1 == 17):
      # atomic transitions: 2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 1P* - 1S     : 13 row
      i = 12     #Case: A
      #i = 13    #Case: B
      #---------------------------------------
   elif (loc1 == 18):
      # atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 1P - 1D*     : 15 row
      i = 14     #Case: A
      #i = 15    #Case: B
      #---------------------------------------
   elif (loc1 == 19):
      # atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 1P - 1P*     : 17 row
      i = 16     #Case: A
      #i = 17    #Case: B
      #---------------------------------------
   elif ((loc1 >= 20) & (loc1 <= 25)):
      # atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3D - 3F*     : 19 row
      #i = 18    #Case: A
      i = 19     #Case: B
      #---------------------------------------
   elif ((loc1 >= 26) & (loc1 <= 32)):
      # atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3D - 3D*     : 20 row
      #i = 20    #Case: A
      i = 21     #Case: B
      #---------------------------------------
   elif ((loc1 >= 33) & (loc1 <= 38)):
      # atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3D - 3P*     : 21 row
      #i = 22    #Case: A
      i = 23     #Case: B
      #---------------------------------------
   elif ((loc1 >= 39) & (loc1 <= 44)):
      # atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).4s 3D - 3P*     : 22 row
      #i = 24    #Case: A
      i = 25     #Case: B
      #---------------------------------------
   elif ((loc1 >= 45) & (loc1 <= 47)):
      # atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3S - 3P*     : 24 row
      #i = 26    #Case: A
      i = 27     #Case: B
      #---------------------------------------
   elif ((loc1 >= 48) & (loc1 <= 50)):
      # atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).4s 3S - 3P*     : 26 row
      #i = 28    #Case: A
      i = 29     #Case: B
      #---------------------------------------
   elif ((loc1 >= 51) & (loc1 <= 56)):
      # atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3P - 3D*     : 28 row
      #i = 30    #Case: A
      i = 31     #Case: B
      #---------------------------------------
   elif ((loc1 >= 57) & (loc1 <= 62)):
      # atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3P - 3P*     : 29 row
      #i = 32    #Case: A
      i = 33     #Case: B
      #---------------------------------------
   elif ((loc1 >= 63) & (loc1 <= 68)):
      # atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).4s 3P - 3P*     : 30 row
      #i = 34    #Case: A
      i = 35     #Case: B
      #---------------------------------------
   elif (loc1 == 69):
      # atomic transitions: 2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 1D - 1F*     : 31 row
      i = 36     #Case: A
      #i = 37    #Case: B
      #---------------------------------------
   elif ((loc1 >= 70) & (loc1 <= 75)):
      # atomic transitions: 2s2.2p.(2P*).3d - 2s2.2p.(2P*).4p 3F* - 3D     : 36 row
      #i = 38     #Case: A
      i = 39     #Case: B
      #---------------------------------------
   elif ((loc1 >= 76) & (loc1 <= 81)):
      # atomic transitions: 2s2.2p.(2P*).3d - 2s2.2p.(2P*<3/2>).4f 3F* - 3G     : 39 row
      #i = 40    #Case: A
      i = 41     #Case: B
      #---------------------------------------
   elif (loc1 == 82):
      # atomic transitions: 2s2.2p.(2P*).3d - 2s2.2p.(2P*<3/2>).4f 1F* - 1G     : 58 row
      i = 42     #Case: A
      #i = 43    #Case: B
      #---------------------------------------
   elif ((loc1 >= 83) & (loc1 <= 88)):
      # atomic transitions: 3d 3D* - 4f 3F 4242     : 48 row
      #i = 44    #Case: A
      i = 45     #Case: B
      #---------------------------------------
   elif ((loc1 >= 89) & (loc1 <= 94)):
      # atomic transitions: 3d 3P* - 4f 3D 4435     : 55 row
      #i = 46   #Case: A
      i = 47     #Case: B
      #---------------------------------------
   elif (loc1 == 95):
      # atomic transitions: 3d 1D* - 4f 1F 4176     : 43 (RMT 42) row
      i = 48     #Case: A
      #i = 49     #Case: B
      #---------------------------------------
   elif (loc1 == 96):
      # atomic transitions: 3d 1P* - 4f 1D 4677     : 61 (rmt 62) row
      i = 49     #Case: A
      #i = 51    #Case: B
      #---------------------------------------
   elif (loc1 == 97):
      # Adopted from MOCASSIN, Ercolano et al. 2005MNRAS.362.1038E
      # atomic transitions: 3d 3F* - 4f 1G 4026     : 39b row
      # Case: A
      a0 = 0.108
      b0 = -0.754
      c0 = 2.587
      d0 = 0.719
      z0 = 2.
      br0 = 0.350
      aeff = 1.e-13 * z0 * a0 * (temp4 / z0 ** 2) ** (b0)
      aeff = aeff / (1. + c0 * (temp4 / z0 ** 2) ** (d0)) * br0
      #---------------------------------------
   elif (loc1 == 98):
      # Adopted from MOCASSIN, Ercolano et al. 2005MNRAS.362.1038E
      # atomic transitions: 3d 1F* - 4f 3G 4552     : 58a row
      # Case: A
      a0 = 0.326
      b0 = -0.754
      c0 = 2.587
      d0 = 0.719
      z0 = 2.
      br0 = 0.074
      aeff = 1.e-13 * z0 * a0 * (temp4 / z0 ** 2) ** (b0)
      aeff = aeff / (1. + c0 * (temp4 / z0 ** 2) ** (d0)) * br0
   else:   
      print('wavelength has an illegal value.')
   
   if ((loc1 >= 0) & (loc1 <= 96)):
      a0=np.float64(n_ii_rc_data['a'][i - 1])
      b0=np.float64(n_ii_rc_data['b'][i - 1])
      c0=np.float64(n_ii_rc_data['c'][i - 1])
      aeff = 10. ** (a0 + b0 * np.log10(temp4) + c0 * np.log10(temp4) ** 2)
   
   #niiRLs_Int = 100.0 * emissivity / hbeta_ems * abund
   #abund=line_flux/niiRLs_Int
   emissivity = (np.float64(aeff * rl_br) / np.float64(wave)) * np.float64(h_planck * c_speed * 1.e8)
   
   return emissivity

def calc_emiss_n_iii_rl(temperature=None, density=None,
                        wavelength=None, n_iii_rc_data=None):
   """
         This function calculates the emissivity
         for the given wavelength of N III recombination line
         by using the recombination coefficients from
         Pequignot et al. 1991A&A...251..680P.

     :Returns:
        type=double. This function returns the line emissivity.

     :Keywords:
         temperature   :     in, required, type=float
                             electron temperature
         density       :     in, required, type=float
                             electron density
         wavelength    :     in, required, type=float
                             Line Wavelength in Angstrom
         n_iii_rc_data  :     in, required, type=array/object
                             N III recombination coefficients

    :Examples:
       For example::

         >>> import pyequib
         >>> import atomneb
         >>> import os
         >>> base_dir = '../externals/atomneb/'
         >>> data_rc_dir = os.path.join('atomic-data-rc')
         >>> atom_rc_ppb91_file=os.path.join(base_dir,data_rc_dir, 'rc_PPB91.fits')
         >>> atom_rc_sh95_file= os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
         >>>
         >>> atom='n'
         >>> ion='iv' # N III
         >>> n_iii_rc_data=atomneb.read_aeff_ppb91(atom_rc_ppb91_file, atom, ion)
         >>> temperature=np.float64(10000.0)
         >>> density=np.float64(5000.0)
         >>> wavelength=4640.64
         >>> emiss_n_iii=pyequib.calc_abund_n_iii_rl(temperature=temperature, density=density,
         >>>                                 wavelength=wavelength,
         >>>                                 n_iii_rc_data=n_iii_rc_data)
         >>> print('Emissivity:', emiss_n_iii)
            Emissivity:   4.7908644e-24

     :Categories:
       Abundance Analysis, Recombination Lines, Emissivity

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.3.0

     :History:
         Based on  effective radiative recombination coefficients for N III lines from
         Pequignot, Petitjean, Boisson, C. 1991A&A...251..680P.

         10/05/2013, A. Danehkar, IDL code written.

         20/04/2017, A. Danehkar, Integration with AtomNeb.

         10/07/2019, A. Danehkar, Made a new function calc_emiss_n_iii_rl()
                          for calculating line emissivities and separated it
                          from calc_abund_n_iii_rl().

         03/10/2020, A. Danehkar, Transferred from IDL to Python.
   """
   # niiiRLstructure ={Wave:np.float64(0.0), Int:np.float64(0.0), Obs:np.float64(0.0), Abundance:np.float64(0.0)}
   h_planck = 6.62606957e-27 # erg s
   c_speed = 2.99792458e10 # cm/s 
   
   if (temperature is not None) == 0:   
      print('Temperature is not set')
      return 0
   if (density is not None) == 0:   
      print('Density is not set')
      return 0
   if (n_iii_rc_data is not None) == 0:   
      print('N III recombination coefficients (n_iii_rc_data) are not set')
      return 0
   if (wavelength is not None) == 0:   
      print('Wavelength is not given')
      return 0
   if ((temperature <= 0.e0) | (density <= 0.e0)):
      print('temperature = ', temperature, ', density = ', density)
      return 0
   
   lamb = np.float64(0.0)
   a = np.float64(0.0)
   b = np.float64(0.0)
   c = np.float64(0.0)
   d = np.float64(0.0)
   br = np.float64(0.0)
   aeff = np.float64(0.0)
   ion = ''
   
   z = 3.0 # ion level c^3+
   # equation (1) in 1991A&A...251..680P
   temp4 = 1.0e-4 * temperature / z ** 2
   loc1 = np.where(abs(n_iii_rc_data['wavelength'] - wavelength) <= 0.01)
   loc1=np.asarray(loc1[0])
   temp2 = len(loc1)
   if temp2 != 1:
      wavelength_min = np.amin(n_iii_rc_data['wavelength'][loc1])
      loc1 = np.where(n_iii_rc_data['wavelength'] == wavelength_min)
      loc1 = np.asarray(loc1[0])
      temp2 = len(loc1)
      if temp2 != 1:
         loc1 = np.where((n_iii_rc_data['wavelength'] == wavelength_min)&(n_iii_rc_data['case1'] == 'B'))
         loc1 = np.asarray(loc1[0])
   lamb = np.float64(n_iii_rc_data['wavelength'][loc1])
   a = np.float64(n_iii_rc_data['a'][loc1])
   b = np.float64(n_iii_rc_data['b'][loc1])
   c = np.float64(n_iii_rc_data['c'][loc1])
   d = np.float64(n_iii_rc_data['d'][loc1])
   br = np.float64(n_iii_rc_data['br'][loc1])
   # equation (1) in 1991A&A...251..680P
   aeff = 1.0e-13 * z * br
   aeff = aeff * (a * (temp4 ** b)) / (1. + c * (temp4 ** d))
   #niiiRLs_Int = 100.0 * (aeff/hbeta_aeff) * (4861.33/lamb) * abund
   #abund=line_flux/niiiRLs_Int
   emissivity = (np.float64(aeff) / np.float64(lamb)) * np.float64(h_planck * c_speed * 1.e8)
   
   return emissivity

def calc_emiss_o_ii_rl(temperature=None, density=None,
                       wavelength=None,
                       o_ii_rc_br=None, o_ii_rc_data=None):
   """
        This function calculates the emissivity
        for the given wavelength of O II recombination line
        by using the recombination coefficients from
        Storey 1994A&A...282..999S and Liu et al. 1995MNRAS.272..369L.

    :Returns:
       type=double. This function returns the line emissivity.

    :Keywords:
        temperature   :     in, required, type=float
                            electron temperature
        density       :     in, required, type=float
                            electron density
        wavelength    :     in, required, type=float
                            Line Wavelength in Angstrom
        o_ii_rc_br    :     in, required, type=array/object
                            O II branching ratios (Br)
        o_ii_rc_data  :     in, required, type=array/object
                            O II recombination coefficients

    :Examples:
       For example::

        >>> import pyequib
        >>> import atomneb
        >>> import os
        >>> base_dir = '../externals/atomneb/'
        >>> data_rc_dir = os.path.join('atomic-data-rc')
        >>> atom_rc_all_file= os.path.join(base_dir,data_rc_dir, 'rc_collection.fits')
        >>> atom_rc_sh95_file= os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
        >>>
        >>> atom='o'
        >>> ion='iii' # O II
        >>> o_ii_rc_data=atomneb.read_aeff_collection(atom_rc_all_file, atom, ion)
        >>> o_ii_rc_data_br=atomneb.read_aeff_collection(atom_rc_all_file, atom, ion, br=True)
        >>> temperature=np.float64(10000.0)
        >>> density=np.float64(5000.0)
        >>> wavelength=4613.68
        >>> emiss_o_ii=pyequib.calc_emiss_o_ii_rl(temperature=temperature, density=density,
        >>>                               wavelength=wavelength,
        >>>                               o_ii_rc_br=o_ii_rc_data_br, o_ii_rc_data=o_ii_rc_data,
        >>>                               h_i_aeff_data=h_i_aeff_data)
        >>> print('Emissivity:', emiss_o_ii)
           Emissivity:   5.9047319e-27

    :Categories:
      Abundance Analysis, Recombination Lines, Emissivity

    :Dirs:
     ./
         Main routines

    :Author:
      Ashkbiz Danehkar

    :Copyright:
      This library is released under a GNU General Public License.

    :Version:
      0.3.0

    :History:
        Based on recombination coefficients for O II lines from
        Storey 1994A&A...282..999S and Liu et al. 1995MNRAS.272..369L.

        Adopted from MIDAS script Roii.prg written by X.W.Liu.

        Revised based on scripts by Yong Zhang added to MOCASSIN, 02/2003
                          Ercolano et al. 2005MNRAS.362.1038E.

        10/05/2013, A. Danehkar, Translated to IDL code.

        25/04/2017, A. Danehkar, Integration with AtomNeb.

        10/07/2019, A. Danehkar, Made a new function calc_emiss_o_ii_rl()
                         for calculating line emissivities and separated it
                         from calc_abund_o_ii_rl().

        03/10/2020, A. Danehkar, Transferred from IDL to Python.
   """

#  oiiRLstructure ={Wave:np.float64(0.0), #REAL*8
#              Int:np.float64(0.0),
#              Obs:np.float64(0.0),
#              abundance:np.float64(0.0),
#              g1:long(0), #INTEGER
#              g2:long(0), #INTEGER
#              Mult1:'', #CHARACTER*7
#              Term1:'', #CHARACTER*9
#              Term2:'' #CHARACTER*9
#              }
   h_planck = np.float64(6.62606957e-27) # erg s
   c_speed = np.float64(2.99792458e10) # cm/s
   
   if (temperature is not None) == 0:   
      print('Temperature is not set')
      return 0
   if (density is not None) == 0:   
      print('Density is not set')
      return 0
   if (o_ii_rc_data is not None) == 0:   
      print('O II recombination coefficients (o_ii_rc_data) are not set')
      return 0
   if (o_ii_rc_br is not None) == 0:   
      print('O II branching ratios (o_ii_rc_br) are not set')
      return 0
   if (wavelength is not None) == 0:   
      print('Wavelength is not given')
      return 0
   if ((temperature <= 0.e0) | (density <= 0.e0)):
      print('temperature = ', temperature, ', density = ', density)
      return 0
   
   wave = np.float64(0.0)
   br_a = np.float64(0.0)
   br_b = np.float64(0.0)
   br_c = np.float64(0.0)
   g1 = np.float64(0.0)
   g2 = np.float64(0.0)
   temp4 = temperature / 10000.0
   loc1 = np.where(abs(o_ii_rc_br['wavelength'] - wavelength) <= 0.01)
   loc1=np.asarray(loc1[0])
   temp2 = len(loc1)
   if temp2 != 1:
      wavelength_min = np.amin(o_ii_rc_br['wavelength'][loc1])
      loc1 = np.where(o_ii_rc_br['wavelength'] == wavelength_min)
      loc1 = np.asarray(loc1[0])
   temp2 = len(loc1)
   if temp2 != 1:
      loc1 = np.amin(loc1)
   
   wave = np.float64(o_ii_rc_br['wavelength'][loc1])
   br_a = np.float64(o_ii_rc_br['br_a'][loc1])
   br_b = np.float64(o_ii_rc_br['br_b'][loc1])
   br_c = np.float64(o_ii_rc_br['br_c'][loc1])
   g1 = np.float64(o_ii_rc_br['g1'][loc1])
   g2 = np.float64(o_ii_rc_br['g2'][loc1])
   densi = np.float64(density)
   log10ne = np.log10(densi)

   #---------------------------------------
   if (((loc1 >= 0) & (loc1 <= 182))):
      # 4f-3d transitions
      a = np.float64(o_ii_rc_data['a4'][0]) # 0.232
      b = np.float64(o_ii_rc_data['b'][0]) #-0.92009
      c = np.float64(o_ii_rc_data['c'][0]) # 0.15526
      d = np.float64(o_ii_rc_data['d'][0]) # 0.03442
      # an = [0.236, 0.232, 0.228, 0.222]
      an = np.array([np.float64(o_ii_rc_data['a2'][0]),
                     np.float64(o_ii_rc_data['a4'][0]),
                     np.float64(o_ii_rc_data['a5'][0]),
                     np.float64(o_ii_rc_data['a6'][0])])
      if (log10ne <= 2):   
         a = an[0]
      else:   
         if ((log10ne > 2) & (log10ne <= 4)):
            a = an[0] + (an[1] - an[0]) / 2. * (log10ne - 2.)
         else:   
            if ((log10ne > 4) & (log10ne <= 5)):
               a = an[1] + (an[2] - an[1]) * (log10ne - 2.)
            else:   
               if ((log10ne > 5) & (log10ne <= 6)):
                  a = an[2] + (an[3] - an[2]) * (log10ne - 2.)
               else:   
                  a = an[3]
      
      temp4 = temperature / 10000.0
      aeff = 1.0e-14 * a * temp4 ** (b)
      aeff = aeff * (1. + c * (1. - temp4) + d * (1. - temp4) ** 2)
      
      # data for 1000 < T < 5000 K, Ne = 100/cm3
      if (temp4 <= 0.5):   
         a = 0.236
         b = -1.07552
         c = -0.04843
         aeff = 1.0e-14 * a * temp4 ** (b + c * np.log(temp4))
      br = br_b
      #---------------------------------------
   elif ((loc1 >= 183) & (loc1 <= 218)):
      # 3d-3p ^4F transitions (Case A=B=C for a,b,c,d# Br diff. slightly, adopt Case B)
      a = np.float64(o_ii_rc_data['a4'][1]) # 0.876
      b = np.float64(o_ii_rc_data['b'][1]) # -0.73465
      c = np.float64(o_ii_rc_data['c'][1]) # 0.13689
      d = np.float64(o_ii_rc_data['d'][1]) # 0.06220
      # an = [0.876, 0.876, 0.877, 0.880] #a for logNe = 2,4,5,6 (LSBC95, Tab.5a)
      an = np.array([np.float64(o_ii_rc_data['a2'][1]),
                     np.float64(o_ii_rc_data['a4'][1]),
                     np.float64(o_ii_rc_data['a5'][1]),
                     np.float64(o_ii_rc_data['a6'][1])])
      if (log10ne <= 2):   
         a = an[0]
      else:   
         if ((log10ne > 2) & (log10ne <= 4)):
            a = an[0] + (an[1] - an[0]) / 2. * (log10ne - 2.)
         else:   
            if ((log10ne > 4) & (log10ne <= 5)):
               a = an[1] + (an[2] - an[1]) * (log10ne - 2.)
            else:   
               if ((log10ne > 5) & (log10ne <= 6)):
                  a = an[2] + (an[3] - an[2]) * (log10ne - 2.)
               else:   
                  a = an[3]
      
      aeff = 1.e-14 * a * temp4 ** (b)
      aeff = aeff * (1. + c * (1. - temp4) + d * (1. - temp4) ** 2)
      # for 1000 < T < 5000 K, Ne = 100/cm3
      if (temp4 <= 0.5):   
         a = 0.878
         b = -0.86175
         c = -0.02470
         aeff = 1.e-14 * a * temp4 ** (b + c * np.log(temp4))
      br = br_b
      #---------------------------------------
   elif ((loc1 >= 219) & (loc1 <= 309)):
      # 3d-3p ^4D, ^4P transitions
      a = np.float64(o_ii_rc_data['a4'][3]) # 0.745
      b = np.float64(o_ii_rc_data['b'][3]) # -0.74621
      c = np.float64(o_ii_rc_data['c'][3]) # 0.15710
      d = np.float64(o_ii_rc_data['d'][3]) # 0.07059
      #an = [0.727,0.726,0.725,0.726] # Case: A
      # an = [0.747, 0.745, 0.744, 0.745] # Case: B
      #an = [0.769,0.767,0.766,0.766] # Case: C
      an = np.array([np.float64(o_ii_rc_data['a2'][3]),
                     np.float64(o_ii_rc_data['a4'][3]),
                     np.float64(o_ii_rc_data['a5'][3]),
                     np.float64(o_ii_rc_data['a6'][3])])
      if (log10ne <= 2):   
         a = an[0]
      else:   
         if ((log10ne > 2) & (log10ne <= 4)):
            a = an[0] + (an[1] - an[0]) / 2. * (log10ne - 2.)
         else:   
            if ((log10ne > 4) & (log10ne <= 5)):
               a = an[1] + (an[2] - an[1]) * (log10ne - 2.)
            else:   
               if ((log10ne > 5) & (log10ne <= 6)):
                  a = an[2] + (an[3] - an[2]) * (log10ne - 2.)
               else:   
                  a = an[3]
      aeff = 1.e-14 * a * temp4 ** (b)
      aeff = aeff * (1. + c * (1. - temp4) + d * (1. - temp4) ** 2)
      # for 1000 < T < 5000 K, Ne = 100/cm3
      if (temp4 <= 0.5):   
         a = 0.747
         b = -0.89382
         c = -0.02906
         aeff = 1.e-14 * a * temp4 ** (b + c * np.log(temp4))
      br = br_b
      #---------------------------------------
   elif ((loc1 >= 310) & (loc1 <= 327)):
      # 3d-3p ^2F transitions
      a = np.float64(o_ii_rc_data['a4'][5]) # 0.745
      b = np.float64(o_ii_rc_data['b'][5]) # -0.74621
      c = np.float64(o_ii_rc_data['c'][5]) # 0.15710
      d = np.float64(o_ii_rc_data['d'][5]) # 0.07059
      #an = [0.727, 0.726, 0.725, 0.726] # Case: A
      #an = [0.747,0.745,0.744,0.745] # Case: B
      #an = [0.769,0.767,0.766,0.766] # Case: C
      an = np.array([o_ii_rc_data['a2'][5], o_ii_rc_data['a4'][5], o_ii_rc_data['a5'][5], o_ii_rc_data['a6'][5]])
      if (log10ne <= 2):   
         a = an[0]
      else:   
         if ((log10ne > 2) & (log10ne <= 4)):
            a = an[0] + (an[1] - an[0]) / 2. * (log10ne - 2.)
         else:   
            if ((log10ne > 4) & (log10ne <= 5)):
               a = an[1] + (an[2] - an[1]) * (log10ne - 2.)
            else:   
               if ((log10ne > 5) & (log10ne <= 6)):
                  a = an[2] + (an[3] - an[2]) * (log10ne - 2.)
               else:   
                  a = an[3]
      aeff = 1.e-14 * a * temp4 ** (b)
      aeff = aeff * (1. + c * (1. - temp4) + d * (1. - temp4) ** 2)
      # for 1000 < T < 5000 K, Ne = 100/cm3
      if (temp4 <= 0.5):   
         a = 0.747
         b = -0.89382
         c = -0.02906
         aeff = 1.e-14 * a * temp4 ** (b + c * np.log(temp4))
      br = br_a
      #---------------------------------------
   elif ((loc1 >= 328) & (loc1 <= 357)):
      # 3d-3p ^2D transitions
      a = np.float64(o_ii_rc_data['a4'][11]) # 0.601
      b = np.float64(o_ii_rc_data['b'][11]) # -0.79533
      c = np.float64(o_ii_rc_data['c'][11]) # 0.15314
      d = np.float64(o_ii_rc_data['d'][11]) # 0.05322
      #an = [0.603, 0.601, 0.600, 0.599] # Case: A
      #an = [0.620,0.618,0.616,0.615] # Case: C
      an = np.array([np.float64(o_ii_rc_data['a2'][11]),
                     np.float64(o_ii_rc_data['a4'][11]),
                     np.float64(o_ii_rc_data['a5'][11]),
                     np.float64(o_ii_rc_data['a6'][11])])
      if (log10ne <= 2):
         a = an[0]
      else:   
         if ((log10ne > 2) & (log10ne <= 4)):
            a = an[0] + (an[1] - an[0]) / 2. * (log10ne - 2.)
         else:   
            if ((log10ne > 4) & (log10ne <= 5)):
               a = an[1] + (an[2] - an[1]) * (log10ne - 2.)
            else:   
               if ((log10ne > 5) & (log10ne <= 6)):
                  a = an[2] + (an[3] - an[2]) * (log10ne - 2.)
               else:   
                  a = an[3]
      aeff = 1.e-14 * a * temp4 ** (b)
      aeff = aeff * (1. + c * (1. - temp4) + d * (1. - temp4) ** 2)
      # for 1000 < T < 5000 K, Ne = 100/cm3
      if (temp4 <= 0.5):   
         a = 0.603
         b = -0.94025
         c = -0.03467
         aeff = 1.e-14 * a * temp4 ** (b + c * np.log(temp4))
      br = br_a
      #---------------------------------------
   elif ((loc1 >= 358) & (loc1 <= 387)):
      # 3d-3p ^2P transitions
      a = np.float64(o_ii_rc_data['a4'][13]) # 0.524
      b = np.float64(o_ii_rc_data['b'][13]) # -0.78448
      c = np.float64(o_ii_rc_data['c'][13]) # 0.13681
      d = np.float64(o_ii_rc_data['d'][13]) # 0.05608
      #an = [0.526, 0.524, 0.523, 0.524] # Case: A
      #an = [0.538,0.536,0.535,0.536] # Case: C
      an = np.array([np.float64(o_ii_rc_data['a2'][13]),
                     np.float64(o_ii_rc_data['a4'][13]),
                     np.float64(o_ii_rc_data['a5'][13]),
                     np.float64(o_ii_rc_data['a6'][13])])
      if (log10ne <= 2):   
         a = an[0]
      else:   
         if ((log10ne > 2) & (log10ne <= 4)):
            a = an[0] + (an[1] - an[0]) / 2. * (log10ne - 2.)
         else:   
            if ((log10ne > 4) & (log10ne <= 5)):
               a = an[1] + (an[2] - an[1]) * (log10ne - 2.)
            else:   
               if ((log10ne > 5) & (log10ne <= 6)):
                  a = an[2] + (an[3] - an[2]) * (log10ne - 2.)
               else:   
                  a = an[3]
      aeff = 1.e-14 * a * temp4 ** (b)
      aeff = aeff * (1. + c * (1. - temp4) + d * (1. - temp4) ** 2)
      # for 1000 < T < 5000 K, Ne = 100/cm3
      if (temp4 <= 0.5):   
         a = 0.526
         b = -0.91758
         c = -0.03120
         aeff = 1.e-14 * a * temp4 ** (b + c * np.log(temp4))
      br = br_a
      #---------------------------------------
   elif ((loc1 >= 388) & (loc1 <= 395)):
      # 3p-3s ^4D - ^4P transitions
      #an = [34.7,34.9,35.1,35.0] #a  Case: A
      #a =  36.2
      #b =  -0.749
      #c =  0.023
      #d =  0.074
      an = np.array([36.0, 36.2, 36.4, 36.3]) #a  Case: B
      a = 36.2
      b = -0.736
      c = 0.033
      d = 0.077
      if (log10ne <= 2):   
         a = an[0]
      else:   
         if ((log10ne > 2) & (log10ne <= 4)):
            a = an[0] + (an[1] - an[0]) / 2. * (log10ne - 2.)
         else:   
            if ((log10ne > 4) & (log10ne <= 5)):
               a = an[1] + (an[2] - an[1]) * (log10ne - 2.)
            else:   
               if ((log10ne > 5) & (log10ne <= 6)):
                  a = an[2] + (an[3] - an[2]) * (log10ne - 2.)
               else:   
                  a = an[3]
      aeff = 1.e-14 * a * temp4 ** (b)
      aeff = aeff * (1. + c * (1. - temp4) + d * (1. - temp4) ** 2)
      # for 1000 < T < 5000 K, Ne = 100/cm3
      if (temp4 <= 0.5):   
         a = 36.288
         b = -0.75421
         c = 0.02883
         d = 0.01213
         aeff = 1.e-14 * a * temp4 ** (b + c * np.log(temp4) + d * np.log(temp4) ** 2)
      br = br_b
      #---------------------------------------
   elif ((loc1 >= 396) & (loc1 <= 402)):
      # 3p-3s ^4P - ^4P transitions
      #an = [10.4,10.4,10.5,10.4] # Case: A
      #a =  10.4
      #b =  -0.721
      #c =  0.073
      #d =  0.072
      an = np.array([14.6, 14.6, 14.7, 14.6]) # Case: B
      a = 14.6
      b = -0.732
      c = 0.081
      d = 0.066
      if (log10ne <= 2):   
         a = an[0]
      else:   
         if ((log10ne > 2) & (log10ne <= 4)):
            a = an[0] + (an[1] - an[0]) / 2. * (log10ne - 2.)
         else:   
            if ((log10ne > 4) & (log10ne <= 5)):
               a = an[1] + (an[2] - an[1]) * (log10ne - 2.)
            else:   
               if ((log10ne > 5) & (log10ne <= 6)):
                  a = an[2] + (an[3] - an[2]) * (log10ne - 2.)
               else:   
                  a = an[3]
      aeff = 1.e-14 * a * temp4 ** (b)
      aeff = aeff * (1. + c * (1. - temp4) + d * (1. - temp4) ** 2)
      # for 1000 < T < 5000 K, Ne = 100/cm3
      if (temp4 <= 0.5):   
         a = 14.656
         b = -0.80449
         c = 0.00018
         d = 0.00517
         aeff = 1.e-14 * a * temp4 ** (b + c * np.log(temp4) + d * np.log(temp4) ** 2)
      br = br_b
      #---------------------------------------
   elif ((loc1 >= 403) & (loc1 <= 405)):
      # 3p-3s ^4S - ^4P transitions
      #an = [0.90,0.90,0.90,1.00] # Case: A
      #a =  0.90
      #b =  -0.485
      #c =  -0.047
      #d =  0.140
      an = np.array([4.80, 4.90, 4.90, 4.90]) # Case: B
      a = 4.90
      b = -0.730
      c = -0.003
      d = 0.057
      if (log10ne <= 2):   
         a = an[0]
      else:   
         if ((log10ne > 2) & (log10ne <= 4)):
            a = an[0] + (an[1] - an[0]) / 2. * (log10ne - 2.)
         else:   
            if ((log10ne > 4) & (log10ne <= 5)):
               a = an[1] + (an[2] - an[1]) * (log10ne - 2.)
            else:   
               if ((log10ne > 5) & (log10ne <= 6)):
                  a = an[2] + (an[3] - an[2]) * (log10ne - 2.)
               else:   
                  a = an[3]
      aeff = 1.e-14 * a * temp4 ** (b)
      aeff = aeff * (1. + c * (1. - temp4) + d * (1. - temp4) ** 2)
      # for 1000 < T < 5000 K, Ne = 100/cm3
      if (temp4 <= 0.5):   
         a = 4.8340
         b = -0.71947
         c = 0.02544
         d = 0.00936
         aeff = 1.e-14 * a * temp4 ** (b + c * np.log(temp4) + d * np.log(temp4) ** 2)
      br = br_b
      #---------------------------------------
   elif ((loc1 >= 406) & (loc1 <= 408)):
      # 3p-3s ^2D - ^2P transitions
      an = np.array([2.40, 2.40, 2.50, 2.60]) # Case: A
      a = 2.40
      b = -0.550
      c = -0.051
      d = 0.178
      #an = [14.5,14.6,14.5,14.3] # Case: C
      #a =  14.6
      #b =  -0.736
      #c =  0.068
      #d =  0.066
      if (log10ne <= 2):   
         a = an[0]
      else:   
         if ((log10ne > 2) & (log10ne <= 4)):
            a = an[0] + (an[1] - an[0]) / 2. * (log10ne - 2.)
         else:   
            if ((log10ne > 4) & (log10ne <= 5)):
               a = an[1] + (an[2] - an[1]) * (log10ne - 2.)
            else:   
               if ((log10ne > 5) & (log10ne <= 6)):
                  a = an[2] + (an[3] - an[2]) * (log10ne - 2.)
               else:   
                  a = an[3]
      aeff = 1.e-14 * a * temp4 ** (b)
      aeff = aeff * (1. + c * (1. - temp4) + d * (1. - temp4) ** 2)
      # for 1000 < T < 5000 K, Ne = 100/cm3
      if (temp4 <= 0.5):   
         a = 2.3616
         b = -0.46263
         c = 0.14697
         d = 0.03856
         aeff = 1.e-14 * a * temp4 ** (b + c * np.log(temp4) + d * np.log(temp4) ** 2)
      br = br_a
      #---------------------------------------
   elif ((loc1 >= 409) & (loc1 <= 412)):
      # 3p-3s ^2P - ^2P transitions
      an = np.array([1.10, 1.20, 1.20, 1.20]) # Case: A
      a = 1.20
      b = -0.523
      c = -0.044
      d = 0.173
      #an = [1.30,1.40,1.40,1.40] # Case: C
      #a =  1.40
      #b =  -0.565
      #c =  -0.042
      #d =  0.158
      if (log10ne <= 2):   
         a = an[0]
      else:   
         if ((log10ne > 2) & (log10ne <= 4)):
            a = an[0] + (an[1] - an[0]) / 2. * (log10ne - 2.)
         else:   
            if ((log10ne > 4) & (log10ne <= 5)):
               a = an[1] + (an[2] - an[1]) * (log10ne - 2.)
            else:   
               if ((log10ne > 5) & (log10ne <= 6)):
                  a = an[2] + (an[3] - an[2]) * (log10ne - 2.)
               else:   
                  a = an[3]
      aeff = 1.e-14 * a * temp4 ** (b)
      aeff = aeff * (1. + c * (1. - temp4) + d * (1. - temp4) ** 2)
      # for 1000 < T < 5000 K, Ne = 100/cm3
      if (temp4 <= 0.5):   
         a = 1.1198
         b = -0.44147
         c = 0.13837
         d = 0.03191
         aeff = 1.e-14 * a * temp4 ** (b + c * np.log(temp4) + d * np.log(temp4) ** 2)
      br = br_a
      #---------------------------------------
   elif ((loc1 >= 413) & (loc1 <= 414)):
      # 3p-3s ^2S - ^2P transitions
      an = np.array([0.40, 0.40, 0.40, 0.40]) # Case: A
      a = 0.40
      b = -0.461
      c = -0.083
      d = 0.287
      #an = [0.50,0.50,0.50,0.60] # Case: C
      #a =  0.50
      #b =  -0.547
      #c =  -0.074
      #d =  0.244
      if (log10ne <= 2):   
         a = an[0]
      else:   
         if ((log10ne > 2) & (log10ne <= 4)):
            a = an[0] + (an[1] - an[0]) / 2. * (log10ne - 2.)
         else:   
            if ((log10ne > 4) & (log10ne <= 5)):
               a = an[1] + (an[2] - an[1]) * (log10ne - 2.)
            else:   
               if ((log10ne > 5) & (log10ne <= 6)):
                  a = an[2] + (an[3] - an[2]) * (log10ne - 2.)
               else:   
                  a = an[3]
      aeff = 1.e-14 * a * temp4 ** (b)
      aeff = aeff * (1. + c * (1. - temp4) + d * (1. - temp4) ** 2)
      # for 1000 < T < 5000 K, Ne = 100/cm3
      if (temp4 <= 0.5):   
         a = 0.3922
         b = -0.35043
         c = 0.26366
         d = 0.06666
         aeff = 1.e-14 * a * temp4 ** (b + c * np.log(temp4) + d * np.log(temp4) ** 2)
      br = br_a
   else:   
      print('wavelength has an illegal value.')
   
   # Ems1 = aeff * (h_Planck*c_Speed*1.e8) /Wave * g2 * Br_A#[loc1]
   # oiiRLs_Int = 100. * Ems1 / hbeta_ems * abund
   #abund=line_flux/oiiRLs_Int
   
   emissivity = (np.float64(aeff * g2 * br) / np.float64(wave)) * np.float64(h_planck * c_speed * 1.e8)
   
   return emissivity

def calc_emiss_ne_ii_rl(temperature=None, density=None,
                        wavelength=None, ne_ii_rc_data=None):
   """
         This function calculates the emissivity
         for the given wavelength of Ne II recombination line
         by using the recombination coefficients from
         Kisielius et al. (1998) & Storey (unpublished).

     :Returns:
        type=double. This function returns the line emissivity.

     :Keywords:
         temperature   :     in, required, type=float
                             electron temperature
         density       :     in, required, type=float
                             electron density
         wavelength    :     in, required, type=float
                             Line Wavelength in Angstrom
         ne_ii_rc_data  :    in, required, type=array/object
                             Ne II recombination coefficients

    :Examples:
       For example::

         >>> import pyequib
         >>> import atomneb
         >>> import os
         >>> base_dir = '../externals/atomneb/'
         >>> data_rc_dir = os.path.join('atomic-data-rc')
         >>> atom_rc_all_file= os.path.join(base_dir,data_rc_dir, 'rc_collection.fits')
         >>> atom_rc_sh95_file= os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
         >>>
         >>> atom='ne'
         >>> ion='iii' # Ne II
         >>> ne_ii_rc_data=atomneb.read_aeff_collection(atom_rc_all_file, atom, ion)
         >>> temperature=np.float64(10000.0)
         >>> density=np.float64(5000.0)
         >>> wavelength=3777.14
         >>> emiss_ne_ii=pyequib.calc_emiss_ne_ii_rl(temperature=temperature, density=density,
         >>>                                 wavelength=wavelength,
         >>>                                 ne_ii_rc_data=ne_ii_rc_data, h_i_aeff_data=h_i_aeff_data)
         >>> print('Emissivity:', emiss_ne_ii)
            Emissivity:   1.5996881e-25

     :Categories:
       Abundance Analysis, Recombination Lines, Emissivity

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.3.0

     :History:
         Based on effective radiative recombination coefficients for Ne II lines
         from Kisielius et al. 1998A&AS..133..257K & Storey (unpublished).

         Adopted from MOCASSIN, Ercolano et al. 2005MNRAS.362.1038E.

         02/2003, Yong Zhang, scripts added to MOCASSIN.

         14/05/2013, A. Danehkar, Translated to IDL code.

         10/04/2017, A. Danehkar, Integration with AtomNeb.

         10/07/2019, A. Danehkar, Made a new function calc_emiss_ne_ii_rl()
                          for calculating line emissivities and separated it
                          from calc_abund_ne_ii_rl().

         03/10/2020, A. Danehkar, Transferred from IDL to Python.
   """
   # neiiRLstructure ={Wave:np.float64(0.0), Int:np.float64(0.0), Obs:np.float64(0.0), Abundance:np.float64(0.0)}
   h_planck = 6.62606957e-27 # erg s
   c_speed = 2.99792458e10 # cm/s 
   
   if (temperature is not None) == 0:   
      print('Temperature is not set')
      return 0
   if (density is not None) == 0:   
      print('Density is not set')
      return 0
   if (ne_ii_rc_data is not None) == 0:   
      print('Ne II recombination coefficients (ne_ii_rc_data) are not set')
      return 0
   if (wavelength is not None) == 0:   
      print('Wavelength is not given')
      return 0
   if ((temperature <= 0.e0) | (density <= 0.e0)):
      print('temperature = ', temperature, ', density = ', density)
      return 0
   
   lamb = np.float64(0.0)
   a = np.float64(0.0)
   b = np.float64(0.0)
   c = np.float64(0.0)
   d = np.float64(0.0)
   f = np.float64(0.0)
   br = np.float64(0.0)
   aeff = np.float64(0.0)
   ion = ''
   
   z = 3.0 # ion level c^3+
   # equation (1) in 1991A&A...251..680P
   temp4 = temperature / 10000.0
   loc1 = np.where(abs(ne_ii_rc_data['wavelength'] - wavelength) <= 0.01)
   loc1=np.asarray(loc1[0])
   temp2 = len(loc1)
   if temp2 != 1:
      wavelength_min = np.amin(ne_ii_rc_data['wavelength'][loc1])
      loc1 = np.where(ne_ii_rc_data['wavelength'] == wavelength_min)
      loc1 = np.asarray(loc1[0])
   lamb = np.float64(ne_ii_rc_data['wavelength'][loc1])
   a = np.float64(ne_ii_rc_data['a'][loc1])
   b = np.float64(ne_ii_rc_data['b'][loc1])
   c = np.float64(ne_ii_rc_data['c'][loc1])
   d = np.float64(ne_ii_rc_data['d'][loc1])
   f = np.float64(ne_ii_rc_data['f'][loc1])
   br = np.float64(ne_ii_rc_data['br'][loc1])
   # equation (1) in 1991A&A...251..680P
   aeff = 1.0e-14 * (a * (temp4 ** f)) * br
   aeff = aeff * (1. + b * (1.0 - temp4) + c * (1.0 - temp4) ** 2 + d * (1.0 - temp4) ** 3)
   #neiiRLs_Int = 100.0 * (aeff/hbeta_aeff) * (4861.33/lamb) * abund
   #abund=line_flux/neiiRLs_Int
   emissivity = (np.float64(aeff) / np.float64(lamb)) * np.float64(h_planck * c_speed * 1.e8)
   
   return emissivity

def calc_abund_he_i_rl(temperature=None, density=None,
                       linenum=None, line_flux=None,
                       he_i_aeff_data=None, h_i_aeff_data=None):
   """
         This function determines the ionic abundance from the observed
         flux intensity for the given wavelength of He I recombination line
         by using the recombination coefficients from Porter et al.
         2012MNRAS.425L..28P.

     :Returns:
        type=double. This function returns the ionic abundanc.

     :Keywords:
         temperature    :    in, required, type=float
                             electron temperature
         density        :    in, required, type=float
                             electron density
         linenum        :    in, required, type=int
                             Line Number for Wavelength

                             Wavelength=4120.84:linenum=7,

                             Wavelength=4387.93: linenum=8,

                             Wavelength=4437.55: linenum=9,

                             Wavelength=4471.50: linenum=10,

                             Wavelength=4921.93: linenum=12,

                             Wavelength=5015.68: linenum=13,

                             Wavelength=5047.74: linenum=14,

                             Wavelength=5875.66: linenum=15,

                             Wavelength=6678.16: linenum=16,

                             Wavelength=7065.25: linenum=17,

                             Wavelength=7281.35: linenum=18.

         line_flux      :    in, required, type=float
                             line flux intensity
         he_i_aeff_data :    in, required, type=array/object
                             He I recombination coefficients
         h_i_aeff_data  :    in, required, type=array/object
                             H I recombination coefficients

    :Examples:
       For example::

         >>> import pyequib
         >>> import atomneb
         >>> import os
         >>> base_dir = '../externals/atomneb/'
         >>> data_rc_dir = os.path.join('atomic-data-rc')
         >>> atom_rc_he_i_file= os.path.join(base_dir,data_rc_dir, 'rc_he_ii_PFSD12.fits')
         >>> atom_rc_sh95_file= os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
         >>> atom='h'
         >>> ion='ii' # H I
         >>> h_i_rc_data=atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)
         >>> h_i_aeff_data=h_i_rc_data['aeff'][0]
         >>> atom='he'
         >>> ion='ii' # He I
         >>> he_i_rc_data=atomneb.read_aeff_he_i_pfsd12(atom_rc_he_i_file, atom, ion)
         >>> he_i_aeff_data=he_i_rc_data['aeff'][0]
         >>> temperature=np.float64(10000.0)
         >>> density=np.float64(5000.0)
         >>> he_i_4471_flux= 2.104
         >>> linenum=10# 4471.50
         >>> abund_he_i=pyequib.calc_abund_he_i_rl(temperature=temperature, density=density,
         >>>                                  linenum=linenum, line_flux=he_i_4471_flux,
         >>>                                  he_i_aeff_data=he_i_aeff_data, h_i_aeff_data=h_i_aeff_data)
         >>> print('N(He^+)/N(H^+):', abund_he_i)
            N(He^+)/N(H^+):     0.040848393

     :Categories:
       Abundance Analysis, Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.3.0

     :History:
         Based on improved He I emissivities in the case B
         from Porter et al. 2012MNRAS.425L..28P

         15/12/2013, A. Danehkar, IDL code written.

         20/03/2017, A. Danehkar, Integration with AtomNeb.

         10/07/2019, A. Danehkar, Made a new function calc_emiss_he_i_rl()
                          for calculating line emissivities and separated it
                          from calc_abund_he_i_rl().

         03/10/2020, A. Danehkar, Transferred from IDL to Python.
   """
   if (temperature is not None) == 0:   
      print('Temperature is not set')
      return 0
   if (density is not None) == 0:   
      print('Density is not set')
      return 0
   if (he_i_aeff_data is not None) == 0:   
      print('He I recombination coefficients (he_i_aeff_data) are not set')
      return 0
   if (h_i_aeff_data is not None) == 0:   
      print('H I recombination coefficients (h_i_aeff_data) are not set')
      return 0
   if (linenum is not None) == 0:   
      print('Line Number for Wavelength is not given')
      return 0
   if (line_flux is not None) == 0:   
      print('Line flux intensity (line_flux) is not given')
      return 0
   if ((temperature <= 0.e0) | (density <= 0.e0)):
      print('temperature = ', temperature, ', density = ', density)
      return 0
   emissivity_hbeta = calc_emiss_h_beta(temperature=temperature, density=density, h_i_aeff_data=h_i_aeff_data)
   
   emissivity = calc_emiss_he_i_rl(temperature=temperature, density=density, linenum=linenum, he_i_aeff_data=he_i_aeff_data)
   abund = (emissivity_hbeta / emissivity) * np.float64(line_flux / 100.0)
   
   return abund

def calc_abund_he_ii_rl(temperature=None, density=None,
                        line_flux=None,
                        he_ii_aeff_data=None, h_i_aeff_data=None):
   """
         This function determines the ionic abundance from the observed
         flux intensity for the He II recombination line 4686 A
         by using the helium emissivities from
         Storey & Hummer, 1995MNRAS.272...41S.

     :Returns:
        type=double. This function returns the ionic abundanc.

     :Keywords:
         temperature     :   in, required, type=float
                             electron temperature
         density         :   in, required, type=float
                             electron density
         line_flux       :   in, required, type=float
                             line flux intensity
         he_ii_aeff_data :   in, required, type=array/object
                             He II recombination coefficients
         h_i_aeff_data   :   in, required, type=array/object
                             H I recombination coefficients

    :Examples:
       For example::

         >>> import pyequib
         >>> import atomneb
         >>> import os
         >>> base_dir = '../externals/atomneb/'
         >>> data_rc_dir = os.path.join('atomic-data-rc')
         >>> atom_rc_he_i_file= os.path.join(base_dir,data_rc_dir, 'rc_he_ii_PFSD12.fits')
         >>> atom_rc_sh95_file= os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
         >>> atom='h'
         >>> ion='ii' # H I
         >>> h_i_rc_data=atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)
         >>> h_i_aeff_data=h_i_rc_data['aeff'][0]
         >>> atom='he'
         >>> ion='iii' # He II
         >>> he_ii_rc_data=atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)
         >>> he_ii_aeff_data=he_ii_rc_data['aeff'][0]
         >>> temperature=np.float64(10000.0)
         >>> density=np.float64(5000.0)
         >>> he_ii_4686_flux = 135.833
         >>> abund_he_ii=pyequib.calc_abund_he_ii_rl(temperature=temperature, density=density,
         >>>                                 line_flux=he_ii_4686_flux,
         >>>                                 he_ii_aeff_data=he_ii_aeff_data, h_i_aeff_data=h_i_aeff_data)
         >>> print('N(He^2+)/N(H^+):', abund_he_ii)
            N(He^2+)/N(H^+):      0.11228817

     :Categories:
       Abundance Analysis, Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.3.0

     :History:
         Based on He II emissivities
         from Storey & Hummer, 1995MNRAS.272...41S.

         15/12/2013, A. Danehkar, IDL code written.

         02/04/2017, A. Danehkar, Integration with AtomNeb.

         10/07/2019, A. Danehkar, Made a new function calc_emiss_he_ii_rl()
                          for calculating line emissivities and separated it
                          from calc_abund_he_ii_rl().

         03/10/2020, A. Danehkar, Transferred from IDL to Python.
   """
   if (temperature is not None) == 0:   
      print('Temperature is not set')
      return 0
   if (density is not None) == 0:   
      print('Density is not set')
      return 0
   if (he_ii_aeff_data is not None) == 0:   
      print('He II recombination coefficients (he_ii_aeff_data) are not set')
      return 0
   if (h_i_aeff_data is not None) == 0:   
      print('H I recombination coefficients (h_i_aeff_data) are not set')
      return 0
   if (line_flux is not None) == 0:   
      print('Line flux intensity (line_flux) is not given')
      return 0
   if ((temperature <= 0.e0) | (density <= 0.e0)):
      print('temperature = ', temperature, ', density = ', density)
      return 0
   
   emissivity_hbeta = calc_emiss_h_beta(temperature=temperature, density=density, h_i_aeff_data=h_i_aeff_data)
   
   emissivity = calc_emiss_he_ii_rl(temperature=temperature, density=density, he_ii_aeff_data=he_ii_aeff_data)
   abund = (emissivity_hbeta / emissivity) * np.float64(line_flux / 100.0)
   
   return abund

def calc_abund_c_ii_rl(temperature=None, density=None,
                       wavelength=None, line_flux=None,
                       c_ii_rc_data=None, h_i_aeff_data=None):
   """
        This function determines the ionic abundance from the observed
        flux intensity for the given wavelength of C II recombination line
        by using the recombination coefficients from
        from Davey et al. (2000) 2000A&AS..142...85D.

    :Returns:
       type=double. This function returns the ionic abundanc.

    :Keywords:
        temperature   :     in, required, type=float
                            electron temperature
        density       :     in, required, type=float
                            electron density
        wavelength    :     in, required, type=float
                            Line Wavelength in Angstrom
        line_flux     :     in, required, type=float
                            line flux intensity
        c_ii_rc_data  :     in, required, type=array/object
                            C II recombination coefficients
        h_i_aeff_data :     in, required, type=array/object
                            H I recombination coefficients

    :Examples:
       For example::

        >>> import pyequib
        >>> import atomneb
        >>> import os
        >>> base_dir = '../externals/atomneb/'
        >>> data_rc_dir = os.path.join('atomic-data-rc')
        >>> atom_rc_all_file= os.path.join(base_dir,data_rc_dir, 'rc_collection.fits')
        >>> atom_rc_sh95_file= os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
        >>> atom='h'
        >>> ion='ii' # H I
        >>> h_i_rc_data=atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)
        >>> h_i_aeff_data=h_i_rc_data['aeff'][0]
        >>> atom='c'
        >>> ion='iii' # C II
        >>> c_ii_rc_data=atomneb.read_aeff_collection(atom_rc_all_file, atom, ion)
        >>> temperature=np.float64(10000.0)
        >>> density=np.float64(5000.0)
        >>> c_ii_6151_flux = 0.028
        >>> wavelength=6151.43
        >>> abund_c_ii=pyequib.calc_abund_c_ii_rl(temperature=temperature, density=density,
        >>>                               wavelength=wavelength, line_flux=c_ii_6151_flux,
        >>>                               c_ii_rc_data=c_ii_rc_data, h_i_aeff_data=h_i_aeff_data)
        >>> print('N(C^2+)/N(H+):', abund_c_ii)
           N(C^2+)/N(H+):    0.00063404650

    :Categories:
      Abundance Analysis, Recombination Lines

    :Dirs:
     ./
         Main routines

    :Author:
      Ashkbiz Danehkar

    :Copyright:
      This library is released under a GNU General Public License.

    :Version:
      0.3.0

    :History:
        Based on recombination coefficients for C II lines from
        Davey et al. 2000A&AS..142...85D.

        Adopted from MOCASSIN, Ercolano et al. 2005MNRAS.362.1038E.

        02/2003, Yong Zhang, added to MOCASSIN.

        10/05/2013, A. Danehkar, Translated to IDL code.

        15/04/2017, A. Danehkar, Integration with AtomNeb.

        10/07/2019, A. Danehkar, Made a new function calc_emiss_c_ii_rl()
                         for calculating line emissivities and separated it
                         from calc_abund_c_ii_rl().

        03/10/2020, A. Danehkar, Transferred from IDL to Python. 
   """
   # ciiRLstructure ={Wave:np.float64(0.0), Int:np.float64(0.0), Obs:np.float64(0.0), Abundance:np.float64(0.0)}
   if (temperature is not None) == 0:   
      print('Temperature is not set')
      return 0
   if (density is not None) == 0:   
      print('Density is not set')
      return 0
   if (c_ii_rc_data is not None) == 0:   
      print('C II recombination coefficients (c_ii_rc_data) are not set')
      return 0
   if (h_i_aeff_data is not None) == 0:   
      print('H I recombination coefficients (h_i_aeff_data) are not set')
      return 0
   if (wavelength is not None) == 0:   
      print('Wavelength is not given')
      return 0
   if (line_flux is not None) == 0:   
      print('Line flux intensity (line_flux) is not given')
      return 0
   if ((temperature <= 0.e0) | (density <= 0.e0)):
      print('temperature = ', temperature, ', density = ', density)
      return 0
   
   emissivity_hbeta = calc_emiss_h_beta(temperature=temperature, density=density, h_i_aeff_data=h_i_aeff_data)
   emissivity = calc_emiss_c_ii_rl(temperature=temperature, density=density, wavelength=wavelength, c_ii_rc_data=c_ii_rc_data)
   abund = (emissivity_hbeta / emissivity) * np.float64(line_flux / 100.0)
   
   return abund

def calc_abund_c_iii_rl(temperature=None, density=None,
                        wavelength=None, line_flux=None,
                        c_iii_rc_data=None, h_i_aeff_data=None):
   """
         This function determines the ionic abundance from the observed
         flux intensity for the given wavelength of C III recombination line
         by using the recombination coefficients from
         Pequignot et al. 1991A&A...251..680P.

     :Returns:
        type=double. This function returns the ionic abundanc.

     :Keywords:
         temperature   :     in, required, type=float
                             electron temperature
         density       :     in, required, type=float
                             electron density
         wavelength    :     in, required, type=float
                             Line Wavelength in Angstrom
         line_flux     :     in, required, type=float
                             line flux intensity
         c_iii_rc_data :     in, required, type=array/object
                             C III recombination coefficients
         h_i_aeff_data :     in, required, type=array/object
                             H I recombination coefficients

    :Examples:
       For example::

         >>> import pyequib
         >>> import atomneb
         >>> import os
         >>> base_dir = '../externals/atomneb/'
         >>> data_rc_dir = os.path.join('atomic-data-rc')
         >>> atom_rc_ppb91_file=os.path.join(base_dir,data_rc_dir, 'rc_PPB91.fits')
         >>> atom_rc_sh95_file= os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
         >>> atom='h'
         >>> ion='ii' # H I
         >>> h_i_rc_data=atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)
         >>> h_i_aeff_data=h_i_rc_data['aeff'][0]
         >>> atom='c'
         >>> ion='iv' # C III
         >>> c_iii_rc_data=atomneb.read_aeff_ppb91(atom_rc_ppb91_file, atom, ion)
         >>> temperature=np.float64(10000.0)
         >>> density=np.float64(5000.0)
         >>> c_iii_4647_flux = 0.107
         >>> wavelength=4647.42
         >>> abund_c_iii=pyequib.calc_abund_c_iii_rl(temperature=temperature, density=density,
         >>>                                 wavelength=wavelength, line_flux=c_iii_4647_flux,
         >>>                                 c_iii_rc_data=c_iii_rc_data, h_i_aeff_data=h_i_aeff_data)
         >>> print('N(C^3+)/N(H+):', abund_c_iii)
            N(C^3+)/N(H+):    0.00017502840

     :Categories:
       Abundance Analysis, Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.3.0

     :History:
         Based on effective radiative recombination coefficients for C III lines from
         Pequignot, Petitjean, Boisson, C. 1991A&A...251..680P.

         18/05/2013, A. Danehkar, Translated to IDL code.

         06/04/2017, A. Danehkar, Integration with AtomNeb.

         10/07/2019, A. Danehkar, Made a new function calc_emiss_c_iii_rl()
                          for calculating line emissivities and separated it
                          from calc_abund_c_iii_rl().

         03/10/2020, A. Danehkar, Transferred from IDL to Python.
   """
   # ciiiRLstructure ={Wave:np.float64(0.0), Int:np.float64(0.0), Obs:np.float64(0.0), Abundance:np.float64(0.0)}
   if (temperature is not None) == 0:   
      print('Temperature is not set')
      return 0
   if (density is not None) == 0:   
      print('Density is not set')
      return 0
   if (c_iii_rc_data is not None) == 0:   
      print('C III recombination coefficients (c_iii_rc_data) are not set')
      return 0
   if (h_i_aeff_data is not None) == 0:   
      print('H I recombination coefficients (h_i_aeff_data) are not set')
      return 0
   if (wavelength is not None) == 0:   
      print('Wavelength is not given')
      return 0
   if (line_flux is not None) == 0:   
      print('Line flux intensity (line_flux) is not given')
      return 0
   if ((temperature <= 0.e0) | (density <= 0.e0)):
      print('temperature = ', temperature, ', density = ', density)
      return 0
   
   emissivity_hbeta = calc_emiss_h_beta(temperature=temperature, density=density, h_i_aeff_data=h_i_aeff_data)
   emissivity = calc_emiss_c_iii_rl(temperature=temperature, density=density, wavelength=wavelength, c_iii_rc_data=c_iii_rc_data)
   abund = (emissivity_hbeta / emissivity) * np.float64(line_flux / 100.0)
   
   return abund

def calc_abund_n_ii_rl(temperature=None, density=None,
                       wavelength=None, line_flux=None,
                       n_ii_rc_br=None, n_ii_rc_data=None,
                       h_i_aeff_data=None):
   """
        This function determines the ionic abundance from the observed
        flux intensity for the given wavelength of N II recombination line
        by using the recombination coefficients from
        Escalante & Victor 1990ApJS...73..513E.

    :Returns:
       type=double. This function returns the ionic abundanc.

    :Keywords:
        temperature   :     in, required, type=float
                            electron temperature
        density       :     in, required, type=float
                            electron density
        wavelength    :     in, required, type=float
                            Line Wavelength in Angstrom
        line_flux     :     in, required, type=float
                            line flux intensity
        n_ii_rc_br    :     in, required, type=array/object
                            N II branching ratios (Br)
        n_ii_rc_data  :     in, required, type=array/object
                            N II recombination coefficients
        h_i_aeff_data :     in, required, type=array/object
                            H I recombination coefficients

    :Examples:
       For example::

        >>> import pyequib
        >>> import atomneb
        >>> import os
        >>> base_dir = '../externals/atomneb/'
        >>> data_rc_dir = os.path.join('atomic-data-rc')
        >>> atom_rc_all_file= os.path.join(base_dir,data_rc_dir, 'rc_collection.fits')
        >>> atom_rc_sh95_file= os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
        >>> atom='h'
        >>> ion='ii' # H I
        >>> h_i_rc_data=atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)
        >>> h_i_aeff_data=h_i_rc_data['aeff'][0]
        >>> atom='n'
        >>> ion='iii' # N II
        >>> n_ii_rc_data=atomneb.read_aeff_collection(atom_rc_all_file, atom, ion)
        >>> n_ii_rc_data_br=atomneb.read_aeff_collection(atom_rc_all_file, atom, ion, br=True)
        >>> temperature=np.float64(10000.0)
        >>> density=np.float64(5000.0)
        >>> n_ii_4442_flux = 0.017
        >>> wavelength=4442.02
        >>> abund_n_ii=pyequib.calc_abund_n_ii_rl(temperature=temperature, density=density,
        >>>                               wavelength=wavelength, line_flux=n_ii_4442_flux,
        >>>                               n_ii_rc_br=n_ii_rc_data_br, n_ii_rc_data=n_ii_rc_data,
        >>>                               h_i_aeff_data=h_i_aeff_data)
        >>> print('N(N^2+)/N(H+):', abund_n_ii)
           N(N^2+)/N(H+):   0.00069297541

    :Categories:
      Abundance Analysis, Recombination Lines

    :Dirs:
     ./
         Main routines

    :Author:
      Ashkbiz Danehkar

    :Copyright:
      This library is released under a GNU General Public License.

    :Version:
      0.3.0

    :History:
        Based on Effective recombination coefficients for N II lines from
        Escalante & Victor 1990ApJS...73..513E.

        Adopted from MIDAS Rnii script written by X.W.Liu.

        Revised based on scripts by Yong Zhang added to MOCASSIN, 02/2003
                          Ercolano et al. 2005MNRAS.362.1038E.

        10/05/2013, A. Danehkar, Translated to IDL code.

        25/04/2017, A. Danehkar, Integration with AtomNeb.

        10/07/2019, A. Danehkar, Made a new function calc_emiss_n_ii_rl()
                         for calculating line emissivities and separated it
                         from calc_abund_n_ii_rl().

        03/10/2020, A. Danehkar, Transferred from IDL to Python.
   """

#  niiRLstructure ={Wave:np.float64(0.0), #REAL*8
#              Int:np.float64(0.0),
#              Obs:np.float64(0.0),
#              abundance:np.float64(0.0),
#              g1:long(0), #INTEGER
#              g2:long(0), #INTEGER
#              Mult1:'', #CHARACTER*7
#              Term1:'', #CHARACTER*9
#              Term2:'' #CHARACTER*9
#              }
   if (temperature is not None) == 0:   
      print('Temperature is not set')
      return 0
   if (density is not None) == 0:   
      print('Density is not set')
      return 0
   if (n_ii_rc_data is not None) == 0:   
      print('N II recombination coefficients (n_ii_rc_data) are not set')
      return 0
   if (n_ii_rc_br is not None) == 0:   
      print('N II branching ratios (n_ii_rc_br) are not set')
      return 0
   if (h_i_aeff_data is not None) == 0:   
      print('H I recombination coefficients (h_i_aeff_data) are not set')
      return 0
   if (wavelength is not None) == 0:   
      print('Wavelength is not given')
      return 0
   if (line_flux is not None) == 0:   
      print('Line flux intensity (line_flux) is not given')
      return 0
   if ((temperature <= 0.e0) | (density <= 0.e0)):
      print('temperature = ', temperature, ', density = ', density)
      return 0
   
   emissivity_hbeta = calc_emiss_h_beta(temperature=temperature, density=density, h_i_aeff_data=h_i_aeff_data)
   emissivity = calc_emiss_n_ii_rl(temperature=temperature, density=density, wavelength=wavelength, n_ii_rc_br=n_ii_rc_br, n_ii_rc_data=n_ii_rc_data)
   abund = (emissivity_hbeta / emissivity) * np.float64(line_flux / 100.0)
   
   return abund

def calc_abund_n_iii_rl(temperature=None, density=None,
                        wavelength=None, line_flux=None,
                        n_iii_rc_data=None, h_i_aeff_data=None):
   """
         This function determines the ionic abundance from the observed
         flux intensity for the given wavelength of N III recombination line
         by using the recombination coefficients from
         Pequignot et al. 1991A&A...251..680P.

     :Returns:
        type=double. This function returns the ionic abundanc.

     :Keywords:
         temperature   :     in, required, type=float
                             electron temperature
         density       :     in, required, type=float
                             electron density
         wavelength    :     in, required, type=float
                             Line Wavelength in Angstrom
         line_flux     :     in, required, type=float
                             line flux intensity
         n_iii_rc_data  :     in, required, type=array/object
                             N III recombination coefficients
         h_i_aeff_data :     in, required, type=array/object
                             H I recombination coefficients

    :Examples:
       For example::

         >>> import pyequib
         >>> import atomneb
         >>> import os
         >>> base_dir = '../externals/atomneb/'
         >>> data_rc_dir = os.path.join('atomic-data-rc')
         >>> atom_rc_ppb91_file=os.path.join(base_dir,data_rc_dir, 'rc_PPB91.fits')
         >>> atom_rc_sh95_file= os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
         >>> atom='h'
         >>> ion='ii' # H I
         >>> h_i_rc_data=atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)
         >>> h_i_aeff_data=h_i_rc_data['aeff'][0]
         >>> atom='n'
         >>> ion='iv' # N III
         >>> n_iii_rc_data=atomneb.read_aeff_ppb91(atom_rc_ppb91_file, atom, ion)
         >>> temperature=np.float64(10000.0)
         >>> density=np.float64(5000.0)
         >>> n_iii_4641_flux = 0.245
         >>> wavelength=4640.64
         >>> abund_n_iii=pyequib.calc_abund_n_iii_rl(temperature=temperature, density=density,
         >>>                                 wavelength=wavelength, line_flux=n_iii_4641_flux,
         >>>                                 n_iii_rc_data=n_iii_rc_data, h_i_aeff_data=h_i_aeff_data)
         >>> print('N(N^3+)/N(H+):', abund_n_iii)
            N(N^3+)/N(H+):    6.3366175e-05

     :Categories:
       Abundance Analysis, Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.3.0

     :History:
         Based on  effective radiative recombination coefficients for N III lines from
         Pequignot, Petitjean, Boisson, C. 1991A&A...251..680P.

         10/05/2013, A. Danehkar, IDL code written.

         20/04/2017, A. Danehkar, Integration with AtomNeb.

         10/07/2019, A. Danehkar, Made a new function calc_emiss_n_iii_rl()
                          for calculating line emissivities and separated it
                          from calc_abund_n_iii_rl().

         03/10/2020, A. Danehkar, Transferred from IDL to Python.
   """
   # niiiRLstructure ={Wave:np.float64(0.0), Int:np.float64(0.0), Obs:np.float64(0.0), Abundance:np.float64(0.0)}

   if (temperature is not None) == 0:   
      print('Temperature is not set')
      return 0
   if (density is not None) == 0:   
      print('Density is not set')
      return 0
   if (n_iii_rc_data is not None) == 0:   
      print('N III recombination coefficients (n_iii_rc_data) are not set')
      return 0
   if (h_i_aeff_data is not None) == 0:   
      print('H I recombination coefficients (h_i_aeff_data) are not set')
      return 0
   if (wavelength is not None) == 0:   
      print('Wavelength is not given')
      return 0
   if (line_flux is not None) == 0:   
      print('Line flux intensity (line_flux) is not given')
      return 0
   if ((temperature <= 0.e0) | (density <= 0.e0)):
      print('temperature = ', temperature, ', density = ', density)
      return 0
   
   emissivity_hbeta = calc_emiss_h_beta(temperature=temperature, density=density, h_i_aeff_data=h_i_aeff_data)
   emissivity = calc_emiss_n_iii_rl(temperature=temperature, density=density, wavelength=wavelength, n_iii_rc_data=n_iii_rc_data)
   abund = (emissivity_hbeta / emissivity) * np.float64(line_flux / 100.0)
   
   return abund

def calc_abund_o_ii_rl(temperature=None, density=None,
                       wavelength=None, line_flux=None,
                       o_ii_rc_br=None, o_ii_rc_data=None,
                       h_i_aeff_data=None):
   """
        This function determines the ionic abundance from the observed
        flux intensity for the given wavelength of O II recombination line
        by using the recombination coefficients from
        Storey 1994A&A...282..999S and Liu et al. 1995MNRAS.272..369L.

    :Returns:
       type=double. This function returns the ionic abundanc.

    :Keywords:
        temperature   :     in, required, type=float
                            electron temperature
        density       :     in, required, type=float
                            electron density
        wavelength    :     in, required, type=float
                            Line Wavelength in Angstrom
        line_flux     :     in, required, type=float
                            line flux intensity
        o_ii_rc_br    :     in, required, type=array/object
                            O II branching ratios (Br)
        o_ii_rc_data  :     in, required, type=array/object
                            O II recombination coefficients
        h_i_aeff_data :     in, required, type=array/object
                            H I recombination coefficients

    :Examples:
       For example::

        >>> import pyequib
        >>> import atomneb
        >>> import os
        >>> base_dir = '../externals/atomneb/'
        >>> data_rc_dir = os.path.join('atomic-data-rc')
        >>> atom_rc_all_file= os.path.join(base_dir,data_rc_dir, 'rc_collection.fits')
        >>> atom_rc_sh95_file= os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
        >>> atom='h'
        >>> ion='ii' # H I
        >>> h_i_rc_data=atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)
        >>> h_i_aeff_data=h_i_rc_data['aeff'][0]
        >>> atom='o'
        >>> ion='iii' # O II
        >>> o_ii_rc_data=atomneb.read_aeff_collection(atom_rc_all_file, atom, ion)
        >>> o_ii_rc_data_br=atomneb.read_aeff_collection(atom_rc_all_file, atom, ion, br=True)
        >>> temperature=np.float64(10000.0)
        >>> density=np.float64(5000.0)
        >>> o_ii_4614_flux = 0.009
        >>> wavelength=4613.68
        >>> abund_o_ii=pyequib.calc_abund_o_ii_rl(temperature=temperature, density=density,
        >>>                               wavelength=wavelength, line_flux=o_ii_4614_flux,
        >>>                               o_ii_rc_br=o_ii_rc_data_br, o_ii_rc_data=o_ii_rc_data,
        >>>                               h_i_aeff_data=h_i_aeff_data)
        >>> print('N(O^2+)/N(H+):', abund_o_ii)
           N(O^2+)/N(H+):    0.0018886330

    :Categories:
      Abundance Analysis, Recombination Lines

    :Dirs:
     ./
         Main routines

    :Author:
      Ashkbiz Danehkar

    :Copyright:
      This library is released under a GNU General Public License.

    :Version:
      0.3.0

    :History:
        Based on recombination coefficients for O II lines from
        Storey 1994A&A...282..999S and Liu et al. 1995MNRAS.272..369L.

        Adopted from MIDAS script Roii.prg written by X.W.Liu.

        Revised based on scripts by Yong Zhang added to MOCASSIN, 02/2003
                          Ercolano et al. 2005MNRAS.362.1038E.

        10/05/2013, A. Danehkar, Translated to IDL code.

        25/04/2017, A. Danehkar, Integration with AtomNeb.

        10/07/2019, A. Danehkar, Made a new function calc_emiss_o_ii_rl()
                         for calculating line emissivities and separated it
                         from calc_abund_o_ii_rl().

        03/10/2020, A. Danehkar, Transferred from IDL to Python.
   """

#  oiiRLstructure ={Wave:np.float64(0.0), #REAL*8
#              Int:np.float64(0.0),
#              Obs:np.float64(0.0),
#              abundance:np.float64(0.0),
#              g1:long(0), #INTEGER
#              g2:long(0), #INTEGER
#              Mult1:'', #CHARACTER*7
#              Term1:'', #CHARACTER*9
#              Term2:'' #CHARACTER*9
#              }
   if (temperature is not None) == 0:   
      print('Temperature is not set')
      return 0
   if (density is not None) == 0:   
      print('Density is not set')
      return 0
   if (o_ii_rc_data is not None) == 0:   
      print('O II recombination coefficients (o_ii_rc_data) are not set')
      return 0
   if (o_ii_rc_br is not None) == 0:   
      print('O II branching ratios (o_ii_rc_br) are not set')
      return 0
   if (h_i_aeff_data is not None) == 0:   
      print('H I recombination coefficients (h_i_aeff_data) are not set')
      return 0
   if (wavelength is not None) == 0:   
      print('Wavelength is not given')
      return 0
   if (line_flux is not None) == 0:   
      print('Line flux intensity (line_flux) is not given')
      return 0
   if ((temperature <= 0.e0) | (density <= 0.e0)):
      print('temperature = ', temperature, ', density = ', density)
      return 0
   
   emissivity_hbeta = calc_emiss_h_beta(temperature=temperature, density=density, h_i_aeff_data=h_i_aeff_data)
   emissivity = calc_emiss_o_ii_rl(temperature=temperature, density=density, wavelength=wavelength, o_ii_rc_br=o_ii_rc_br, o_ii_rc_data=o_ii_rc_data)
   abund = (emissivity_hbeta / emissivity) * np.float64(line_flux / 100.0)
   
   return abund

def calc_abund_ne_ii_rl(temperature=None, density=None,
                        wavelength=None, line_flux=None,
                        ne_ii_rc_data=None, h_i_aeff_data=None):
   """
         This function determines the ionic abundance from the observed
         flux intensity for the given wavelength of Ne II recombination line
         by using the recombination coefficients from
         Kisielius et al. (1998) & Storey (unpublished).

     :Returns:
        type=double. This function returns the ionic abundanc.

     :Keywords:
         temperature   :     in, required, type=float
                             electron temperature
         density       :     in, required, type=float
                             electron density
         wavelength    :     in, required, type=float
                             Line Wavelength in Angstrom
         line_flux     :     in, required, type=float
                             line flux intensity
         ne_ii_rc_data  :    in, required, type=array/object
                             Ne II recombination coefficients
         h_i_aeff_data :     in, required, type=array/object
                             H I recombination coefficients

    :Examples:
       For example::

         >>> import pyequib
         >>> import atomneb
         >>> import os
         >>> base_dir = '../externals/atomneb/'
         >>> data_rc_dir = os.path.join('atomic-data-rc')
         >>> atom_rc_all_file= os.path.join(base_dir,data_rc_dir, 'rc_collection.fits')
         >>> atom_rc_sh95_file= os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
         >>> atom='h'
         >>> ion='ii' # H I
         >>> h_i_rc_data=atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)
         >>> h_i_aeff_data=h_i_rc_data['aeff'][0]
         >>> atom='ne'
         >>> ion='iii' # Ne II
         >>> ne_ii_rc_data=atomneb.read_aeff_collection(atom_rc_all_file, atom, ion)
         >>> temperature=np.float64(10000.0)
         >>> density=np.float64(5000.0)
         >>> ne_ii_3777_flux = 0.056
         >>> wavelength=3777.14
         >>> abund_ne_ii=pyequib.calc_abund_ne_ii_rl(temperature=temperature, density=density,
         >>>                                 wavelength=wavelength, line_flux=ne_ii_3777_flux,
         >>>                                 ne_ii_rc_data=ne_ii_rc_data, h_i_aeff_data=h_i_aeff_data)
         >>> print('N(Ne^2+)/N(H+):', abund_ne_ii)
            N(Ne^2+)/N(H+):    0.00043376850

     :Categories:
       Abundance Analysis, Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.3.0

     :History:
         Based on effective radiative recombination coefficients for Ne II lines
         from Kisielius et al. 1998A&AS..133..257K & Storey (unpublished).

         Adopted from MOCASSIN, Ercolano et al. 2005MNRAS.362.1038E.

         02/2003, Yong Zhang, scripts added to MOCASSIN.

         14/05/2013, A. Danehkar, Translated to IDL code.

         10/04/2017, A. Danehkar, Integration with AtomNeb.

         10/07/2019, A. Danehkar, Made a new function calc_emiss_ne_ii_rl()
                          for calculating line emissivities and separated it
                          from calc_abund_ne_ii_rl().

         03/10/2020, A. Danehkar, Transferred from IDL to Python.
   """
   # neiiRLstructure ={Wave:np.float64(0.0), Int:np.float64(0.0), Obs:np.float64(0.0), Abundance:np.float64(0.0)}
   if (temperature is not None) == 0:   
      print('Temperature is not set')
      return 0
   if (density is not None) == 0:   
      print('Density is not set')
      return 0
   if (ne_ii_rc_data is not None) == 0:   
      print('Ne II recombination coefficients (ne_ii_rc_data) are not set')
      return 0
   if (h_i_aeff_data is not None) == 0:   
      print('H I recombination coefficients (h_i_aeff_data) are not set')
      return 0
   if (wavelength is not None) == 0:   
      print('Wavelength is not given')
      return 0
   if (line_flux is not None) == 0:   
      print('Line flux intensity (line_flux) is not given')
      return 0
   if ((temperature <= 0.e0) | (density <= 0.e0)):
      print('temperature = ', temperature, ', density = ', density)
      return 0
   
   emissivity_hbeta = calc_emiss_h_beta(temperature=temperature, density=density, h_i_aeff_data=h_i_aeff_data)
   emissivity = calc_emiss_ne_ii_rl(temperature=temperature, density=density, wavelength=wavelength, ne_ii_rc_data=ne_ii_rc_data)
   abund = (emissivity_hbeta / emissivity) * np.float64(line_flux / 100.0)
   
   return abund

def redlaw(wavelength, ext_law=None, rv=None, fmlaw=None):
   """
         This function determines the reddening law function of the line at the given wavelength
         for the used extinction law.

     :Returns:
        type=double/array. This function returns the reddening law function value for the given wavelength.

     :Params:
         wavelength :  in, required, type=float/array
                       Wavelength in Angstrom

     :Keywords:
        ext_law  :  in, optional, type=string, default='GAL'
                    the extinction law:

                    'GAL' for Howarth Galactic.

                    'GAL2' for Savage and Mathis.

                    'CCM' for CCM galactic.

                    'JBK' for Whitford, Seaton, Kaler.

                    'FM' for Fitxpatrick.

                    'SMC' for Prevot SMC.

                    'LMC' for Howarth LMC.

        rv       :  in, optional, type=float, default=3.1
                    the optical total-to-selective extinction ratio, RV = A(V)/E(B-V).

        fmlaw    :  in, optional, type=string, default='GAL'
                    the fmlaw keyword is used only in the redlaw_fm function:

                    'GAL' for  the default fit parameters for the R-dependent
                               Galactic extinction curve from Fitzpatrick & Massa
                               (Fitzpatrick, 1999, PASP, 111, 63).

                    'LMC2' for the fit parameters are those determined for
                                  reddening the LMC2 field (inc. 30 Dor)
                                  from Misselt et al.  (1999, ApJ, 515, 128).

                    'AVGLMC' for  the fit parameters are those determined for
                                  reddening in the general Large Magellanic Cloud (LMC)
                                  field by Misselt et al.  (1999, ApJ, 515, 128).

    :Examples:
       For example::

         >>> import pyequib
         >>> wavelength=6563.0
         >>> r_v=3.1
         >>> fl=pyequib.redlaw(wavelength, rv=r_v)
         >>> print('fl(6563)', fl)
            fl(6563)     -0.32013816

     :Categories:
       Interstellar Extinction

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.3.0

     :History:
         Originally from IRAF STSDAS SYNPHOT redlaw.x, ebmvxfunc.x

         31/08/2012, A. Danehkar, Converted to IDL code.

         03/10/2020, A. Danehkar, Transferred from IDL to Python.
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
      fl = redlaw_jbk(wavelength)
   elif _expr == 'FM':
      fl = redlaw_fm(wavelength, rv=rv, fmlaw=fmlaw)
   elif _expr == 'SMC':
      fl = redlaw_smc(wavelength)
   elif _expr == 'LMC':
      fl = redlaw_lmc(wavelength)
   else:   
      print("ext_law cannnot find")
   
   return fl

def redlaw_gal(wavelength, rv=None):
   """
       This function determines the reddening law function of the line at the given wavelength
       for Galactic Seaton1979+Howarth1983+CCM1983.

    :Returns:
       type=double/array. This function returns the reddening law function value(s) for the given wavelength(s).

    :Params:
       wavelength :  in, required, type=float
                      Wavelength in Angstrom

    :Keywords:
       rv       :  in, optional, type=float, default=3.1
                   the optical total-to-selective extinction ratio, RV = A(V)/E(B-V).

    :Examples:
       For example::

        >>> import pyequib
        >>> wavelength=6563.0
        >>> r_v=3.1
        >>> fl=pyequib.redlaw_gal(wavelength, rv=r_v)
        >>> print('fl(6563)', fl)
           fl(6563)     -0.32013816

    :Categories:
      Interstellar Extinction

    :Dirs:
     ./
         Subroutines

    :Author:
      Ashkbiz Danehkar

    :Copyright:
      This library is released under a GNU General Public License.

    :Version:
      0.3.0

    :History:
        Based on the UV Formulae from Seaton 1979, MNRAS, 187, 73
        1979MNRAS.187P..73S, the opt/NIR from Howarth 1983, MNRAS, 203, 301
        the FIR from Cardelli, Clayton and Mathis 1989, ApJ, 345, 245
        1989ApJ...345..245C

        Originally from IRAF STSDAS SYNPHOT ebmvxfunc.x, pyneb.extinction

        31/08/2012, A. Danehkar, Converted to IDL code.

        03/10/2020, A. Danehkar, Transferred from IDL to Python.
   """
   # Tabulated inverse wavelengths in microns:
   xtable = np.array([0., 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                      2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7])
   etable = np.array([0., 1.36, 1.64, 1.84, 2.04, 2.24, 2.44, 2.66, 2.88, 3.14,
                      3.36, 3.56, 3.77, 3.96, 4.15, 4.26, 4.40, 4.52, 4.64])

   if hasattr(wavelength, "__len__"):
      npts = len(wavelength)
      extl = np.zeros(npts)
   else:
      npts = 1
      extl = np.int32(0)
   if (rv is not None):   
      r_v = rv
   else:   
      r_v = 3.1
   for pix in range(0, npts):
   # Convert wavelength in angstroms to 1/microns
      if hasattr(wavelength, "__len__"):
          wavel = wavelength[pix]
      else:
          wavel = wavelength
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

def redlaw_gal2(wavelength):
   """
       This function determines the reddening law function of the line at the given wavelength
       for Galactic Savage & Mathis 1979.

    :Returns:
       type=double/array. This function returns the reddening law function value(s) for the given wavelength(s).

    :Params:
        wavelength :  in, required, type=float
                      Wavelength in Angstrom

    :Examples:
       For example::

        >>> import pyequib
        >>> wavelength=6563.0
        >>> fl=pyequib.redlaw_gal2(wavelength)
        >>> print('fl(6563)', fl)
           fl(6563)     -0.30925984

    :Categories:
      Interstellar Extinction

    :Dirs:
     ./
        Subroutines

    :Author:
      Ashkbiz Danehkar

    :Copyright:
      This library is released under a GNU General Public License.

    :Version:
      0.3.0

    :History:
        Based on Savage & Mathis 1979, ARA&A, vol. 17, 73-111

        Originally from IRAF STSDAS SYNPHOT ebmvxfunc.x

        20/09/1994, R. A. Shaw, Initial IRAF implementation.

        04/03/1995, R. A. Shaw, Return A(lambda)/A(V) instead.

        31/08/2012, A. Danehkar, Converted to IDL code.

        03/10/2020, A. Danehkar, Transferred from IDL to Python.
   """
   # Tabulated inverse wavelengths in microns:
   xtable = np.array([0.00, 0.29, 0.45, 0.80, 1.11, 1.43, 1.82, 2.27, 2.50, 2.91, 3.65,
                      4.00, 4.17, 4.35, 4.57, 4.76, 5.00, 5.26, 5.56, 5.88, 6.25, 6.71,
                      7.18, 8.00, 8.50, 9.00, 9.50, 10.00])
   
   #  Tabulated extinction function, A(lambda)/E(B-V):
   etable = np.array([0.00, 0.16, 0.38, 0.87, 1.50, 2.32, 3.10, 4.10, 4.40, 4.90, 6.20,
                      7.29, 8.00, 8.87, 9.67, 9.33, 8.62, 8.00, 7.75, 7.87, 8.12, 8.15,
                      8.49, 9.65, 10.55, 11.55, 12.90, 14.40])

   if hasattr(wavelength, "__len__"):
     npts = len(wavelength)
     extl = np.zeros(npts)
   else:
     npts = 1
     extl = np.int32(0)
   for pix in range(0, npts):
      if hasattr(wavelength, "__len__"):
         wavel=wavelength[pix]
      else:
         wavel=wavelength
      if (wavel < 1000.0):   
         print("redlaw_smc: Invalid wavelength")
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

def redlaw_ccm(wavelength, rv=None):
   """
        This function determines the reddening law function of Cardelli, Clayton & Mathis.

     :Returns:
        type=double/array. This function returns the reddening law function value for the given wavelength.

     :Params:
         wavelength :  in, required, type=float/array
                       Wavelength in Angstrom

     :Keywords:
        RV       :  in, optional, type=float, default=3.1
                    the optical total-to-selective extinction ratio, RV = A(V)/E(B-V).

    :Examples:
       For example::

         >>> import pyequib
         >>> wavelength=6563.0
         >>> r_v=3.1
         >>> fl=pyequib..redlaw_ccm(wavelength, rv=r_v)
         >>> print('fl(6563)', fl)
            fl(6563)     -0.29756615

     :Categories:
       Interstellar Extinction

     :Dirs:
      ./
          Subroutines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.3.0

     :History:
         Based on Formulae by Cardelli, Clayton & Mathis 1989, ApJ 345, 245-256.
         1989ApJ...345..245C

         Originally from IRAF STSDAS SYNPHOT redlaw.x

         18/05/1993, R. A. Shaw, Initial IRAF implementation, based upon CCM module
             in onedspec.deredden.

         31/08/2012, A. Danehkar, Converted to IDL code.

         03/10/2020, A. Danehkar, Transferred from IDL to Python.
   """
   if hasattr(wavelength, "__len__"):
     npts = len(wavelength)
     extl = np.zeros(npts)
   else:
     npts = 1
     extl = np.int32(0)
   if (rv is not None):   
      r_v = rv
   else:   
      r_v = 3.1
   for pix in range(0, npts):
      if hasattr(wavelength, "__len__"):
         wavel=wavelength[pix]
      else:
         wavel=wavelength
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
               a = 1 + y * (0.17699 + y * (-0.50447 +
                            y * (-0.02427 + y * (0.72085 + y * (0.01979 + y * (-0.77530 + y * 0.32999))))))
               b = y * (1.41338 + y * (2.28305 +
                            y * (1.07233 + y * (-5.38434 + y * (-0.62251 + y * (5.30260 - y * 2.09002))))))
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

def redlaw_jbk(wavelength):
   """
       This function determines the reddening law function for Galactic Whitford1958 + Seaton1977 + Kaler1976.

    :Returns:
       type=double/array. This function returns the reddening law function value(s) for the given wavelength(s).

    :Params:
        wavelength :  in, required, type=float
                      Wavelength in Angstrom

    :Examples:
       For example::

        >>> import pyequib
        >>> wavelength=6563.0
        >>> fl=pyequib.redlaw_jbk(wavelength)
        >>> print('fl(6563)', fl)
           fl(6563)     -0.33113684

    :Categories:
      Interstellar Extinction

    :Dirs:
     ./
         Subroutines

    :Author:
      Ashkbiz Danehkar

    :Copyright:
      This library is released under a GNU General Public License.

    :Version:
      0.3.0

    :History:
        Based on Whitford (1958), extended to the UV by Seaton (1977),
        adapted by Kaler (1976).

        Originally from IRAF STSDAS SYNPHOT redlaw.x

        13/05/1993, R. A. Shaw, Initial IRAF implementation.

        31/08/2012, A. Danehkar, Converted to IDL code.

        03/10/2020, A. Danehkar, Transferred from IDL to Python.
   """
   # Tabulated wavelengths, Angstroms:
   refw = np.array([1150., 1200., 1250., 1300., 1350., 1400., 1450., 1500., 1550.,
                    1600., 1650., 1700., 1750., 1800., 1850., 1900., 1950., 2000.,
                    2050., 2100., 2150., 2200., 2250., 2300., 2350., 2400., 2450.,
                    2500., 2550., 2600., 2650., 2700., 2750., 2800., 2850., 2900.,
                    2950., 3000., 3050., 3100., 3333., 3500., 3600., 3700., 3800.,
                    3900., 4000., 4100., 4200., 4300., 4400., 4500., 4600., 4700.,
                    4800., 4861.3, 5000., 5100., 5200., 5300., 5400., 5500., 5600.,
                    5700., 5800., 5900., 6000., 6100., 6200., 6300., 6400., 6500.,
                    6600., 6700., 6800., 6900., 7000., 7200., 7400., 7600., 7800.,
                    8000., 8200., 8400., 8600., 8800., 9000., 9500., 10000., 11000.,
                    12000., 14000., 16000., 20000., 1.e+6])
   
   #  Tabulated extinction function:
   extab = np.array([1.96, 1.78, 1.61, 1.49, 1.37, 1.29, 1.24, 1.20, 1.20, 1.20, 1.17,
                     1.13, 1.11, 1.10, 1.12, 1.17, 1.25, 1.35, 1.45, 1.53, 1.60, 1.62,
                     1.52, 1.40, 1.28, 1.17, 1.06, 0.98, 0.9, 0.84, 0.77, 0.72, 0.68,
                     0.64, 0.60, 0.57, 0.53, 0.51, 0.48, 0.46, 0.385, 0.358, 0.33,
                     0.306, 0.278, 0.248, 0.220, 0.195, 0.168, 0.143, 0.118, 0.095,
                     0.065, 0.040, 0.015, 0.000, -0.030, -0.055, -0.078, -0.10, -0.121,
                     -0.142, -0.164, -0.182, -0.201, -0.220, -0.238, -0.254, -0.273,
                     -0.291, -0.306, -0.321, -0.337, -0.351, -0.365, -0.377, -0.391,
                     -0.416, -0.441, -0.465, -0.490, -0.510, -0.529, -0.548, -0.566,
                     -0.582, -0.597, -0.633, -0.663, -0.718, -0.763, -0.840, -0.890,
                     -0.960, -1.000])
   
   xtable = 10000.e+0 / refw
   if hasattr(wavelength, "__len__"):
     npts = len(wavelength)
     extl = np.zeros(npts)
   else:
     npts = 1
     extl = np.int32(0)
   for pix in range(0, npts):
      if hasattr(wavelength, "__len__"):
         wavel=wavelength[pix]
      else:
         wavel=wavelength
      if (wavel < 1000.0):   
         print("redlaw_smc: Invalid wavelength")
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

def redlaw_fm(wavelength, rv=None, fmlaw=None):
   """
        This function determines the reddening law function by Fitzpatrick & Massa
        for the line at the given wavelength.

     :Returns:
        type=double/array. This function returns the reddening law function value for the given wavelength.

     :Params:
         wavelength :  in, required, type=float/array
                       Wavelength in Angstrom

     :Keywords:
        RV       :  in, optional, type=float, default=3.1
                    the optical total-to-selective extinction ratio, RV = A(V)/E(B-V).

        fmlaw    :  in, optional, type=string, default='GAL'
                    the fmlaw keyword is used only in the redlaw_fm function:

                    'GAL' for  the default fit parameters for the R-dependent
                               Galactic extinction curve from Fitzpatrick & Massa
                               (Fitzpatrick, 1999, PASP, 111, 63).

                    'LMC2' for the fit parameters are those determined for
                                  reddening the LMC2 field (inc. 30 Dor)
                                  from Misselt et al.  (1999, ApJ, 515, 128).

                    'AVGLMC' for  the fit parameters are those determined for
                                  reddening in the general Large Magellanic Cloud (LMC)
                                  field by Misselt et al.  (1999, ApJ, 515, 128).

    :Examples:
       For example::

         >>> import pyequib
         >>> wavelength=6563.0
         >>> r_v=3.1
         >>> fl=pyequib.redlaw_fm(wavelength, rv=r_v)
         >>> print('fl(6563)', fl)
            fl(6563)     -0.35054942

     :Categories:
       Interstellar Extinction

     :Dirs:
      ./
          Subroutines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.3.0

     :History:
         Based on Formulae by Fitzpatrick 1999, PASP, 11, 63
         1999PASP..111...63F, Fitzpatrick & Massa 1990,
         ApJS, 72, 163, 1990ApJS...72..163F

         Adopted from NASA IDL Library & PyAstronomy.

         30/12/2016, A. Danehkar, Revised in IDL code.

         03/10/2020, A. Danehkar, Transferred from IDL to Python.
   """
   # Tabulated inverse wavelengths in microns:
   if hasattr(wavelength, "__len__"):
     npts = len(wavelength)
     extl = np.zeros(npts)
   else:
     npts = 1
     extl = np.int32(0)
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
   for pix in range(0, npts):
   # Convert input wavelength to inverse microns
      if hasattr(wavelength, "__len__"):
         wavel=wavelength[pix]
      else:
         wavel=wavelength
      x = 10000.e+0 / wavel
      curve = x * 0.
      
      # Compute UV portion of A(lambda)/E(B-V) curve using FM fitting function and
      # R-dependent coefficients
      xcutuv = np.array([10000.0/2700.0])
      xspluv = 10000.0/np.array([2700.0,2600.0])
      
      iuv = np.where(x >= xcutuv)[0]
      n_uv = len(iuv)
      iopir = np.where(x <= xcutuv)[0]
      nopir = len(iopir)
      if (n_uv > 0):   
         xuv = np.concatenate((xspluv,x[iuv]))
      else:   
         xuv = xspluv
      
      yuv = c1 + c2 * xuv
      yuv = yuv + c3 * xuv ** 2 / ((xuv ** 2 - x0 ** 2) ** 2
                                   + (xuv * gamma1) ** 2)
      yuv = yuv + c4 * (0.5392 * ((np.maximum(xuv, 5.9)) - 5.9) ** 2
                        + 0.05644 * ((np.maximum(xuv, 5.9)) - 5.9) ** 3)
      yuv = yuv + r_v
      yspluv = yuv[0:2]                  # save spline points
      
      if (n_uv > 0):   
         curve[iuv] = yuv[2:]      # remove spline points
      
      # Compute optical portion of A(lambda)/E(B-V) curve
      # using cubic spline anchored in UV, optical, and IR
      xsplopir = np.concatenate(([0],10000.0/np.array([26500.0,12200.0,6000.0,5470.0,4670.0,4110.0])))
      ysplir   = np.array([0.0,0.26469,0.82925])*r_v/3.1 
      ysplop   = np.array((np.polyval([-4.22809e-01, 1.00270, 2.13572e-04][::-1],r_v), 
            np.polyval([-5.13540e-02, 1.00216, -7.35778e-05][::-1],r_v), 
            np.polyval([ 7.00127e-01, 1.00184, -3.32598e-05][::-1],r_v), 
            np.polyval([ 1.19456, 1.01707, -5.46959e-03, 7.97809e-04, -4.45636e-05][::-1],r_v) ))
      
      ysplopir = np.concatenate((ysplir, ysplop))
      
      if (nopir > 0):   
         tck = interpolate.splrep(np.concatenate((xsplopir,xspluv)),np.concatenate((ysplopir,yspluv)),s=0)
         if hasattr(extl, "__len__"):
            curve[iopir] = interpolate.splev(x[iopir], tck)
         else:
            curve = interpolate.splev(x, tck)
      if hasattr(extl, "__len__"):
         extl[pix] = curve
      else:
         extl=curve
   return (extl / 3.63) - 1.0

def redlaw_smc(wavelength):
   """
       This function determines the reddening law function of the line at the given wavelength
       for Small Magellanic Cloud.

    :Returns:
       type=double/array. This function returns the reddening law function value(s) for the given wavelength(s).

    :Params:
        wavelength :  in, required, type=float
                      Wavelength in Angstrom

    :Examples:
       For example::

        >>> import pyequib
        >>> wavelength=6563.0
        >>> fl=pyequib.redlaw_smc(wavelength)
        >>> print('fl(6563)', fl)
           fl(6563)     -0.22659261

    :Categories:
      Interstellar Extinction

    :Dirs:
     ./
         Subroutines

    :Author:
      Ashkbiz Danehkar

    :Copyright:
      This library is released under a GNU General Public License.

    :Version:
      0.3.0

    :History:
        Based on Prevot et al. (1984), A&A, 132, 389-392
        1984A%26A...132..389P

        Originally from IRAF STSDAS SYNPHOT redlaw.x, ebmvxfunc.x

        20/09/1994, R. A. Shaw, Initial IRAF implementation.

        04/03/1995, R. A. Shaw, Return A(lambda)/A(V) instead.

        31/08/2012, A. Danehkar, Converted to IDL code.

        03/10/2020, A. Danehkar, Transferred from IDL to Python.
   """
   # Tabulated inverse wavelengths in microns:
   xtab = np.array([0.00, 0.29, 0.45, 0.80, 1.11, 1.43, 1.82, 2.35, 2.70, 3.22, 3.34,
                    3.46, 3.60, 3.75, 3.92, 4.09, 4.28, 4.50, 4.73, 5.00, 5.24, 5.38,
                    5.52, 5.70, 5.88, 6.07, 6.27, 6.48, 6.72, 6.98, 7.23, 7.52, 7.84])
   
   # Tabulated extinction function, E(lambda-V)/E(B-V):
   extab = np.array([-3.10, -2.94, -2.72, -2.23, -1.60, -0.78, 0.00, 1.00, 1.67, 2.29,
                     2.65, 3.00, 3.15, 3.49, 3.91, 4.24, 4.53, 5.30, 5.85, 6.38, 6.76,
                     6.90, 7.17, 7.71, 8.01, 8.49, 9.06, 9.28, 9.84, 10.80, 11.51, 12.52, 13.54])
   
   if hasattr(wavelength, "__len__"):
     npts = len(wavelength)
     extl = np.zeros(npts)
   else:
     npts = 1
     extl = np.int32(0)
   for pix in range(0, npts):
      if hasattr(wavelength, "__len__"):
         wavel=wavelength[pix]
      else:
         wavel=wavelength
      if (wavel < 1000.0):
         print("redlaw_smc: Invalid wavelength")
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

def redlaw_lmc(wavelength):
   """
       This function determines the reddening law function of the line at the given wavelength
       for the Large Magellanic Cloud.

    :Returns:
       type=double/array. This function returns the reddening law function value(s) for the given wavelength(s).

    :Params:
        wavelength :  in, required, type=float
                      Wavelength in Angstrom

    :Examples:
       For example::

        >>> import pyequib
        >>> wavelength=6563.0
        >>> fl=pyequib.redlaw_lmc(wavelength)
        >>> print('fl(6563)', fl)
           fl(6563)     -0.30871187

    :Categories:
      Interstellar Extinction

    :Dirs:
     ./
         Subroutines

    :Author:
      Ashkbiz Danehkar

    :Copyright:
      This library is released under a GNU General Public License.

    :Version:
      0.3.0

    :History:
        Based on Formulae by Howarth 1983, MNRAS, 203, 301
        1983MNRAS.203..301H

        Originally from IRAF STSDAS SYNPHOT ebmvlfunc.x, redlaw.x

        18/10/1994, R. A. Shaw, Initial IRAF implementation.

        14/03/1995, R. A. Shaw, Return A(lambda)/A(V) instead.

        31/08/2012, A. Danehkar, Converted to IDL code.

        03/10/2020, A. Danehkar, Transferred from IDL to Python.
   """
   # Tabulated inverse wavelengths in microns:
   xtab = np.array([0.00, 0.29, 0.45, 0.80, 1.11, 1.43, 1.82])
   
   # Tabulated extinction function, A(lambda)/E(B-V), from Savage & Mathis:
   extab = np.array([0.00, 0.16, 0.38, 0.87, 1.50, 2.32, 3.10])
   
   if hasattr(wavelength, "__len__"):
     npts = len(wavelength)
     extl = np.zeros(npts)
   else:
     npts = 1
     extl = np.int32(0)
   for pix in range(0, npts):
      if hasattr(wavelength, "__len__"):
         wavel=wavelength[pix]
      else:
         wavel=wavelength
      if (wavel < 1000.0):   
         print('redlaw_lmc: Invalid wavelength')
      
      # Convert input wavelength to inverse microns
      x = 10000.e+0 / wavel
      
      # Infrared - optical
      if (x <= 1.82):   
         # linear interpolation of Savage & Mathis 1979
         #  val = lin_interp(extab, xtab,  x)
         #  extl[pix] = val #+ 3.1
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

def deredden_flux(wavelength, flux, m_ext, ext_law=None, rv=None, fmlaw=None):
   """
         This function dereddens absolute flux intensity based on the reddening law.

     :Returns:
        type=double. This function returns the deredden flux intensity.

     :Params:
         wavelength :  in, required, type=float/array
                       Wavelength in Angstrom
         flux       :  in, required, type=float,
                       absolute flux intensity
         m_ext      :  in, required, type=float,
                       logarithmic extinction

     :Keywords:
        ext_law  :  in, optional, type=string, default='GAL'
                    the extinction law:

                    'GAL' for Howarth Galactic.

                    'GAL2' for Savage and Mathis.

                    'CCM' for CCM galactic.

                    'JBK' for Whitford, Seaton, Kaler.

                    'FM' for Fitxpatrick.

                    'SMC' for Prevot SMC.

                    'LMC' for Howarth LMC.

        rv       :  in, optional, type=float, default=3.1
                    the optical total-to-selective extinction ratio, RV = A(V)/E(B-V).

        fmlaw    :  in, optional, type=string, default='GAL'
                    the fmlaw keyword is used only in the redlaw_fm function:

                    'GAL' for  the default fit parameters for the R-dependent
                               Galactic extinction curve from Fitzpatrick & Massa
                               (Fitzpatrick, 1999, PASP, 111, 63).

                    'LMC2' for the fit parameters are those determined for
                                  reddening the LMC2 field (inc. 30 Dor)
                                  from Misselt et al.  (1999, ApJ, 515, 128).

                    'AVGLMC' for  the fit parameters are those determined for
                                  reddening in the general Large Magellanic Cloud (LMC)
                                  field by Misselt et al.  (1999, ApJ, 515, 128).

    :Examples:
       For example::

         >>> import pyequib
         >>> wavelength=6563.0
         >>> ext_law='GAL'
         >>> r_v=3.1
         >>> m_ext=1.0
         >>> flux=1.0
         >>> flux_deredden=pyequib.deredden_flux(wavelength, flux, m_ext, ext_law=ext_law, rv=r_v) # deredden absolute flux intensity
         >>> print('dereddened flux(6563):', flux_deredden)
            dereddened flux(6563):       4.7847785

     :Categories:
       Interstellar Extinction

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.3.0

     :History:
         31/08/2012, A. Danehkar, IDL code.

         03/10/2020, A. Danehkar, Transferred from IDL to Python.
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

def deredden_relflux(wavelength, relflux, m_ext, ext_law=None, rv=None, fmlaw=None):
   """
         This function dereddens flux intensity relative to Hb=100,  based on the reddening law.

     :Returns:
        type=double. This function returns the deredden flux intensity relative to Hb=100.

     :Params:
         wavelength :  in, required, type=float/array
                       Wavelength in Angstrom
         relflux       :  in, required, type=float,
                       flux intensity relative to Hb=100
         m_ext      :  in, required, type=float,
                       logarithmic extinction

     :Keywords:
        ext_law  :  in, optional, type=string, default='GAL'
                    the extinction law:

                    'GAL' for Howarth Galactic.

                    'GAL2' for Savage and Mathis.

                    'CCM' for CCM galactic.

                    'JBK' for Whitford, Seaton, Kaler.

                    'FM' for Fitxpatrick.

                    'SMC' for Prevot SMC.

                    'LMC' for Howarth LMC.

        rv       :  in, optional, type=float, default=3.1
                    the optical total-to-selective extinction ratio, RV = A(V)/E(B-V).

        fmlaw    :  in, optional, type=string, default='GAL'
                    the fmlaw keyword is used only in the redlaw_fm function:

                    'GAL' for  the default fit parameters for the R-dependent
                               Galactic extinction curve from Fitzpatrick & Massa
                               (Fitzpatrick, 1999, PASP, 111, 63).

                    'LMC2' for the fit parameters are those determined for
                                  reddening the LMC2 field (inc. 30 Dor)
                                  from Misselt et al.  (1999, ApJ, 515, 128).

                    'AVGLMC' for  the fit parameters are those determined for
                                  reddening in the general Large Magellanic Cloud (LMC)
                                  field by Misselt et al.  (1999, ApJ, 515, 128).

    :Examples:
       For example::

         >>> import pyequib
         >>> wavelength=6563.0
         >>> ext_law='GAL'
         >>> r_v=3.1
         >>> m_ext=1.0
         >>> flux=1.0
         >>> flux_deredden=pyequib.deredden_relflux(wavelength, flux, m_ext, ext_law=ext_law, rv=r_v) # deredden absolute flux intensity
         >>> print('dereddened relative flux(6563):', flux_deredden)
            dereddened relative flux(6563):       0.47847785

     :Categories:
       Interstellar Extinction

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.3.0

     :History:
         31/08/2012, A. Danehkar, IDL code.

         03/10/2020, A. Danehkar, Transferred from IDL to Python.
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

def lin_interp(vv, xx, xout):
   """
          This function perfoms a linear interpolation/extrapolaton.

      :Private:

      :Returns:
         type=double. This function returns the interpolated/extrapolated value.

      :Params:
          vv  :  in, required, type=float
                        VV array to interpolate

          xx  :  in, required, type=float
                        X array that correspond to x(0), x(1), ...

          xout :  in, required, type=float
                        X values at which vv should be interpolated
   """
   # Make a copy so we don't overwrite the input arguments.
   v = vv
   x = xx
   # vout=np.interp(xout, xx, vv)
   interpfunc = interpolate.interp1d(xx,vv, kind='linear')
   vout=interpfunc(xout)
   
   return vout

def lin_interp_v2(vv, xx, xout):
    """
         This function perfoms a linear interpolation/extrapolaton.

     :Private:

     :Returns:
        type=double. This function returns the interpolated/extrapolated value.

     :Params:
         vv  :  in, required, type=float
                       VV array to interpolate

         xx  :  in, required, type=float
                       X array that correspond to x(0), x(1), ...

         xout :  in, required, type=float
                       X values at which vv should be interpolated
    """
    v = vv
    x = xx
    m = len(v) ## of input pnts

    type1 = type(v[0]).__name__

    #s = min(max(value_locate(x, xout), 0), (m - 2))  # Subscript intervals.
    s = (np.abs(x - xout)).argmin()
    # Linear, not regular
    _expr = (type1)
    if _expr == 'int':
        diff = v[s + 1] - np.int32(v[s])
    elif _expr == 'long':
        diff = v[s + 1] - np.int64(v[s])
    else:
        diff = v[s + 1] - v[s]

    p = (xout - x[s]) * diff / (x[s + 1] - x[s]) + v[s]
    return p

def check_sign(a, b):
    """
    NAME:
        equib_sign
    PURPOSE:

    EXPLANATION:

    CALLING SEQUENCE:
        ret= equib_sign(A, B)

    INPUTS:
        A -     A parameter
        B -     B parameter
    RETURN:  value
    """
    if b < 0:
        return -abs(a)
    else:
        return abs(a)


def do_str2int(str1):
    """Converts the list of string of miles into a list of integers of miles"""
    try:
        integer = np.int32(str1)
        return integer
    except ValueError:
        print("here")
        print("Please try again and enter a list of integers.")
        exit()


def do_strsplit(s, delim, escapech='/'):
    ret = []
    current = []
    itr = iter(s)
    for ch in itr:
        if ch == escapech:
            try:
                # skip the next character# it has been escaped!
                current.append('')
                current.append(next(itr))
            except StopIteration:
                pass
        elif ch == delim:
            # split! (add current to the list and reset it)
            ret.append(''.join(current))
            current = []
        else:
            current.append(ch)
    ret.append(''.join(current))
    return ret


def do_strtrim(s):
    if s.endswith(" "): s = s[:-1]
    if s.startswith(" "): s = s[1:]
    return s
