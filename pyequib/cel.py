"""
This module contains functions for plasma diagnostics and abundance analysis 
from collisionally excited lines (CELs)
"""

# A. Danehkar
#
# Version 0.1, 15/08/2016
# First Release
#

import numpy, os
import array, math
from scipy import interpolate

def calc_emissivity(temperature=None, density=None, ion=None, levels=None):
   """
    NAME:
        calc_emissivity
    PURPOSE:
        calculate line emissivities for specified ion with level(s) by
        solving atomic level populations and in statistical equilibrium
        for given electron density and temperature.
   
    EXPLANATION:
   
    CALLING SEQUENCE:
        path='proEQUIB/atomic-data/'
        set_atomic_data_path, path
        
        ion='oiii'
        temperature=double(10000.0)
        density=double(5000.0)
        levels5007='3,4/'
        
        emiss5007=calc_emissivity(temperature=temperature, density=density, 
                                  ion=ion, levels=levels5007)
        print, emiss5007
   
    INPUTS:
        temperature   -     electron temperature
        density       -     electron density
        atomic_levels -     level(s) e.g '1,2/', '1,2,1,3/'
        elj_data      -     energy levels (Ej) data
        omij_data     -     collision strengths (omega_ij) data
        aij_data      -     transition probabilities (Aij) data
   
    RETURN:  ionic abundance
   
    REVISION HISTORY:
        Converted from FORTRAN to IDL code by A. Danehkar, 15/09/2013
        Replaced str2int with strnumber, A. Danehkar, 20/10/2016
        Replaced CFY, SPLMAT, and CFD with
             IDL function INTERPOL( /SPLINE), A. Danehkar, 20/10/2016
        Replaced LUSLV with IDL LAPACK function
                          LA_LINEAR_EQUATION, A. Danehkar, 20/10/2016
        Replaced LA_LINEAR_EQUATION (not work in GDL)
              with IDL function LUDC & LUSOL, A. Danehkar, 15/11/2016
        Replaced INTERPOL (not accurate) with
                       SPL_INIT & SPL_INTERP, A. Danehkar, 19/11/2016
        Made a new function calc_populations() for solving atomic
          level populations and separated it from
          calc_abundance(), calc_density() and calc_temperature(), A. Danehkar, 20/11/2016
        Made a new function calc_emissivity() for calculating
                         line emissivities and separated it
                         from calc_abundance(), A. Danehkar, 21/11/2016
        Integration with AtomNeb, now uses atomic data input elj_data,
                         omij_data, aij_data, A. Danehkar, 10/03/2017
        Cleaned the function, and remove unused varibales
              new function calc_emissivity(), A. Danehkar, 12/06/2017
   
    FORTRAN EQUIB HISTORY (F77/F90):
    1981-05-03 I.D.Howarth  Version 1
    1981-05-05 I.D.Howarth  Minibug fixed!
    1981-05-07 I.D.Howarth  Now takes collision rates or strengths
    1981-08-03 S.Adams      Interpolates collision strengths
    1981-08-07 S.Adams      Input method changed
    1984-11-19 R.E.S.Clegg  SA files entombed in scratch disk. Logical
                            filenames given to SA's data files.
    1995-08    D.P.Ruffle   Changed input file format. Increased matrices.
    1996-02    X.W.Liu      Tidy up. SUBROUTINES SPLMAT, HGEN, CFY and CFD
                            modified such that matrix sizes (i.e. maximum
                            of Te and maximum no of levels) can now be cha
                            by modifying the parameters NDIM1, NDIM2 and N
                            in the Main program. EASY!
                            Now takes collision rates as well.
                            All variables are declared explicitly
                            Generate two extra files (ionpop.lis and ionra
                            of plain stream format for plotting
    1996-06    C.J.Pritchet Changed input data format for cases IBIG=1,2.
                            Fixed readin bug for IBIG=2 case.
                            Now reads reformatted upsilons (easier to see
                            and the 0 0 0 data end is excluded for these c
                            The A values have a different format for IBIG=
    2006       B.Ercolano   Converted to F90
   """
   global atomic_data_path
   
   h_planck = 6.62606957e-27 # erg s
   c_speed = 2.99792458e10 # cm/s 
   
   gx = numpy.int32(0)
   id1=(2 + 1)*[numpy.int32(0)]
   jd=(2 + 1)*[numpy.int32(0)]
   
   i = numpy.int32(0)
   j = numpy.int32(0)
   k = numpy.int32(0)
   l = numpy.int32(0)
   nlines = numpy.int32(0)
   level_num = numpy.int32(0)
   temp_num = numpy.int32(0)
   ibig = numpy.int32(0)
   irats = numpy.int32(0)
   ntra = numpy.int32(0)
   itemp = numpy.int32(0)
   in1 = numpy.int32(0)
   kp1 = numpy.int32(0)
   it = numpy.int32(0)
   ip1 = numpy.int32(0)
   ikt = numpy.int32(0)
   ic = numpy.int32(0)
   ic1 = numpy.int32(0)
   ic2 = numpy.int32(0)
   
   #dens = numpy.float64(0)
   #temp = numpy.float64(0)
   
   eji = numpy.float64(0)
   wav = numpy.float64(0)
   
   qx = numpy.float64(0)
   ax = numpy.float64(0)
   ex = numpy.float64(0)
   ltext = ''#
   
   abund =  numpy.float64(0)
   
   i = numpy.int32(0)
   j = numpy.int32(0)
   k = numpy.int32(0)
   ip1 = numpy.int32(0)
   
   levels_str = strsplit(levels, ',',escapech='/')
   
   levels_num = len(levels_str)
   levels_num = int(levels_num/2)
   itranc = numpy.zeros((2 + 1, levels_num + 1))
   itranc[:,:] = 0
   levels_i = int(0)
   for i in range(1, levels_num+1):
      itranc[1][i] = equib_str2int(levels_str[levels_i])
      itranc[2][i] = equib_str2int(levels_str[levels_i + 1])
      levels_i = levels_i + 2
      if levels_i >= levels_num:   
         break
   ion1 = strtrim(ion)
   #atomic_filename = atomic_data_path + '/' + ion1 + '.dat'
   modelpath=getmodelpath()
   atomic_filename = modelpath+'/atomic-data/' + ion1 + '.dat'
   fp = open(atomic_filename, 'r')
   line = fp.readline()
   tempval  = list(map(int, line.split()))
   nlines=tempval[0]
   for i in range(1, (nlines)+(1)):
      line = fp.readline()
   # Read no. of levels (max=NDIM2) level_num,
   line = fp.readline()
   tempval  = list(map(int, line.split()))
   level_num=tempval[0]
   temp_num=tempval[1]
   
   glj = numpy.zeros(level_num+1)
   
   nlj = numpy.zeros(level_num)
   omij = numpy.zeros((temp_num+1, level_num+1, level_num+1))
   aij = numpy.zeros((level_num + 1, level_num + 1))
   elj = numpy.zeros(level_num+1)
   temp_list = numpy.zeros(temp_num+1)
   
   label1 = (level_num + 1)*['']
   
   glj[:] = 0
   # no. of Te (max=NDIM1) temp_num and the
   for i in range(1, (level_num)+(1)):
   # input format (cf Readme)
      ltext = fp.readline()
      label1[i] = ltext
   # be
   ibig = 0
   # Read in Te's where coll. strengths are tabulated
   for i in range(1, (temp_num)+(1)):
      ddtemp = numpy.float64(0)
      line = fp.readline()
      tempval= list(map(float, line.split()))
      ddtemp=tempval[0]
      temp_list[i] = ddtemp
      temp_list[i] = math.log10(temp_list[i])
   # If IRATS=0, what tabulated are collision strengths
   line = fp.readline()
   tempval= list(map(int, line.split())) 
   irats=tempval[0]
   # Else Coll. rates = tabulated values * 10 ** IRATS
   
   if (ibig == 0):   
      qx = 1.0
   while (qx != 0.e0):
      lontemp1 = numpy.int32(0)
      lontemp2 = numpy.int32(0)
      ddtemp = numpy.float64(0)
      line = fp.readline()
      tempval= list(map(float, line.split()))
      lontemp1=tempval[0]
      lontemp2=tempval[1]
      ddtemp=tempval[2]
      id1[2] = lontemp1
      jd[2] = lontemp2
      qx = ddtemp
      if qx == 0:   
         break
      if (id1[2] == 0):   
         id1[2] = id1[1]
         k = int(k + 1)
      else:   
         id1[1] = id1[2]
         k = int(1)
      if (jd[2] == 0):   
         jd[2] = jd[1]
      else:   
         jd[1] = jd[2]
      i = int(id1[2])
      j = int(jd[2])
      omij[k,i,j] = qx
   
   if ((ibig == 1) or (ibig == 2)):  
       line = fp.readline()
       tempval= list(map(float, line.split()))
       ntra= tempval[0]
       for in1 in range(1, (ntra)+(1)):
           #readf(lun1, i, j, qom[j,i,1:(temp_num)+1])
           line = fp.readline()
           tempval= list(map(float, line.split()))
           i= tempval[0]
           j= tempval[1]
           omij[1,i,j]= tempval[2]
           #READ(1,*) I,J,(QOM(ITEMP,I,J),ITEMP=1,temp_num)
   # Read transition probabilities
   level_num1 = level_num - 1
   if (ibig == 1):   
       line = fp.readline()
       tempval= list(map(float, line.split()))
       i= tempval[0]
       j= tempval[1]
       aij[j,i]= tempval[2]
       #readf(lun1, i, j, a[i,j])#,L=K+1,level_num),K=1,level_num1
       #READ(1,7000) ((I,J,A(J,I),L=K+1,level_num),K=1,level_num1)
   else:   
     for k in range(1, (level_num1)+(1)):
        kp1 = k + 1
        for l in range(kp1, (level_num)+(1)):
               line = fp.readline()
               tempval= list(map(float, line.split()))
               i= int(tempval[0])
               j= int(tempval[1])
               ax= tempval[2]
               aij[j,i] = ax
   # Read statistical weights, energy levels (cm-1)
   for j in range(1, (level_num)+(1)):
       line = fp.readline()
       tempval= list(map(float, line.split()))
       i= int(tempval[0])
       gx= int(tempval[1])
       ex= tempval[2]
       glj[i] = gx
       elj[i] = ex
   fp.close()
   
   if ((temperature <= 0.e0) or (density <= 0.e0)):  
	   print('temperature = ', temperature, ', density = ', density) 
	   return 0
   
   nlj = calc_populations(temperature=temperature, density=density, temp_list=temp_list, omij=omij, aij=aij, elj=elj, glj=glj, level_num=level_num, temp_num=temp_num, irats=irats)
   emissivity_all = numpy.float64(0)
   for ikt in range(1, (levels_num)+(1)):
      i = int(itranc[1,ikt])
      j = int(itranc[2,ikt])
      emissivity_line = numpy.float64(0)
      if (aij[j,i] != 0.e0):   
         eji = elj[j] - elj[i]
         wav = 1.e8 / eji
         emissivity_line = nlj[j] * aij[j,i] * h_planck * c_speed * 1.e8 / (wav * density)
         emissivity_all = emissivity_all + emissivity_line
   return emissivity_all

def calc_populations(temperature=None, density=None, temp_list=None, omij=None, aij=None, elj=None, glj=None, level_num=None, temp_num=None, irats=None):
   """
    NAME:
        calc_populations
    PURPOSE:
        solve atomic level populations in statistical equilibrium
        for given electron temperature and density.
   
    EXPLANATION:
   
    CALLING SEQUENCE:
        Nlj=calc_populations(temperature=temperature, density=density, temp_list=temp_list, 
                             omij=omij, aij=aij, elj=elj, glj=glj, level_num=level_num, 
                             temp_num=temp_num, irats=irats)
   
    INPUTS:
        temperature -     electron temperature
        density -     electron density
        temp_list -   temperature intervals (array)
        Omij - Collision Strengths (Omega_ij)
        Aij - Transition Probabilities (A_ij)
        Elj - Energy Levels (E_j)
        Glj - Ground Levels (G_j)
        level_num -Number of levels
        temp_num - Number of temperature intervals
        IRATS - Else Coll. rates = tabulated values * 10 ** IRATS
    RETURN:  N_j (array): atomic level populations
   
    REVISION HISTORY:
        Converted from FORTRAN to IDL code by A. Danehkar, 15/09/2013
        Replaced str2int with strnumber, A. Danehkar, 20/10/2016
        Replaced CFY, SPLMAT, and CFD with
             IDL function INTERPOL( /SPLINE), A. Danehkar, 20/10/2016
        Replaced LUSLV with IDL LAPACK function
                          LA_LINEAR_EQUATION, A. Danehkar, 20/10/2016
        Replaced LA_LINEAR_EQUATION (not work in GDL) with IDL function
                                LUDC & LUSOL, A. Danehkar, 15/11/2016
        Replaced INTERPOL (not accurate) with
                       SPL_INIT & SPL_INTERP, A. Danehkar, 19/11/2016
        Make a new function calc_populations() and separated from
          calc_abundance(), calc_density() and calc_temperature(), A. Danehkar, 20/11/2016
   
    FORTRAN EQUIB HISTORY (F77/F90):
    1981-05-03 I.D.Howarth  Version 1
    1981-05-05 I.D.Howarth  Minibug fixed!
    1981-05-07 I.D.Howarth  Now takes collision rates or strengths
    1981-08-03 S.Adams      Interpolates collision strengths
    1981-08-07 S.Adams      Input method changed
    1984-11-19 R.E.S.Clegg  SA files entombed in scratch disk. Logical
                            filenames given to SA's data files.
    1995-08    D.P.Ruffle   Changed input file format. Increased matrices.
    1996-02    X.W.Liu      Tidy up. SUBROUTINES SPLMAT, HGEN, CFY and CFD
                            modified such that matrix sizes (i.e. maximum
                            of Te and maximum no of levels) can now be cha
                            by modifying the parameters NDIM1, NDIM2 and N
                            in the Main program. EASY!
                            Now takes collision rates as well.
                            All variables are declared explicitly
                            Generate two extra files (ionpop.lis and ionra
                            of plain stream format for plotting
    1996-06    C.J.Pritchet Changed input data format for cases IBIG=1,2.
                            Fixed readin bug for IBIG=2 case.
                            Now reads reformatted upsilons (easier to see
                            and the 0 0 0 data end is excluded for these c
                            The A values have a different format for IBIG=
    2006       B.Ercolano   Converted to F90
   """
   
   dd = numpy.float64(0)
   i = numpy.int32(0)
   j = numpy.int32(0)
   k = numpy.int32(0)
   im1 = numpy.int32(0)
   jm1 = numpy.int32(0)
   deltek = numpy.float64(0)
   expe = numpy.float64(0)
   pop_sum = numpy.float64(0)
   value = numpy.float64(0)
   
   cs =  numpy.zeros((level_num + 1, level_num + 1))
   qq = numpy.zeros(temp_num + 1)
   qeff = numpy.zeros((level_num + 1, level_num + 1))
   x = numpy.zeros((level_num + 1, level_num + 1))
   y = numpy.zeros(level_num + 1)
   nlj = numpy.zeros(level_num+1)
   
   x[:,:] = 0
   cs[:,:] = 0
   qeff[:,:] = 0
   y[:] = 0
   
   level_num1 = level_num - 1
   
   tlogt = math.log10(temperature)
   temp2 = math.sqrt(temperature)
   
   #IOPT=0
   if (temp_num == 1):   
      print('Coll. strengths available for 1 Te only - assuming const')
   else:   
      if (temp_num == 2):   
         print('Coll. strengths available for 2 Te only - linear interp')
   
   for i in range(2, (level_num)+(1)):
      for j in range(i, (level_num)+(1)):
      #Negative!
         deltek = (elj[i - 1] - elj[j]) * 1.4388463e0
         expe = math.exp(deltek / temperature)
         for it in range(1, (temp_num)+(1)):
            if (irats == 0.e+00):   
               qq[it] = omij[it,i - 1,j]
            else:   
               #Take out the exp. depend.
               qq[it] = omij[it,i - 1,j] / expe
               # before interpolation
         if (temp_num == 1):   
            dd = qq[1]
         else:   
            if (temp_num == 2):   
               dd = qq[1] + (qq[2] - qq[1]) / (temp_list[2] - temp_list[1]) * (tlogt - temp_list[1])
            else:   
               interpfunc = interpolate.interp1d(temp_list[1:temp_num+1],qq[1:temp_num+1], kind='cubic')
               dd=interpfunc(tlogt)
         if (irats == 0.e+00):   
            cs[i - 1,j] = dd
         else:   
            cs[i - 1,j] = dd * expe
         if (irats == 0.e+00):   
            qeff[i - 1,j] = 8.63e-06 * cs[i - 1,j] * expe / (glj[i - 1] * temp2)
            qeff[j,i - 1] = 8.63e-06 * cs[i - 1,j] / (glj[j] * temp2)
         else:   
            qeff[i - 1,j] = cs[i - 1,j] * 10. ** irats
            # Be careful
            qeff[j,i - 1] = glj[i - 1] * qeff[i - 1,j] / (expe * glj[j])
            # G integer!
   for i in range(2, (level_num)+(1)):
      for j in range(1, (level_num)+(1)):
         if (j != i):   
            x[i,j] = x[i,j] + density * qeff[j,i]
            x[i,i] = x[i,i] - density * qeff[i,j]
            if (j > i):   
               x[i,j] = x[i,j] + aij[j,i]
            else:   
               x[i,i] = x[i,i] - aij[i,j]
   for i in range(2, (level_num)+(1)):
      im1 = i - 1
      value = 0.e0 - x[i,1]
      y[im1] = value
      for j in range(2, (level_num)+(1)):
         jm1 = j - 1
         value = x[i,j]
         x[im1,jm1] = value
   # Solve matrices for populations
   # YY=la_linear_equation(transpose(X[1:level_num1,1:level_num1]), Y[1:level_num1]); not work in GDL
   yy = numpy.linalg.solve(x[1:level_num1+1,1:level_num1+1],y[1:level_num1+1])
   y[1:level_num1+1]=yy[0:level_num1]
   for i in range(level_num, 1, -1):
      nlj[i] = y[i - 1]
   pop_sum = 1.e0
   for i in range(2, (level_num)+(1)):
      pop_sum = pop_sum + nlj[i]
   for i in range(2, (level_num)+(1)):
      nlj[i] = nlj[i] / pop_sum
   nlj[1] = 1.e0 / pop_sum
   return nlj
   
def calc_abundance(temperature=None, density=None, line_flux=None, ion=None, atomic_levels=None):
   """
    NAME:
        calc_abundance
    PURPOSE:
        determine ionic abundance from observed
        flux intensity for specified ion with level(s)
        by solving atomic level populations and
        line emissivities in statistical equilibrium
        for given electron density and temperature.
   
    EXPLANATION:
   
    CALLING SEQUENCE:
        path='proEQUIB/atomic-data/'
        set_atomic_data_path, path
   
        ion='oiii'
        temperature=double(10000.0)
        density=double(5000.0)
        levels5007='3,4/'
        iobs5007=double(1200.0)
        Abb5007=double(0.0)
        Abb5007=calc_abundance(temperature=temperature, density=density, 
                              line_flux=iobs5007, ion=ion, atomic_levels=levels5007)
        print, Abb5007
   
    INPUTS:
        ion -       ion name e.g. 'oii', 'oiii'
        levels -    level(s) e.g '1,2/', '1,2,1,3/'
        tempi -     electron temperature
        densi -     electron density
        iobs -      observed flux intensity
    RETURN:  ionic abundance
   
    REVISION HISTORY:
        Converted from FORTRAN to IDL code by A. Danehkar, 15/09/2013
        Replaced str2int with strnumber,      A. Danehkar, 20/10/2016
        Replaced CFY, SPLMAT, and CFD with
             IDL function INTERPOL( /SPLINE), A. Danehkar, 20/10/2016
        Replaced LUSLV with IDL LAPACK function
                          LA_LINEAR_EQUATION, A. Danehkar, 20/10/2016
        Replaced LA_LINEAR_EQUATION (not work in GDL)
              with IDL function LUDC & LUSOL, A. Danehkar, 15/11/2016
        Replaced INTERPOL (not accurate) with
                       SPL_INIT & SPL_INTERP, A. Danehkar, 19/11/2016
        Made a new function calc_populations() for solving atomic
          level populations and separated it from
          calc_abundance(), calc_density() and calc_temperature(), A. Danehkar, 20/11/2016
        Made a new function calc_emissivity() for calculating
                         line emissivities and separated it
                         from calc_abundance(), A. Danehkar, 21/11/2016
        Cleaned the function, and remove unused varibales
               new function calc_abundance(), A. Danehkar, 12/06/2017
   
    FORTRAN EQUIB HISTORY (F77/F90):
    1981-05-03 I.D.Howarth  Version 1
    1981-05-05 I.D.Howarth  Minibug fixed!
    1981-05-07 I.D.Howarth  Now takes collision rates or strengths
    1981-08-03 S.Adams      Interpolates collision strengths
    1981-08-07 S.Adams      Input method changed
    1984-11-19 R.E.S.Clegg  SA files entombed in scratch disk. Logical
                            filenames given to SA's data files.
    1995-08    D.P.Ruffle   Changed input file format. Increased matrices.
    1996-02    X.W.Liu      Tidy up. SUBROUTINES SPLMAT, HGEN, CFY and CFD
                            modified such that matrix sizes (i.e. maximum
                            of Te and maximum no of levels) can now be cha
                            by modifying the parameters NDIM1, NDIM2 and N
                            in the Main program. EASY!
                            Now takes collision rates as well.
                            All variables are declared explicitly
                            Generate two extra files (ionpop.lis and ionra
                            of plain stream format for plotting
    1996-06    C.J.Pritchet Changed input data format for cases IBIG=1,2.
                            Fixed readin bug for IBIG=2 case.
                            Now reads reformatted upsilons (easier to see
                            and the 0 0 0 data end is excluded for these c
                            The A values have a different format for IBIG=
    2006       B.Ercolano   Converted to F90
   """
   
   ahb = numpy.float64(0)
   
   h_planck = 6.62606957e-27 # erg s
   c_speed = 2.99792458e10 # cm/s 
   
   if ((temperature <= 0.e0) or (density <= 0.e0)):   
      print('Temperature = ', temperature, ', Density = ', density)
      return 0
   
   # Eff. recombination coef. of Hb
   t4 = temperature * 1.0e-4
   ahb = 3.036e-14 * t4 ** (-0.87e0) # Brocklehurt 1971; Aller (1984), Physics of Thermal Gaseous Nebulae, p. 76
   wavhb = 4861.33e0 #4861.D0
   emissivity_hbeta = ahb * h_planck * c_speed * 1.e8 / wavhb # N(H+) * N(e-) (erg/s) 
   # emissivity_Hbeta=1.387D-25*T4^(-0.983D0)* 10.D0^(-0.0424D0/T4) ;  Brocklehurst (1971); Aller (1984)
   
   emissivity_all = numpy.float64(0)
   emissivity_all = calc_emissivity(temperature=temperature, density=density, ion=ion, levels=atomic_levels)
   
   abund = (emissivity_hbeta / emissivity_all) * (line_flux / 100.0)
   return abund

def calc_density(line_flux_ratio=None, temperature=None, ion=None, upper_levels=None, lower_levels=None):
   """
    NAME:
        calc_density
    PURPOSE:
        determine electron density from given
        flux intensity ratio for specified ion with upper level(s)
        lower level(s) by solving atomic level populations and
        line emissivities in statistical equilibrium
        for given electron temperature.
    
    EXPLANATION:
   
    CALLING SEQUENCE:
        import pyequib
        ion='sii'
        upper_levels='1,2/'
        lower_levels='1,3/'
        temperature = 7000.0
        siiNratio=1.50
        density=calc_density(line_flux_ratio=siiNratio, density=density, $
                                     ion = ion, upper_levels=upper_levels, lower_levels=lower_levels)
        print(density)
    INPUTS:
        line_flux_ratio  -     flux intensity ratio
        temperature      -     electron temperature
        ion -       ion name e.g. 'sii', 'nii'
        upper_levels -      upper level(s) e.g '1,2/', '1,2,1,3/'
        lower_levels -      lower level(s) e.g '1,2/', '1,2,1,3/'
    RETURN:  electron density
    
    REVISION HISTORY:
        Converted from FORTRAN to Python code by A. Danehkar, 15/09/2013
        Replaced CFY, SPLMAT, and CFD with
              Python scipy.interpolate.interp1d, A. Danehkar, 20/10/2016
        Replaced LUSLV with Python numpy.linalg.solve, A. Danehkar, 20/10/2016
        Cleaning the function, and remove unused varibales
                     new function calc_density(), A. Danehkar, 12/06/2017
                           
    FORTRAN EQUIB HISTORY (F77/F90):
    1981-05-03 I.D.Howarth  Version 1
    1981-05-05 I.D.Howarth  Minibug fixed!
    1981-05-07 I.D.Howarth  Now takes collision rates or strengths
    1981-08-03 S.Adams      Interpolates collision strengths
    1981-08-07 S.Adams      Input method changed
    1984-11-19 R.E.S.Clegg  SA files entombed in scratch disk. Logical
                            filenames given to SA's data files.
    1995-08    D.P.Ruffle   Changed input file format. Increased matrices.
    1996-02    X.W.Liu      Tidy up. SUBROUTINES SPLMAT, HGEN, CFY and CFD
                            modified such that matrix sizes (i.e. maximum
                            of Te and maximum no of levels) can now be cha
                            by modifying the parameters NDIM1, NDIM2 and N
                            in the Main program. EASY!
                            Now takes collision rates as well.
                            All variables are declared explicitly
                            Generate two extra files (ionpop.lis and ionra
                            of plain stream format for plotting
    1996-06    C.J.Pritchet Changed input data format for cases IBIG=1,2.
                            Fixed readin bug for IBIG=2 case.
                            Now reads reformatted upsilons (easier to see
                            and the 0 0 0 data end is excluded for these c
                            The A values have a different format for IBIG=
    2006       B.Ercolano   Converted to F90
   """
   h_planck = 6.62606957e-27 # erg s
   c_speed = 2.99792458e10 # cm/s 
   
   gx = numpy.int32(0)
   id1=(2 + 1)*[numpy.int32(0)]
   jd=(2 + 1)*[numpy.int32(0)]
   
   iteration = numpy.int32(0)
   
   check_value = numpy.zeros(3 + 1)
   
   i = numpy.int32(0)
   i1 = numpy.int32(0)
   i2 = numpy.int32(0)
   j = numpy.int32(0)
   k = numpy.int32(0)
   l = numpy.int32(0)
   kk = numpy.int32(0)
   jt = numpy.int32(0)
   jjd = numpy.int32(0)
   ionl = numpy.int32(0)
   nlines = numpy.int32(0)
   nlev = numpy.int32(0)
   ntemp = numpy.int32(0)
   ibig = numpy.int32(0)
   irats = numpy.int32(0)
   ntra = numpy.int32(0)
   itemp = numpy.int32(0)
   in1 = numpy.int32(0)
   nlev1 = numpy.int32(0)
   kp1 = numpy.int32(0)
   int1 = numpy.int32(0)
   ind = numpy.int32(0)
   iopt = numpy.int32(0)
   it = numpy.int32(0)
   im1 = numpy.int32(0)
   jm1 = numpy.int32(0)
   ip1 = numpy.int32(0)
   iapr = numpy.int32(0)
   ibpr = numpy.int32(0)
   ikt = numpy.int32(0)
   ia = numpy.int32(0)
   ib = numpy.int32(0)
   ia1 = numpy.int32(0)
   ia2 = numpy.int32(0)
   ib1 = numpy.int32(0)
   ib2 = numpy.int32(0)
   
   tempi = numpy.float64(0)
   tinc = numpy.float64(0)
   densi = numpy.float64(0)
   dinc = numpy.float64(0)
   dens = numpy.float64(0)
   dlogd = numpy.float64(0)
   tlogt = numpy.float64(0)
   temp2 = numpy.float64(0)
   dd = numpy.float64(0)
   deltek = numpy.float64(0)
   expe = numpy.float64(0)
   value = numpy.float64(0)
   eji = numpy.float64(0)
   wav = numpy.float64(0)
   qx = numpy.float64(0)
   ax = numpy.float64(0)
   ex = numpy.float64(0)
   frat = numpy.float64(0)
   dee = numpy.float64(0)
   ltext = ''#
   
   result1 = numpy.float64(0)
   
   i = numpy.int32(0)
   j = numpy.int32(0)
   k = numpy.int32(0)
   ip1 = numpy.int32(0)
   
   upper_levels_str = strsplit(upper_levels, ',',escapech='/')
   lower_levels_str = strsplit(lower_levels, ',',escapech='/')
   
   upper_levels_num = len(upper_levels_str)
   lower_levels_num = len(lower_levels_str)
   
   itrana = numpy.zeros((2 + 1, upper_levels_num + 1))
   itranb = numpy.zeros((2 + 1, lower_levels_num + 1))
   itrana[:,:] = 0
   itranb[:,:] = 0
   upper_levels_i = int(0)
   for i in range(1, upper_levels_num + 1):
      itrana[1][i] = equib_str2int(upper_levels_str[upper_levels_i])
      itrana[2][i] = equib_str2int(upper_levels_str[upper_levels_i + 1])
      upper_levels_i = upper_levels_i + 2
      if upper_levels_i >= upper_levels_num:   
         break
   
   lower_levels_i = int(0)
   for i in range(1, lower_levels_num+1):
      itranb[1][i] = equib_str2int(lower_levels_str[lower_levels_i])
      itranb[2][i] = equib_str2int(lower_levels_str[lower_levels_i + 1])
      lower_levels_i = lower_levels_i + 2
      if lower_levels_i >= lower_levels_num:   
         break#
   
   ion1 = strtrim(ion)
   #atomic_filename = atomic_data_path + '/' + ion1 + '.dat'
   modelpath=getmodelpath()
   atomic_filename = modelpath+'/atomic-data/' + ion1 + '.dat'
   fp = open(atomic_filename, 'r')
   line = fp.readline()
   tempval  = list(map(int, line.split()))
   nlines=tempval[0]
   for i in range(1, (nlines)+(1)):
      line = fp.readline()
   # Read no. of levels (max=NDIM2) NLEV,
   line = fp.readline()
   tempval  = list(map(int, line.split()))
   nlev=tempval[0]
   ntemp=tempval[1]
   
   glj = numpy.zeros(nlev+1)
   
   nlj = numpy.zeros(nlev+1)
   wava = numpy.zeros(nlev + 1)
   wavb = numpy.zeros(nlev + 1)
   omij = numpy.zeros((ntemp+1, nlev+1, nlev+1))
   aij = numpy.zeros((nlev + 1, nlev + 1))
   elj = numpy.zeros(nlev+1)
   telist = numpy.zeros(ntemp+1)
   
   label1 = (nlev + 1)*['']
   
   label1 = (nlev + 1)*['']
   
   glj[:] = 0
   # no. of Te (max=NDIM1) NTEMP and the
   for i in range(1, (nlev)+(1)):
   # input format (cf Readme)
      ltext = fp.readline()
      label1[i] = ltext
   # be
   ibig = 0
   # Read in Te's where coll. strengths are tabulated
   for i in range(1, (ntemp)+(1)):
      ddtemp = numpy.float64(0)
      line = fp.readline()
      tempval= list(map(float, line.split())) 
      ddtemp=tempval[0]
      telist[i] = ddtemp
      telist[i] = math.log10(telist[i])
   # If IRATS=0, what tabulated are collision strengths
   line = fp.readline()
   tempval= list(map(int, line.split()))
   irats=tempval[0]
   # Else Coll. rates = tabulated values * 10 ** IRATS
   
   if (ibig == 0):   
      qx = 1.0
      while (qx != 0.e0):
         lontemp1 = numpy.int32(0)
         lontemp2 = numpy.int32(0)
         ddtemp = numpy.float64(0)
         line = fp.readline()
         tempval= list(map(float, line.split()))
         lontemp1=tempval[0]
         lontemp2=tempval[1]
         ddtemp=tempval[2]
         id1[2] = lontemp1
         jd[2] = lontemp2
         qx = ddtemp
         if qx == 0:   
            break
         if (id1[2] == 0):   
            id1[2] = id1[1]
            k = int(k + 1)
         else:   
            id1[1] = id1[2]
            k = int(1)
         if (jd[2] == 0):   
            jd[2] = jd[1]
         else:   
            jd[1] = jd[2]
         i = int(id1[2])
         j = int(jd[2])
         omij[k,i,j] = qx
   
   if ((ibig == 1) or (ibig == 2)):  
      line = fp.readline()
      tempval= list(map(float, line.split()))
      ntra= tempval[0]
      for in1 in range(1, (ntra)+(1)):
          line = fp.readline() 
          #readf(lun1, i, j, qom[j,i,1:(ntemp)+1])
          tempval= list(map(float, line.split()))
          i= tempval[0]
          j= tempval[1]
          omij[1,i,j]= tempval[2]
          #READ(1,*) I,J,(QOM(ITEMP,I,J),ITEMP=1,NTEMP)
   # Read transition probabilities
   nlev1 = nlev - 1
   if (ibig == 1):   
      line = fp.readline()
      tempval= list(map(float, line.split()))
      i= tempval[0]
      j= tempval[1]
      aij[j,i]= tempval[2]
      #readf(lun1, i, j, a[i,j])#,L=K+1,NLEV),K=1,NLEV1
      #READ(1,7000) ((I,J,A(J,I),L=K+1,NLEV),K=1,NLEV1)
   else:   
      for k in range(1, (nlev1)+(1)):
         kp1 = k + 1
         for l in range(kp1, (nlev)+(1)):
             line = fp.readline()
             tempval= list(map(float, line.split()))
             i= int(tempval[0])
             j= int(tempval[1])
             ax= tempval[2]
             aij[j,i] = ax
   # Read statistical weights, energy levels (cm-1)
   for j in range(1, (nlev)+(1)):
        line = fp.readline()
        tempval= list(map(float, line.split()))
        i= int(tempval[0])
        gx= int(tempval[1])
        ex= tempval[2]
        glj[i] = gx
        elj[i] = ex
   
   fp.close()
   
   #set temperature iterations
   #****************************
   for iteration in range(1, 10):
      if (iteration == 1):   
         densi = 0.0
      else:   
         densi = check_value[2]
      ind = 4
      dinc = (100000.0) / ((ind - 1) ** (iteration))
      tempi = temperature
      tinc = 0
      int1 = 1
      results =  numpy.zeros((3 + 1,ind + 1))
      if (densi <= 0):   
         densi = 1
      if (tempi < 5000):   
         tempi = 5000 
      # Start of density iteration=
      for jt in range(1, (int1)+(1)):
         temp = tempi + (jt - 1) * tinc
         # Start of density iteration=
         for jjd in range(1, (ind)+(1)):
            dens = densi + (jjd - 1) * dinc
            # IF(DENSI.LT.30.D0) THEN
            # DENS=10.D0**DENS
            # ENDIF
            if ((temp <= 0.e0) or (dens <= 0.e0)):   
               print('Temp = ', temp, ', Dens = ', dens)
               return 0
            nlj = calc_populations(temperature=temp, density=dens, temp_list=telist, omij=omij, aij=aij, elj=elj, glj=glj, level_num=nlev, temp_num=ntemp, irats=irats)
            # Search ITRANA & ITRANB  for transitions & sum up
            suma = 0.e0
            sumb = 0.e0
            iapr = 0
            ibpr = 0
            for ikt in range(1, (upper_levels_num)+(1)):
                i = int(itrana[1,ikt])
                j = int(itrana[2,ikt])
                if (aij[j,i] != 0.e0):  
                    eji = elj[j] - elj[i]
                    wav = 1.e8 / eji
                    suma = suma + nlj[j] * aij[j,i] * h_planck * c_speed * 1.e8 / wav
            for ikt in range(1, (lower_levels_num)+(1)):
                i = int(itranb[1,ikt])
                j = int(itranb[2,ikt])
                if (aij[j,i] != 0.e0):  
                    eji = elj[j] - elj[i]
                    wav = 1.e8 / eji
                    sumb = sumb + nlj[j] * aij[j,i] * h_planck * c_speed * 1.e8 / wav
            frat = suma / sumb
              
            results[1,jjd] = temp
            results[2,jjd] = dens
            results[3,jjd] = frat - line_flux_ratio #End of the density iteration
         
         for ia in range(1, (upper_levels_num)+(1)):
            i1 = int(itrana[1,ia])
            i2 = int(itrana[2,ia])
            if (aij[i2,i1] != 0.e0):  
               dee = elj[i2] - elj[i1]
               wava[ia] = 1.e8 / dee
         for ib in range(1, (lower_levels_num)+(1)):
            i1 = int(itranb[1,ib])
            i2 = int(itranb[2,ib])
            if (aij[i2,i1] != 0.e0):  
               dee = elj[i2] - elj[i1]
               wavb[ib] = 1.e8 / dee
         # End of the temperature iteration
      
      int1 = ind
      
      # iteration through array and find out where the sign changes.
      
      for i in range(2, (int1)+(1)):
         check = 0
         if (equib_sign(results[3,i], results[3,1]) != results[3,i]):   
            #if this condition, the values have a different sign
            check_value[:] = results[:,i - 1] # the value before the sign change returned
            check = 1
            break
      
      if ((check == 0) and (iteration < 9)):    #check if there is any change of sign,
         #and checks if it should be upper or lower limit
         if (abs(results[3,1])) < (abs(results[3,int1])):   
            check_value[:] = results[:,1]
         else:   
            if (abs(results[3,int1]) < abs(results[3,1])):   
               check_value[:] = results[:,int1 - 1]
            else:   
               print('check_value is wrong')
               return 0
      else:   
         if ((check == 0) and (iteration == 9)):    #check if no change of sign,
            # and checks if it should be upper or lower limit
            if (abs(results[3,1]) < abs(results[3,int1])):   
               check_value[:] = results[:,1]
            else:   
               if (abs(results[3,int1]) < abs(results[3,1])):   
                  check_value[:] = results[:,int1]
               else:   
                  print('check_value is wrong')
                  return 0
   # end of iterations
   # ****************************  
   result1 = check_value[2]
   return result1

def calc_temperature(line_flux_ratio=None, density=None, ion=None, upper_levels=None, lower_levels=None):
   """
    NAME:
        calc_temperature
    PURPOSE:
        determine electron temperature from given
        flux intensity ratio for specified ion with upper level(s)
        lower level(s) by solving atomic level populations and
        line emissivities in statistical equilibrium
        for given electron density.
   
    EXPLANATION:
   
    CALLING SEQUENCE:
        import pyequib
        ion='s_ii'
        upper_levels='1,2,1,3/'
        lower_levels='1,5/'
        density = 2550.0
        siiTratio=10.753
        temperature=calc_temperature(line_flux_ratio=siiTratio, density=density, $
                                     ion = ion, upper_levels=upper_levels, lower_levels=lower_levels)
        print(temperature)
   
    INPUTS:
        line_flux_ratio  -     flux intensity ratio
        density          -     electron density
        ion -       ion name e.g. 'sii', 'nii'
        upper_levels -      upper level(s) e.g '1,2/', '1,2,1,3/'
        lower_levels -      lower level(s) e.g '1,2/', '1,2,1,3/'
    RETURN:  electron temperature
    
    REVISION HISTORY:
        Converted from FORTRAN to Python code by A. Danehkar, 15/09/2013
        Replaced CFY, SPLMAT, and CFD with
              Python scipy.interpolate.interp1d, A. Danehkar, 20/10/2016
        Replaced LUSLV with Python numpy.linalg.solve, A. Danehkar, 20/10/2016
        Cleaning the function, and remove unused varibales
                     new function cal_temperature(), A. Danehkar, 12/06/2017
   
    FORTRAN EQUIB HISTORY (F77/F90):
    1981-05-03 I.D.Howarth  Version 1
    1981-05-05 I.D.Howarth  Minibug fixed!
    1981-05-07 I.D.Howarth  Now takes collision rates or strengths
    1981-08-03 S.Adams      Interpolates collision strengths
    1981-08-07 S.Adams      Input method changed
    1984-11-19 R.E.S.Clegg  SA files entombed in scratch disk. Logical
                            filenames given to SA's data files.
    1995-08    D.P.Ruffle   Changed input file format. Increased matrices.
    1996-02    X.W.Liu      Tidy up. SUBROUTINES SPLMAT, HGEN, CFY and CFD
                            modified such that matrix sizes (i.e. maximum
                            of Te and maximum no of levels) can now be cha
                            by modifying the parameters NDIM1, NDIM2 and N
                            in the Main program. EASY!
                            Now takes collision rates as well.
                            All variables are declared explicitly
                            Generate two extra files (ionpop.lis and ionra
                            of plain stream format for plotting
    1996-06    C.J.Pritchet Changed input data format for cases IBIG=1,2.
                            Fixed readin bug for IBIG=2 case.
                            Now reads reformatted upsilons (easier to see
                            and the 0 0 0 data end is excluded for these c
                            The A values have a different format for IBIG=
    2006       B.Ercolano   Converted to F90
   """
   h_planck = 6.62606957e-27 # erg s
   c_speed = 2.99792458e10 # cm/s 
   
   gx = numpy.int32(0)
   id1=(2 + 1)*[numpy.int32(0)]
   jd=(2 + 1)*[numpy.int32(0)]
   
   iteration = numpy.int32(0)
   
   check_value = numpy.zeros(3 + 1)
   
   i = numpy.int32(0)
   i1 = numpy.int32(0)
   i2 = numpy.int32(0)
   j = numpy.int32(0)
   k = numpy.int32(0)
   l = numpy.int32(0)
   kk = numpy.int32(0)
   jt = numpy.int32(0)
   jjd = numpy.int32(0)
   ionl = numpy.int32(0)
   nlines = numpy.int32(0)
   nlev = numpy.int32(0)
   ntemp = numpy.int32(0)
   ibig = numpy.int32(0)
   irats = numpy.int32(0)
   ntra = numpy.int32(0)
   itemp = numpy.int32(0)
   in1 = numpy.int32(0)
   nlev1 = numpy.int32(0)
   kp1 = numpy.int32(0)
   int1 = numpy.int32(0)
   ind = numpy.int32(0)
   iopt = numpy.int32(0)
   it = numpy.int32(0)
   im1 = numpy.int32(0)
   jm1 = numpy.int32(0)
   ip1 = numpy.int32(0)
   iapr = numpy.int32(0)
   ibpr = numpy.int32(0)
   ikt = numpy.int32(0)
   ia = numpy.int32(0)
   ib = numpy.int32(0)
   ia1 = numpy.int32(0)
   ia2 = numpy.int32(0)
   ib1 = numpy.int32(0)
   ib2 = numpy.int32(0)
   
   tempi = numpy.float64(0)
   tinc = numpy.float64(0)
   densi = numpy.float64(0)
   dinc = numpy.float64(0)
   dlogd = numpy.float64(0)
   temp = numpy.float64(0)
   tlogt = numpy.float64(0)
   temp2 = numpy.float64(0)
   dd = numpy.float64(0)
   deltek = numpy.float64(0)
   expe = numpy.float64(0)
   value = numpy.float64(0)
   eji = numpy.float64(0)
   wav = numpy.float64(0)
   qx = numpy.float64(0)
   ax = numpy.float64(0)
   ex = numpy.float64(0)
   frat = numpy.float64(0)
   dee = numpy.float64(0)
   ltext = ''#
   
   result1 = numpy.float64(0)
   
   i = numpy.int32(0)
   j = numpy.int32(0)
   k = numpy.int32(0)
   ip1 = numpy.int32(0)
   
   upper_levels_str = strsplit(upper_levels, ',',escapech='/')
   lower_levels_str = strsplit(lower_levels, ',',escapech='/')
   
   upper_levels_num = len(upper_levels_str)
   lower_levels_num = len(lower_levels_str)
   
   itrana = numpy.zeros((2 + 1, upper_levels_num + 1))
   itranb = numpy.zeros((2 + 1, lower_levels_num + 1))
   itrana[:,:] = 0
   itranb[:,:] = 0
   upper_levels_i = int(0)
   for i in range(1, upper_levels_num + 1):
      itrana[1][i] = equib_str2int(upper_levels_str[upper_levels_i])
      itrana[2][i] = equib_str2int(upper_levels_str[upper_levels_i + 1])
      upper_levels_i = upper_levels_i + 2
      if upper_levels_i >= upper_levels_num:   
         break
   
   lower_levels_i = int(0)
   for i in range(1, lower_levels_num+1):
      itranb[1][i] = equib_str2int(lower_levels_str[lower_levels_i])
      itranb[2][i] = equib_str2int(lower_levels_str[lower_levels_i + 1])
      lower_levels_i = lower_levels_i + 2
      if lower_levels_i >= lower_levels_num:   
         break#
   
   ion1 = strtrim(ion)
   #atomic_filename = atomic_data_path + '/' + ion1 + '.dat'
   modelpath=getmodelpath()
   atomic_filename = modelpath+'/atomic-data/' + ion1 + '.dat'
   fp = open(atomic_filename, 'r')
   line = fp.readline()
   tempval  = list(map(int, line.split()))
   nlines=tempval[0]
   for i in range(1, (nlines)+(1)):
      line = fp.readline()
   # Read no. of levels (max=NDIM2) NLEV,
   line = fp.readline()
   tempval  = list(map(int, line.split()))
   nlev=tempval[0]
   ntemp=tempval[1]
   
   glj = numpy.zeros(nlev+1)
   
   nlj = numpy.zeros(nlev+1)
   wava = numpy.zeros(nlev + 1)
   wavb = numpy.zeros(nlev + 1)
   omij = numpy.zeros((ntemp+1, nlev+1, nlev+1))
   aij = numpy.zeros((nlev + 1, nlev + 1))
   elj = numpy.zeros(nlev+1)
   telist = numpy.zeros(ntemp+1)
   
   label1 = (nlev + 1)*['']
   
   label1 = (nlev + 1)*['']
   
   glj[:] = 0
   # no. of Te (max=NDIM1) NTEMP and the
   for i in range(1, (nlev)+(1)):
   # input format (cf Readme)
      ltext = fp.readline()
      label1[i] = ltext
   # be
   ibig = 0
   # Read in Te's where coll. strengths are tabulated
   for i in range(1, (ntemp)+(1)):
      ddtemp = numpy.float64(0)
      line = fp.readline()
      tempval= list(map(float, line.split())) 
      ddtemp=tempval[0]
      telist[i] = ddtemp
      telist[i] = math.log10(telist[i])
   # If IRATS=0, what tabulated are collision strengths
   line = fp.readline()
   tempval= list(map(int, line.split()))
   irats=tempval[0]
   # Else Coll. rates = tabulated values * 10 ** IRATS
   
   if (ibig == 0):   
      qx = 1.0
      while (qx != 0.e0):
         lontemp1 = numpy.int32(0)
         lontemp2 = numpy.int32(0)
         ddtemp = numpy.float64(0)
         line = fp.readline()
         tempval= list(map(float, line.split()))
         lontemp1=tempval[0]
         lontemp2=tempval[1]
         ddtemp=tempval[2]
         id1[2] = lontemp1
         jd[2] = lontemp2
         qx = ddtemp
         if qx == 0:   
            break
         if (id1[2] == 0):   
            id1[2] = id1[1]
            k = int(k + 1)
         else:   
            id1[1] = id1[2]
            k = int(1)
         if (jd[2] == 0):   
            jd[2] = jd[1]
         else:   
            jd[1] = jd[2]
         i = int(id1[2])
         j = int(jd[2])
         omij[k,i,j] = qx
   
   if ((ibig == 1) or (ibig == 2)):  
      line = fp.readline()
      tempval= list(map(float, line.split()))
      ntra= tempval[0]
      for in1 in range(1, (ntra)+(1)):
          line = fp.readline() 
          #readf(lun1, i, j, qom[j,i,1:(ntemp)+1])
          tempval= list(map(float, line.split()))
          i= tempval[0]
          j= tempval[1]
          omij[1,i,j]= tempval[2]
          #READ(1,*) I,J,(QOM(ITEMP,I,J),ITEMP=1,NTEMP)
   # Read transition probabilities
   nlev1 = nlev - 1
   if (ibig == 1):   
      line = fp.readline()
      tempval= list(map(float, line.split()))
      i= tempval[0]
      j= tempval[1]
      aij[j,i]= tempval[2]
      #readf(lun1, i, j, a[i,j])#,L=K+1,NLEV),K=1,NLEV1
      #READ(1,7000) ((I,J,A(J,I),L=K+1,NLEV),K=1,NLEV1)
   else:   
      for k in range(1, (nlev1)+(1)):
         kp1 = k + 1
         for l in range(kp1, (nlev)+(1)):
             line = fp.readline()
             tempval= list(map(float, line.split()))
             i= int(tempval[0])
             j= int(tempval[1])
             ax= tempval[2]
             aij[j,i] = ax
   # Read statistical weights, energy levels (cm-1)
   for j in range(1, (nlev)+(1)):
        line = fp.readline()
        tempval= list(map(float, line.split()))
        i= int(tempval[0])
        gx= int(tempval[1])
        ex= tempval[2]
        glj[i] = gx
        elj[i] = ex
   
   fp.close()
   
   #set temperature iterations
   #****************************
   for iteration in range(1, 10):
      if (iteration == 1):   
         tempi = 5000.0
      else:   
         tempi = check_value[1]
      int1 = 4
      tinc = (15000.0) / ((int1 - 1) ** (iteration))
      densi = density
      dinc = 0
      ind = 1
      results = numpy.zeros((3 + 1,int1 + 1))
      if (densi <= 0):   
         densi = 1
      if (tempi < 5000):   
         tempi = 5000 
      # Start of density iteration=
      for jt in range(1, (int1)+(1)):
         temp = tempi + (jt - 1) * tinc
         # Start of density iteration=
         for jjd in range(1, (ind)+(1)):
            dens = densi + (jjd - 1) * dinc
            # IF(DENSI.LT.30.D0) THEN
            # DENS=10.D0**DENS
            # ENDIF
            if ((temp <= 0.e0) or (dens <= 0.e0)):   
               print('Temp = ', temp, ', Dens = ', dens)
               return 0
            nlj = calc_populations(temperature=temp, density=dens, temp_list=telist, omij=omij, aij=aij, elj=elj, glj=glj, level_num=nlev, temp_num=ntemp, irats=irats)
            # Search ITRANA & ITRANB  for transitions & sum up
            suma = 0.e0
            sumb = 0.e0
            iapr = 0
            ibpr = 0
            for ikt in range(1, (upper_levels_num)+(1)):
                i = int(itrana[1,ikt])
                j = int(itrana[2,ikt])
                if (aij[j,i] != 0.e0):  
                    eji = elj[j] - elj[i]
                    wav = 1.e8 / eji
                    suma = suma + nlj[j] * aij[j,i] * h_planck * c_speed * 1.e8 / wav
            for ikt in range(1, (lower_levels_num)+(1)):
                i = int(itranb[1,ikt])
                j = int(itranb[2,ikt])
                if (aij[j,i] != 0.e0):  
                    eji = elj[j] - elj[i]
                    wav = 1.e8 / eji
                    sumb = sumb + nlj[j] * aij[j,i] * h_planck * c_speed * 1.e8 / wav
            frat = suma / sumb
            
            results[1,jt] = temp
            results[2,jt] = dens
            results[3,jt] = frat - line_flux_ratio
         
         for ia in range(1, (upper_levels_num)+(1)):
            i1 = int(itrana[1,ia])
            i2 = int(itrana[2,ia])
            if (aij[i2,i1] != 0.e0):  
               dee = elj[i2] - elj[i1]
               wava[ia] = 1.e8 / dee
         for ib in range(1, (lower_levels_num)+(1)):
            i1 = int(itranb[1,ib])
            i2 = int(itranb[2,ib])
            if (aij[i2,i1] != 0.e0):  
               dee = elj[i2] - elj[i1]
               wavb[ib] = 1.e8 / dee
         # End of the temperature iteration
      #  iteration and detect the sign change.
      
      for i in range(2, (int1)+(1)):
         check = 0
         if (equib_sign(results[3,i], results[3,1]) != results[3,i]):   
            #if this condition, the values have a different sign
            check_value[:] = results[:,i - 1] # the value before the sign change returned
            check = 1
            break
      
      if ((check == 0) and (iteration < 9)):    #check if there is any change of sign,
         #and checks if it should be upper or lower limit
         if (abs(results[3,1])) < (abs(results[3,int1])):   
            check_value[:] = results[:,1]
         else:   
            if (abs(results[3,int1]) < abs(results[3,1])):   
               check_value[:] = results[:,int1 - 1]
            else:   
               print('check_value is wrong')
               return 0
      else:   
         if ((check == 0) and (iteration == 9)):    #check if no change of sign,
            # and checks if it should be upper or lower limit
            if (abs(results[3,1]) < abs(results[3,int1])):   
               check_value[:] = results[:,1]
            else:   
               if (abs(results[3,int1]) < abs(results[3,1])):   
                  check_value[:] = results[:,int1]
               else:   
                  print('check_value is wrong')
                  return 0
   # end of iterations
   # ****************************
   result1 = check_value[1]
   return result1

def equib_sign(a, b):
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

def equib_str2int(str1):
    """Converts the list of string of miles into a list of integers of miles"""
    try:
        integer = int(str1)
        return integer
    except ValueError:
        print("here")
        sys.exit("Please try again and enter a list of integers.")

def strsplit(s, delim, escapech='/'):
    ret = []
    current = []
    itr = iter(s)
    for ch in itr:
        if ch == escapech:
            try:
                # skip the next character; it has been escaped!
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

def strtrim(s):
    if s.endswith(" "): s = s[:-1]
    if s.startswith(" "): s = s[1:]
    return s
    
def getmodelpath(): 
    path = os.path.dirname(__file__)
    #path = path + '/'
    return path
