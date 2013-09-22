"""
This module contains functions for plasma diagnostics and abundance analysis 
from collisionally excited lines (CELs)
"""

import numpy, os
import array, math

def calc_temp_dens(ion, levu, levl, inratio, diagtype, fixedq):
   """
    NAME:
        calc_temp_dens
    PURPOSE:
        determine electron density or temperature from given
        flux intensity ratio for specified ion with upper level(s)
        lower level(s) by solving atomic level populations and
        line emissivities in statistical equilibrium
        for a fixed electron density or temperature.
   
    EXPLANATION:
   
    CALLING SEQUENCE:
        import pyequib
        ion='sii'
        levu='1,2,1,3/'
        levl='1,5/'
        diagtype='T'
        dens = 2550.0
        niiTratio=10.753
        temp=pyequib.cel.calc_temp_dens(ion, levu, levl, niiTratio, diagtype, dens)
        print temp
   
    INPUTS:
        ion -       ion name e.g. 'sii', 'nii'
        levu -      upper level(s) e.g '1,2/', '1,2,1,3/'
        levl -      lower level(s) e.g '1,2/', '1,2,1,3/'
        inratio -   flux intensity ratio
        diagtype -  diagnostics type
                    'd' or 'D' for electron density
                    't' or 'T' for electron temperature
        fixedq -    fixed quantity
                    electron density when diagtype ='t' or 'T'
                    electron temperature when diagtype ='d' or 'D'
    RETURN:  density or temperature
                    electron density when diagtype ='d' or 'D'
                    electron temperature when diagtype ='t' or 'T'
    REVISION HISTORY:
        Converted from FORTRAN to Python code by A. Danehkar, 15/09/2013
   
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
    2009-05    R.Wesson     Converted to F90, inputs from cmd line, version 
                            written purely to do diagnostics.
   """
   
   #global atomic_data_path
   
   ndim1=numpy.int32(35)
   ndim2=numpy.int32(150)
   # NDIM1T3 should be at least 3*NDIM1
   ndim1t3=numpy.int32(105)
   # Maximum no. of Ne increments
   maxnd=numpy.int32(100)
   
   gx = numpy.int32(0)
   g=(ndim2 + 1)*[numpy.int32(0)]
   id1=(2 + 1)*[numpy.int32(0)]
   jd=(2 + 1)*[numpy.int32(0)]
   wava = numpy.zeros(ndim2 + 1)
   
   itrana = numpy.zeros((2 + 1, ndim2 + 1))
   itranb = numpy.zeros((2 + 1, ndim2 + 1))
   itranc = numpy.zeros((2 + 1, ndim2 + 1))
   loop = numpy.int32(0)
   
   n = numpy.zeros(ndim2 + 1)
   tnij = numpy.zeros((ndim2 + 1,ndim2 + 1))
   fintij = numpy.zeros((ndim2 + 1,ndim2 + 1))
   wava = numpy.zeros(ndim2 + 1)
   wavb = numpy.zeros(ndim2 + 1)
   wavc = numpy.zeros(ndim2 + 1)
   cs = numpy.zeros((ndim2 + 1,ndim2 + 1))
   qeff = numpy.zeros((ndim2 + 1,ndim2 + 1))
   qq = numpy.zeros(ndim1 + 1)
   qom =  numpy.zeros((ndim1 + 1,ndim2 + 1,ndim2 + 1))
   a = numpy.zeros((ndim2 + 1,ndim2 + 1))
   e = numpy.zeros(ndim2 + 1)
   t = numpy.zeros(ndim1 + 1)
   roott = numpy.zeros(ndim1 + 1)
   x = numpy.zeros((ndim2 + 1,ndim2 + 1))
   y = numpy.zeros(ndim2 + 1)
   x2 = numpy.zeros((ndim2 + 1,ndim2 + 1))
   xkeep = numpy.zeros((ndim2 + 1,ndim2 + 1))
   y2 = numpy.zeros(ndim2 + 1)
   ykeep = numpy.zeros(ndim2 + 1)
   hmh = numpy.zeros((ndim1 + 1,ndim1 + 1))
   d = numpy.zeros(ndim1 + 1)
   valtest = numpy.zeros(3 + 1)
   label1 = (ndim2 + 1)*['']
   
   i = numpy.int32(0)
   i1 = numpy.int32(0)
   i2 = numpy.int32(0)
   j = numpy.int32(0)
   k = numpy.int32(0)
   l = numpy.int32(0)
   kk = numpy.int32(0)
   ll = numpy.int32(0)
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
   icpr = numpy.int32(0)
   ikt = numpy.int32(0)
   ia = numpy.int32(0)
   ib = numpy.int32(0)
   ic = numpy.int32(0)
   ia1 = numpy.int32(0)
   ia2 = numpy.int32(0)
   ib1 = numpy.int32(0)
   ib2 = numpy.int32(0)
   ic1 = numpy.int32(0)
   ic2 = numpy.int32(0)
   
   tempi = numpy.float64(0)
   tinc = numpy.float64(0)
   densi = numpy.float64(0)
   dinc = numpy.float64(0)
   dens = numpy.float64(0)
   dlogd = numpy.float64(0)
   temp = numpy.float64(0)
   tlogt = numpy.float64(0)
   temp2 = numpy.float64(0)
   dd = numpy.float64(0)
   deltek = numpy.float64(0)
   expe = numpy.float64(0)
   value = numpy.float64(0)
   sumn = numpy.float64(0)
   ttt = numpy.float64(0)
   ttp = numpy.float64(0)
   ahb = numpy.float64(0)
   eji = numpy.float64(0)
   wav = numpy.float64(0)
   rlint = numpy.float64(0)
   fint = numpy.float64(0)
   suma = numpy.float64(0)
   sumb = numpy.float64(0)
   sumc = numpy.float64(0)
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
   #A=dblarr(NR,NR)
   #FACT=double(0)
   
   itrana[:,:] = 0
   itranb[:,:] = 0
   itranc[:,:] = 0
   
   levu_str = strsplit(levu, ',',escapech='/')
   levl_str = strsplit(levl, ',',escapech='/')
   
   levu_num = len(levu_str)
   levl_num = len(levl_str)
   
   levu_i = int(0)
   for i in range(1, 151):
      itrana[1][i] = equib_str2int(levu_str[levu_i])
      itrana[2][i] = equib_str2int(levu_str[levu_i + 1])
      levu_i = levu_i + 2
      if levu_i >= levu_num:   
         break
   
   levl_i = int(0)
   for i in range(1, 151):
      itranb[1][i] = equib_str2int(levl_str[levl_i])
      itranb[2][i] = equib_str2int(levl_str[levl_i + 1])
      levl_i = levl_i + 2
      if levl_i >= levl_num:   
         break#
   
   #READ(levu,*) ((ITRANA(LL,KK),LL=1,2),KK=1,150)
   #READ(levl,*) ((ITRANB(LL,KK),LL=1,2),KK=1,150)
   
   ion1 = strtrim(ion)
   #atomic_filename = atomic_data_path + '/' + ion1 + '.dat'
   modelpath=getmodelpath()
   atomic_filename = modelpath+'/atomic-data/' + ion1 + '.dat'
   fp = open(atomic_filename, 'r')
   line = fp.readline()
   tempval  = map(int, line.split())
   nlines=tempval[0]
   for i in range(1, (nlines)+(1)):
      line = fp.readline()
   # Read no. of levels (max=NDIM2) NLEV,
   line = fp.readline()
   tempval  = map(int, line.split()) 
   nlev=tempval[0]
   ntemp=tempval[1]
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
      tempval= map(float, line.split()) 
      ddtemp=tempval[0]
      t[i] = ddtemp
      t[i] = math.log10(t[i])
      roott[i] = math.sqrt(t[i])
   # If IRATS=0, what tabulated are collision strengths
   line = fp.readline()
   tempval= map(int, line.split()) 
   irats=tempval[0]
   # Else Coll. rates = tabulated values * 10 ** IRATS
   
   if (ibig == 0):   
      qx = 1.0
      while (qx != 0.e0):
         lontemp1 = numpy.int32(0)
         lontemp2 = numpy.int32(0)
         ddtemp = numpy.float64(0)
         line = fp.readline()
         tempval= map(float, line.split()) 
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
         qom[k,i,j] = qx
   
   if ((ibig == 1) or (ibig == 2)):  
		line = fp.readline()
		tempval= map(float, line.split()) 
		ntra= tempval[0]
		for in1 in range(1, (ntra)+(1)):
			#readf(lun1, i, j, qom[j,i,1:(ntemp)+1])
			line = fp.readline()
			tempval= map(float, line.split()) 
			i= tempval[0]
			j= tempval[1]
			qom[1,i,j]= tempval[2]
			#READ(1,*) I,J,(QOM(ITEMP,I,J),ITEMP=1,NTEMP)
   # Read transition probabilities
   nlev1 = nlev - 1
   if (ibig == 1):   
		line = fp.readline()
		tempval= map(float, line.split()) 
		i= tempval[0]
		j= tempval[1]
		a[j,i]= tempval[2]
		#readf(lun1, i, j, a[i,j])#,L=K+1,NLEV),K=1,NLEV1
		#READ(1,7000) ((I,J,A(J,I),L=K+1,NLEV),K=1,NLEV1)
   else:   
      for k in range(1, (nlev1)+(1)):
         kp1 = k + 1
         for l in range(kp1, (nlev)+(1)):
				line = fp.readline()
				tempval= map(float, line.split()) 
				i= int(tempval[0])
				j= int(tempval[1])
				ax= tempval[2]
				a[j,i] = ax
   # Read statistical weights, energy levels (cm-1)
   for j in range(1, (nlev)+(1)):
		line = fp.readline()
		tempval= map(float, line.split()) 
		i= int(tempval[0])
		gx= int(tempval[1])
		ex= tempval[2]
		g[i] = gx
		e[i] = ex
   
   fp.close()
   
   itranc[:,:] = 0
   
   # newbit
   # set up T and D loops depending on input.
   # Read in Te and Ne where the line
   # ratio is to be calculated
   
   #*****LOOP STARTS HERE*************************
   for loop in range(1, 10):
      if ((diagtype == 't') or (diagtype == 'T')):   
         if (loop == 1):   
            tempi = 5000.0
         else:   
            tempi = valtest[1]
         int1 = 4
         tinc = (15000.0) / ((int1 - 1) ** (loop))
         densi = fixedq
         dinc = 0
         ind = 1
         # ALLOCATE(RESULTS(3,int1))
         results = numpy.zeros((3 + 1,int1 + 1))
      else:   
         if (loop == 1):   
            densi = 0.0
         else:   
            densi = valtest[2]
         ind = 4
         dinc = (100000.0) / ((ind - 1) ** (loop))
         
         tempi = fixedq
         tinc = 0
         int1 = 1
         
         #allocate(results(3,IND))
         results =  numpy.zeros((3 + 1,ind + 1))
      if (densi <= 0):   
         densi = 1
      if (tempi < 5000):   
         tempi = 5000 # add
      # Start of Te loop
      for jt in range(1, (int1)+(1)):
         temp = tempi + (jt - 1) * tinc
         # Start of Ne loop=
         for jjd in range(1, (ind)+(1)):
            dens = densi + (jjd - 1) * dinc
            # IF(DENSI.LT.30.D0) THEN
            # DENS=10.D0**DENS
            # ENDIF
            if ((temp <= 0.e0) or (dens <= 0.e0)):   
               print 'Temp = ', temp, ', Dens = ', dens
               return 0
            dlogd = math.log10(dens)
            tlogt = math.log10(temp)
            temp2 = math.sqrt(temp)
            # Form matrices
            x[:,:] = numpy.float64(0)
            cs[:,:] = numpy.float64(0)
            qeff[:,:] = numpy.float64(0)
            tnij[:,:] = numpy.float64(0)
            y[:] = numpy.float64(0)
            
            iopt = 0
            if (ntemp == 1):   
               print 'Coll. strengths available for 1 Te only - assuming const'
            else:   
               if (ntemp == 2):   
                  print 'Coll. strengths available for 2 Te only - linear interp'
               else:   
                  (t, ntemp, iopt, ndim1, ndim1t3, hmh)=equib_splmat(t, ntemp, iopt, ndim1, ndim1t3, hmh)
                  (tlogt, t, ntemp, ndim1, hmh, d)=equib_cfd(tlogt, t, ntemp, ndim1, hmh, d)
            for i in range(2, (nlev)+(1)):
               for j in range(i, (nlev)+(1)):
               #Negative!
                  deltek = (e[i - 1] - e[j]) * 1.4388463e0
                  expe = math.exp(deltek / temp)
                  for it in range(1, (ntemp)+(1)):
                  
                     if (irats == 0.e+00):   
                        qq[it] = qom[it,i - 1,j]
                     else:   
                        #Take out the exp. depend.
                        qq[it] = qom[it,i - 1,j] / expe
                        # before interpolation
                     
                  
                  if (ntemp == 1):   
                     dd = qq[1]
                  else:   
                     
                     if (ntemp == 2):   
                        dd = qq[1] + (qq[2] - qq[1]) / (t[2] - t[1]) * (tlogt - t[1])
                     else:   
                        (tlogt, dd, t, qq, ntemp, ndim1, hmh, d)=equib_cfy(tlogt, dd, t, qq, ntemp, ndim1, hmh, d)
                  if (irats == 0.e+00):   
                     cs[i - 1,j] = dd
                  else:   
                     cs[i - 1,j] = dd * expe
                  
                  if (irats == 0.e+00):   
                     qeff[i - 1,j] = 8.63e-06 * cs[i - 1,j] * expe / (g[i - 1] * temp2)
                     qeff[j,i - 1] = 8.63e-06 * cs[i - 1,j] / (g[j] * temp2)
                  else:   
                     qeff[i - 1,j] = cs[i - 1,j] * 10. ** irats
                     # Be careful
                     qeff[j,i - 1] = g[i - 1] * qeff[i - 1,j] / (expe * g[j])
                     # G integer!
            for i in range(2, (nlev)+(1)):
               for j in range(1, (nlev)+(1)):
                  if (j != i):   
                     x[i,j] = x[i,j] + dens * qeff[j,i]
                     x[i,i] = x[i,i] - dens * qeff[i,j]
                     if (j > i):   
                        x[i,j] = x[i,j] + a[j,i]
                     else:   
                        x[i,i] = x[i,i] - a[i,j]
            for i in range(2, (nlev)+(1)):
               im1 = i - 1
               value = 0.e0 - x[i,1]
               y[im1] = value
               y2[im1] = value
               ykeep[im1] = value
               for j in range(2, (nlev)+(1)):
                  jm1 = j - 1
                  value = x[i,j]
                  x[im1,jm1] = value
                  x2[im1,jm1] = value
                  xkeep[im1,jm1] = value
            # Solve matrices for populations
            (x, y, nlev1, ndim2)=equib_luslv(x, y, nlev1, ndim2)
            for i in range(nlev, 1, -1):
               n[i] = y[i - 1]
            sumn = 1.e0
            for i in range(2, (nlev)+(1)):
               sumn = sumn + n[i]
            for i in range(2, (nlev)+(1)):
               n[i] = n[i] / sumn
            n[1] = 1.e0 / sumn
            # Output data
            ttt = temp * 1.0e-4
            ttp = ttt ** (-0.87e0)
            # Eff. recombination coef. of Hb
            ahb = 3.036e-14 * ttp
            for i in range(1, (nlev1)+(1)):
               ip1 = i + 1
               for j in range(ip1, (nlev)+(1)):
                  if (a[j,i] != 0.e0):   
                     eji = e[j] - e[i]
                     wav = 1.e8 / eji
                     rlint = a[j,i] * eji
                     rlint = rlint * n[j]
                     tnij[i,j] = rlint
                     fint = n[j] * a[j,i] * 4861.e0 / (dens * ahb * wav)
                     fintij[i,j] = fint
            # Search ITRANA, ITRANB & ITRANC for transitions & sum up
            suma = 0.e0
            sumb = 0.e0
            sumc = 0.e0
            iapr = 0
            ibpr = 0
            icpr = 0
            for ikt in range(1, (ndim2)+(1)):
				ia1 = int(itrana[1,ikt])
				ia2 = int(itrana[2,ikt])
				if ((ia1 != 0) and (ia2 != 0)):   
					suma = suma + tnij[ia1,ia2]
					iapr = iapr + 1
				ib1 = int(itranb[1,ikt])
				ib2 = int(itranb[2,ikt])
				if ((ib1 != 0) and (ib2 != 0)):   
					ibpr = ibpr + 1
					sumb = sumb + tnij[ib1,ib2]
               
				ic1 = int(itranc[1,ikt])
				ic2 = int(itranc[2,ikt])
				if ((ic1 != 0) and (ic2 != 0)):   
					icpr = icpr + 1
					sumc = sumc + fintij[ic1,ic2]
			
            frat = suma / sumb
            # sumc = 1. / sumc
            # TDRAT(1,JJD)=DENS  !are these lines necessary,
            # TDRAT(2,JJD)=FRAT  !TDRAT is now never used again?
            # write(6,*),jd,suma,sumb,sumc,dens,frat
            # WRITE(7,1017) TEMP, DENS, SUMC
            # WRITE(8,1017) TEMP, DENS, FRAT, FRAT-inratio
            if ((diagtype == 't') or (diagtype == 'T')):   
               results[1,jt] = temp
               results[2,jt] = dens
               results[3,jt] = frat - inratio
            else:   
               results[1,jjd] = temp
               results[2,jjd] = dens
               results[3,jjd] = frat - inratio #End of the Ne loop
         
         for ia in range(1, (iapr)+(1)):
            i1 = int(itrana[1,ia])
            i2 = int(itrana[2,ia])
            dee = e[i2] - e[i1]
            wava[ia] = 1.e8 / dee
         for ib in range(1, (ibpr)+(1)):
            i1 = int(itranb[1,ib])
            i2 = int(itranb[2,ib])
            dee = e[i2] - e[i1]
            wavb[ib] = 1.e8 / dee
         for ic in range(1, (icpr)+(1)):
            i1 = int(itranc[1,ic])
            i2 = int(itranc[2,ic])
            dee = e[i2] - e[i1]
            wavc[ic] = 1.e8 / dee
         # End of the Te loop
      # here, find the value in RESULTS which is closest to zero
      # sort values in results to find two lowest values
      
      if ((diagtype == 'D') or (diagtype == 'd')):   
         int1 = ind
      
      # loop through array and find out where the sign changes.
      
      #    for I=1,int1 do begin
      #      if (sign(results[3,I],results[3,1]) ne results[3,I]) then begin
      #when this condition is fulfilled, the values in the array are now a different sign to the first value in the array
      #        valtest[*] = results[*,I-1] ; return the value before the sign change so that the next loop starts at a sensible value
      #        return, 0
      #      endif
      for i in range(2, (int1)+(1)):
         test = 0
         if (equib_sign(results[3,i], results[3,1]) != results[3,i]):   
            #when this condition is fulfilled, the values in the array are now a different sign to the first value in the array
            valtest[:] = results[:,i - 1] # return the value before the sign change so that the next loop starts at a sensible value
            test = 1
            break
      
      if ((test == 0) and (loop < 9)):    #test fails if no change of sign
         #this kicks in then, and checks if it should be upper or lower limit
         if (abs(results[3,1])) < (abs(results[3,int1])):   
            valtest[:] = results[:,1]
         else:   
            if (abs(results[3,int1]) < abs(results[3,1])):   
               valtest[:] = results[:,int1 - 1]
            else:   
               print 'Valtest failed'
               return 0
      else:   
         if ((test == 0) and (loop == 9)):    #test fails if no change of sign
            # this kicks in then, and checks if it should be upper or lower limit
            if (abs(results[3,1]) < abs(results[3,int1])):   
               valtest[:] = results[:,1]
            else:   
               if (abs(results[3,int1]) < abs(results[3,1])):   
                  valtest[:] = results[:,int1]
               else:   
                  print 'Valtest failed'
                  return 0
      
      #LOOP = LOOP + 1
      # DEALLOCATE(RESULTS) ; thanks Bruce!
   # ********LOOP WOULD END HERE**********************
   
   if ((diagtype == 'D') or (diagtype == 'd')):   
      result1 = valtest[2]
      # print*,valtest(2)
   else:   
      result1 = valtest[1]
   # result=1.0d+0
   # print*,result
   return result1

def calc_abundance(ion, levels, tempi, densi, iobs):
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
        import pyequib
        ion='oiii'
        tempi=10000.0
        densi=5000.0
        levels5007='3,4/'
        iobs5007=1200.0
        Abb5007=pyequib.cel.calc_abundance(ion, levels5007, tempi, densi, iobs5007)
        print Abb5007
   
    INPUTS:
        ion -       ion name e.g. 'oii', 'oiii'
        levels -    level(s) e.g '1,2/', '1,2,1,3/'
        tempi -     electron temperature
        densi -     electron density
        iobs -      observed flux intensity
    RETURN:  ionic abundance
   
    REVISION HISTORY:
        Converted from FORTRAN to Python code by A. Danehkar, 15/09/2013
   
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
    2009-05    R.Wesson     Converted to F90. Version written only for
                            calculating ionic abundances. Takes arguments
                            from the command line.
   """
   
   #global atomic_data_path
   
   ndim1=numpy.int32(35)
   ndim2=numpy.int32(150)
   # NDIM1T3 should be at least 3*NDIM1
   ndim1t3=numpy.int32(105)
   # Maximum no. of Ne increments
   maxnd=numpy.int32(100)
   
   gx = numpy.int32(0)
   g=(ndim2 + 1)*[numpy.int32(0)]
   id1=(2 + 1)*[numpy.int32(0)]
   jd=(2 + 1)*[numpy.int32(0)]
   itrana = numpy.zeros((2 + 1, ndim2 + 1))
   itranb = numpy.zeros((2 + 1, ndim2 + 1))
   itranc = numpy.zeros((2 + 1, ndim2 + 1))
   loop = numpy.int32(0)
   
   n = numpy.zeros(ndim2 + 1)
   tdrat = numpy.zeros((2 + 1, maxnd + 1))
   tnij = numpy.zeros((ndim2 + 1,ndim2 + 1))
   fintij = numpy.zeros((ndim2 + 1,ndim2 + 1))
   wava = numpy.zeros(ndim2 + 1)
   wavb = numpy.zeros(ndim2 + 1)
   wavc = numpy.zeros(ndim2 + 1)
   cs = numpy.zeros((ndim2 + 1,ndim2 + 1))
   qeff = numpy.zeros((ndim2 + 1,ndim2 + 1))
   qq = numpy.zeros(ndim1 + 1)
   qom =  numpy.zeros((ndim1 + 1,ndim2 + 1,ndim2 + 1))
   a = numpy.zeros((ndim2 + 1,ndim2 + 1))
   e = numpy.zeros(ndim2 + 1)
   t = numpy.zeros(ndim1 + 1)
   roott = numpy.zeros(ndim1 + 1)
   x = numpy.zeros((ndim2 + 1,ndim2 + 1))
   y = numpy.zeros(ndim2 + 1)
   x2 = numpy.zeros((ndim2 + 1,ndim2 + 1))
   xkeep = numpy.zeros((ndim2 + 1,ndim2 + 1))
   y2 = numpy.zeros(ndim2 + 1)
   ykeep = numpy.zeros(ndim2 + 1)
   hmh = numpy.zeros((ndim1 + 1,ndim1 + 1))
   d = numpy.zeros(ndim1 + 1)
   valtest = numpy.zeros(3 + 1)
   
   label1 = (ndim2 + 1)*['']
   
   i = numpy.int32(0)
   i1 = numpy.int32(0)
   i2 = numpy.int32(0)
   j = numpy.int32(0)
   k = numpy.int32(0)
   l = numpy.int32(0)
   kk = numpy.int32(0)
   ll = numpy.int32(0)
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
   icpr = numpy.int32(0)
   ikt = numpy.int32(0)
   ia = numpy.int32(0)
   ib = numpy.int32(0)
   ic = numpy.int32(0)
   ia1 = numpy.int32(0)
   ia2 = numpy.int32(0)
   ib1 = numpy.int32(0)
   ib2 = numpy.int32(0)
   ic1 = numpy.int32(0)
   ic2 = numpy.int32(0)
   
   #tempi = numpy.float64(0)
   tinc = numpy.float64(0)
   #densi = numpy.float64(0)
   dinc = numpy.float64(0)
   dens = numpy.float64(0)
   dlogd = numpy.float64(0)
   temp = numpy.float64(0)
   tlogt = numpy.float64(0)
   temp2 = numpy.float64(0)
   dd = numpy.float64(0)
   deltek = numpy.float64(0)
   expe = numpy.float64(0)
   value = numpy.float64(0)
   sumn = numpy.float64(0)
   ttt = numpy.float64(0)
   ttp = numpy.float64(0)
   ahb = numpy.float64(0)
   eji = numpy.float64(0)
   wav = numpy.float64(0)
   rlint = numpy.float64(0)
   fint = numpy.float64(0)
   suma = numpy.float64(0)
   sumb = numpy.float64(0)
   sumc = numpy.float64(0)
   qx = numpy.float64(0)
   ax = numpy.float64(0)
   ex = numpy.float64(0)
   frat = numpy.float64(0)
   dee = numpy.float64(0)
   ltext = ''#
   
   abund = numpy.float64(0)
   
   i = numpy.int32(0)
   j = numpy.int32(0)
   k = numpy.int32(0)
   ip1 = numpy.int32(0)
   #A=dblarr(NR,NR)
   #FACT=double(0)
   
   # g[:] = 0
   itrana[:,:] = 0
   itranb[:,:] = 0
   itranc[:,:] = 0
   
   levels_str = strsplit(levels, ',',escapech='/')
   
   levels_num = len(levels_str)

   levels_i = int(0)
   for i in range(1, 151):
      itranc[1][i] = equib_str2int(levels_str[levels_i])
      itranc[2][i] = equib_str2int(levels_str[levels_i + 1])
      levels_i = levels_i + 2
      if levels_i >= levels_num:   
         break
   
   #read(levels,*) ((ITRANC(LL,KK),LL=1,2),KK=1,150)
   
   tinc = 0
   dinc = 0
   int1 = 1
   ind = 1
   
   ion1 = strtrim(ion)
   #atomic_filename = atomic_data_path + '/' + ion1 + '.dat'
   modelpath=getmodelpath()
   atomic_filename = modelpath+'/atomic-data/' + ion1 + '.dat'
   fp = open(atomic_filename, 'r')
   line = fp.readline()
   tempval  = map(int, line.split())
   nlines=tempval[0]
   for i in range(1, (nlines)+(1)):
      line = fp.readline()
   # Read no. of levels (max=NDIM2) NLEV,
   line = fp.readline()
   tempval  = map(int, line.split()) 
   nlev=tempval[0]
   ntemp=tempval[1]
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
      tempval= map(float, line.split()) 
      ddtemp=tempval[0]
      t[i] = ddtemp
      t[i] = math.log10(t[i])
      roott[i] = math.sqrt(t[i])
   # If IRATS=0, what tabulated are collision strengths
   line = fp.readline()
   tempval= map(int, line.split()) 
   irats=tempval[0]
   # Else Coll. rates = tabulated values * 10 ** IRATS
   
   if (ibig == 0):   
      qx = 1.0
      while (qx != 0.e0):
         lontemp1 = numpy.int32(0)
         lontemp2 = numpy.int32(0)
         ddtemp = numpy.float64(0)
         line = fp.readline()
         tempval= map(float, line.split()) 
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
         qom[k,i,j] = qx
   
   if ((ibig == 1) or (ibig == 2)):  
		line = fp.readline()
		tempval= map(float, line.split()) 
		ntra= tempval[0]
		for in1 in range(1, (ntra)+(1)):
			#readf(lun1, i, j, qom[j,i,1:(ntemp)+1])
			line = fp.readline()
			tempval= map(float, line.split()) 
			i= tempval[0]
			j= tempval[1]
			qom[1,i,j]= tempval[2]
			#READ(1,*) I,J,(QOM(ITEMP,I,J),ITEMP=1,NTEMP)
   # Read transition probabilities
   nlev1 = nlev - 1
   if (ibig == 1):   
		line = fp.readline()
		tempval= map(float, line.split()) 
		i= tempval[0]
		j= tempval[1]
		a[j,i]= tempval[2]
		#readf(lun1, i, j, a[i,j])#,L=K+1,NLEV),K=1,NLEV1
		#READ(1,7000) ((I,J,A(J,I),L=K+1,NLEV),K=1,NLEV1)
   else:   
      for k in range(1, (nlev1)+(1)):
         kp1 = k + 1
         for l in range(kp1, (nlev)+(1)):
				line = fp.readline()
				tempval= map(float, line.split()) 
				i= int(tempval[0])
				j= int(tempval[1])
				ax= tempval[2]
				a[j,i] = ax
   # Read statistical weights, energy levels (cm-1)
   for j in range(1, (nlev)+(1)):
		line = fp.readline()
		tempval= map(float, line.split()) 
		i= int(tempval[0])
		gx= int(tempval[1])
		ex= tempval[2]
		g[i] = gx
		e[i] = ex
   
   fp.close()
   
   # Get levels for ratio
   # 150 large enough
   
   # Read in Te and Ne where the line
   # ratio is to be calculated
   
   # Start of Te loop
   for jt in range(1, (int1)+(1)):
         temp = tempi + (jt - 1) * tinc
         # Start of Ne loop=
         for jjd in range(1, (ind)+(1)):
            dens = densi + (jjd - 1) * dinc
            # IF(DENSI.LT.30.D0) THEN
            # DENS=10.D0**DENS
            # ENDIF
            if ((temp <= 0.e0) or (dens <= 0.e0)):   
               print 'Temp = ', temp, ', Dens = ', dens
               return 0
            dlogd = math.log10(dens)
            tlogt = math.log10(temp)
            temp2 = math.sqrt(temp)
            # Form matrices
            x[:,:] = numpy.float64(0)
            cs[:,:] = numpy.float64(0)
            qeff[:,:] = numpy.float64(0)
            tnij[:,:] = numpy.float64(0)
            y[:] = numpy.float64(0)
            
            iopt = 0
            if (ntemp == 1):   
               print 'Coll. strengths available for 1 Te only - assuming const'
            else:   
               if (ntemp == 2):   
                  print 'Coll. strengths available for 2 Te only - linear interp'
               else:   
                  (t, ntemp, iopt, ndim1, ndim1t3, hmh)=equib_splmat(t, ntemp, iopt, ndim1, ndim1t3, hmh)
                  (tlogt, t, ntemp, ndim1, hmh, d)=equib_cfd(tlogt, t, ntemp, ndim1, hmh, d)
            for i in range(2, (nlev)+(1)):
               for j in range(i, (nlev)+(1)):
               #Negative!
                  deltek = (e[i - 1] - e[j]) * 1.4388463e0
                  expe = math.exp(deltek / temp)
                  for it in range(1, (ntemp)+(1)):
                  
                     if (irats == 0.e+00):   
                        qq[it] = qom[it,i - 1,j]
                     else:   
                        #Take out the exp. depend.
                        qq[it] = qom[it,i - 1,j] / expe
                        # before interpolation
                     
                  
                  if (ntemp == 1):   
                     dd = qq[1]
                  else:   
                     
                     if (ntemp == 2):   
                        dd = qq[1] + (qq[2] - qq[1]) / (t[2] - t[1]) * (tlogt - t[1])
                     else:   
                        (tlogt, dd, t, qq, ntemp, ndim1, hmh, d)=equib_cfy(tlogt, dd, t, qq, ntemp, ndim1, hmh, d)
                  if (irats == 0.e+00):   
                     cs[i - 1,j] = dd
                  else:   
                     cs[i - 1,j] = dd * expe
                  
                  if (irats == 0.e+00):   
                     qeff[i - 1,j] = 8.63e-06 * cs[i - 1,j] * expe / (g[i - 1] * temp2)
                     qeff[j,i - 1] = 8.63e-06 * cs[i - 1,j] / (g[j] * temp2)
                  else:   
                     qeff[i - 1,j] = cs[i - 1,j] * 10. ** irats
                     # Be careful
                     qeff[j,i - 1] = g[i - 1] * qeff[i - 1,j] / (expe * g[j])
                     # G integer!
            for i in range(2, (nlev)+(1)):
               for j in range(1, (nlev)+(1)):
                  if (j != i):   
                     x[i,j] = x[i,j] + dens * qeff[j,i]
                     x[i,i] = x[i,i] - dens * qeff[i,j]
                     if (j > i):   
                        x[i,j] = x[i,j] + a[j,i]
                     else:   
                        x[i,i] = x[i,i] - a[i,j]
            for i in range(2, (nlev)+(1)):
               im1 = i - 1
               value = 0.e0 - x[i,1]
               y[im1] = value
               y2[im1] = value
               ykeep[im1] = value
               for j in range(2, (nlev)+(1)):
                  jm1 = j - 1
                  value = x[i,j]
                  x[im1,jm1] = value
                  x2[im1,jm1] = value
                  xkeep[im1,jm1] = value
            # Solve matrices for populations
            (x, y, nlev1, ndim2)=equib_luslv(x, y, nlev1, ndim2)
            for i in range(nlev, 1, -1):
               n[i] = y[i - 1]
            sumn = 1.e0
            for i in range(2, (nlev)+(1)):
               sumn = sumn + n[i]
            for i in range(2, (nlev)+(1)):
               n[i] = n[i] / sumn
            n[1] = 1.e0 / sumn
            # Output data
            ttt = temp * 1.0e-4
            ttp = ttt ** (-0.87e0)
            # Eff. recombination coef. of Hb
            ahb = 3.036e-14 * ttp
            for i in range(1, (nlev1)+(1)):
               ip1 = i + 1
               for j in range(ip1, (nlev)+(1)):
                  if (a[j,i] != 0.e0):   
                     eji = e[j] - e[i]
                     wav = 1.e8 / eji
                     rlint = a[j,i] * eji
                     rlint = rlint * n[j]
                     tnij[i,j] = rlint
                     fint = n[j] * a[j,i] * 4861.e0 / (dens * ahb * wav)
                     fintij[i,j] = fint
            # Search ITRANA, ITRANB & ITRANC for transitions & sum up
            suma = 0.e0
            sumb = 0.e0
            sumc = 0.e0
            iapr = 0
            ibpr = 0
            icpr = 0
            for ikt in range(1, (ndim2)+(1)):
				ia1 = int(itrana[1,ikt])
				ia2 = int(itrana[2,ikt])
				if ((ia1 != 0) and (ia2 != 0)):   
					suma = suma + tnij[ia1,ia2]
					iapr = iapr + 1
				ib1 = int(itranb[1,ikt])
				ib2 = int(itranb[2,ikt])
				if ((ib1 != 0) and (ib2 != 0)):   
					ibpr = ibpr + 1
					sumb = sumb + tnij[ib1,ib2]
               
				ic1 = int(itranc[1,ikt])
				ic2 = int(itranc[2,ikt])
				if ((ic1 != 0) and (ic2 != 0)):   
					icpr = icpr + 1
					sumc = sumc + fintij[ic1,ic2]
			
            # frat = suma / sumb
            sumc = 1. / sumc
            tdrat[1,jjd]=dens  # are these lines necessary,
            tdrat[2,jjd]=frat  # TDRAT is now never used again?
            abund = sumc*iobs/100.0
            # End of the Ne loop
         
   for ia in range(1, (iapr)+(1)):
		i1 = int(itrana[1,ia])
		i2 = int(itrana[2,ia])
		dee = e[i2] - e[i1]
		wava[ia] = 1.e8 / dee
   for ib in range(1, (ibpr)+(1)):
		i1 = int(itranb[1,ib])
		i2 = int(itranb[2,ib])
		dee = e[i2] - e[i1]
		wavb[ib] = 1.e8 / dee
   for ic in range(1, (icpr)+(1)):
		i1 = int(itranc[1,ic])
		i2 = int(itranc[2,ic])
		dee = e[i2] - e[i1]
		wavc[ic] = 1.e8 / dee
   # End of the Te loop
   # here, find the value in RESULTS which is closest to zero
   # sort values in results to find two lowest values
   
   return abund
   
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

def equib_reslv(a, b, n, nr):
   """
   NAME:
       equib_reslv
   PURPOSE:
       Resolve A with B
   EXPLANATION:
  
   CALLING SEQUENCE:
       equib_reslv, A, B, N, NR
  
   INPUTS:
       A -     A parameter
       B -     B parameter
       N -     N parameter
       NR -    NR parameter
   REVISION HISTORY:
       Converted from FORTRAN EQUIB to Python, 15/09/2013
   """
   #N= long(0)
   #NR= long(0)
   
   nm1 = numpy.int32(0)
   i = numpy.int32(0)
   j = numpy.int32(0)
   k = numpy.int32(0)
   l = numpy.int32(0)
   ip1 = numpy.int32(0)
   #A=dblarr(NR+1,NR+1)
   #B=dblarr(NR+1)
   if (n == 1):   
      b[n] = b[n] / a[n,n]
      return (a, b, n, nr)
   nm1 = n - 1
   for i in range(1, (nm1)+(1)):
      ip1 = i + 1
      for j in range(ip1, (n)+(1)):
         b[j] = b[j] - b[i] * a[j,i] / a[i,i]
   b[n] = b[n] / a[n,n]
   for i in range(1, (nm1)+(1)):
      k = n - i
      l = k + 1
      for j in range(l, (n)+(1)):
         b[k] = b[k] - b[j] * a[k,j]
      b[k] = b[k] / a[k,k]
   
   return (a, b, n, nr)

def equib_splmat(xx, npt, iopt, ndim, ndimt3, hmh):
   """
   NAME:
       equib_splmat
   PURPOSE:
  
   EXPLANATION:
  
   CALLING SEQUENCE:
       equib_splmat, XX, NPT, IOPT, NDIM, NDIMT3, HMH
  
   INPUTS:
       XX -     XX parameter
       NPT -    NPT parameter
       IOPT -   IOPT parameter
       NDIM -   NDIM parameter
       NDIMT3 - NDIMT3 parameter
       HMH -    HMH parameter
   REVISION HISTORY:
       Converted from FORTRAN EQUIB to Python, 15/09/2013
   """
   #NDIM= long(0)
   #NDIMT3= long(0)
   #NPT= long(0)
   #IOPT= long(0)
   #NPM= long(0)
   
   nelem = numpy.int32(0)
   #XX=dblarr(NDIM)
   gh = numpy.zeros(ndimt3 + 1)
   y = numpy.zeros(ndim + 1)
   # HMH=dblarr(NDIM+1,NDIM+1)
   npm = npt - 2
   (gh, xx, npt, iopt, ndim, ndimt3)=equib_ghgen(gh, xx, npt, iopt, ndim, ndimt3)
   nelem = 3 * npm - 2
   (gh, npm, ndimt3)=equib_elu(gh, npm, ndimt3)
   (xx, gh, y, npt, iopt, ndim, ndimt3, hmh)=equib_hgen(xx, gh, y, npt, iopt, ndim, ndimt3, hmh)
   
   return (xx, npt, iopt, ndim, ndimt3, hmh)

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
        print "here"
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
