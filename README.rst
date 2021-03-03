======================
pyEQUIB Python Package
======================

.. image:: https://img.shields.io/pypi/v/pyequib.svg?style=flat
    :target: https://pypi.python.org/pypi/pyequib/
    :alt: PyPI Version
    
.. image:: https://travis-ci.org/equib/pyEQUIB.svg?branch=master
    :target: https://travis-ci.org/equib/pyEQUIB
    :alt: Build Status
    
.. image:: https://ci.appveyor.com/api/projects/status/b3gw6vgf8s0vu8nv?svg=true
    :target: https://ci.appveyor.com/project/danehkar/pyequib
    :alt: Build Status
    
.. image:: https://coveralls.io/repos/github/equib/pyEQUIB/badge.svg?branch=master
    :target: https://coveralls.io/github/equib/pyEQUIB?branch=master
    :alt: Coverage Status
    
.. image:: https://img.shields.io/badge/license-GPL-blue.svg
    :target: https://github.com/equib/pyEQUIB/blob/master/LICENSE
    :alt: GitHub license
    
.. image:: https://img.shields.io/conda/vn/conda-forge/pyequib.svg
    :target: https://anaconda.org/conda-forge/pyequib
    :alt: Anaconda Cloud
    
.. image:: https://readthedocs.org/projects/pyequib/badge/?version=latest
    :target: https://pyequib.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
    
.. image:: https://img.shields.io/badge/python-2.7%2C%203.8-blue.svg
    :alt: Support Python versions 2.7 and 3.8
    
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4287576.svg
    :target: https://doi.org/10.5281/zenodo.4287576
    :alt: Zenodo
    
.. image:: http://joss.theoj.org/papers/10.21105/joss.02798/status.svg
    :target: https://doi.org/10.21105/joss.02798
    :alt: JOSS


Description
===========

The **pyEQUIB** library is a collection of `Python <https://www.python.org/>`_ programs developed to perform plasma diagnostics and abundance analysis using emission line fluxes measured in ionzed nebulae. It uses the `AtomNeb Python Package <https://github.com/atomneb/AtomNeb-py>`_ to read collision strengths and transition probabilities for collisionally excited lines (CEL), and recombination coefficients for recombination lines (RL). This Python package can be used to determine interstellar extinctions, electron temperatures, electron densities, and ionic abundances from the measured fluxes of emission lines. It mainly contains the follwing API functions written purely in Python: 

* **API functions for collisionally excited lines (CEL)** have been developed based on the algorithm of the FORTRAN program `EQUIB <http://adsabs.harvard.edu/abs/2016ascl.soft03005H>`_ written in FORTRAN by `Howarth & Adams (1981) <http://adsabs.harvard.edu/abs/1981ucl..rept.....H>`_. The program EQUIB calculates atomic level populations and line emissivities in statistical equilibrium in multi-level atoms for different physical conditions of the stratification layers where the chemical elements are ionized. Using the Python implementation of the program `EQUIB <http://adsabs.harvard.edu/abs/2016ascl.soft03005H>`_, electron temperatures, electron densities, and ionic abundances are determined from the measured fluxes of collisionally excited lines.

* **API functions for recombination lines (RL)** have been developed based on the algorithm of the recombination scripts by X. W. Liu and Y. Zhang included in the FORTRAN program `MOCASSIN <https://github.com/mocassin/MOCASSIN-2.0>`_. These API functiosn are used to determine ionic abundances from recombination lines for some heavy element ions.
 
* **API functions for reddening and extinctions** have been developed according to the methods of the reddening law functions from `STSDAS IRAF Package <http://www.stsci.edu/institute/software_hardware/stsdas>`_, which are used to obtain interstellar extinctions and deredden measured fluxes based on different reddening laws.


Installation
============

Dependent Python Packages
-------------------------

 This package requires the following packages:

    - `NumPy <https://numpy.org/>`_
    - `SciPy <https://scipy.org/scipylib/>`_
    - `AtomNeb <https://github.com/atomneb/AtomNeb-py/>`_
    
* To get this package with the AtomNeb FITS files, you can simply use ``git`` command as follows::

        git clone --recursive https://github.com/equib/pyEQUIB

To install the last version, all you should need to do is

.. code-block::

    $ python setup.py install

To install the stable version, you can use the preferred installer program (pip):

.. code-block::

    $ pip install pyequib

or you can install it from the cross-platform package manager *conda*:

.. code-block::

    $ conda install -c conda-forge pyequib

How to Use
==========

The Documentation of the Python functions provides in detail in the *API Documentation* (`equib.github.io/pyEQUIB/doc <https://equib.github.io/pyEQUIB/doc>`_). There are three main object units:

* **Collision Unit** which have the API functions for plasma diagnostics and abundance analysis of collisionally excited lines. Here are some examples of using *Collision* Unit:

    - *Temperature*::

        import pyequib
        import atomneb
        import os
        base_dir = 'externals/atomneb'
        data_dir = os.path.join('atomic-data', 'chianti70')
        atom_elj_file = os.path.join(base_dir,data_dir, 'AtomElj.fits')
        atom_omij_file = os.path.join(base_dir,data_dir, 'AtomOmij.fits')
        atom_aij_file = os.path.join(base_dir,data_dir, 'AtomAij.fits')
        
        atom = 's'
        ion = 'ii'
        s_ii_elj = atomneb.read_elj(atom_elj_file, atom, ion, level_num=5)
        s_ii_omij = atomneb.read_omij(atom_omij_file, atom, ion)
        s_ii_aij = atomneb.read_aij(atom_aij_file, atom, ion)
        
        upper_levels='1,2,1,3/'
        lower_levels='1,5/'
        density = np.float64(2550)
        line_flux_ratio=np.float64(10.753)
        temperature = pyequib.calc_temperature(line_flux_ratio=line_flux_ratio, density=density, 
                               upper_levels=upper_levels, lower_levels=lower_levels, 
                               elj_data=s_ii_elj, omij_data=s_ii_omij, aij_data=s_ii_aij)
        print("Electron Temperature:", temperature)

      which gives::
    
        Electron Temperature:       7920.2865

    - *Density*::

        import pyequib
        import atomneb
        import os
        base_dir = 'externals/atomneb'
        data_dir = os.path.join('atomic-data', 'chianti70')
        atom_elj_file = os.path.join(base_dir,data_dir, 'AtomElj.fits')
        atom_omij_file = os.path.join(base_dir,data_dir, 'AtomOmij.fits')
        atom_aij_file = os.path.join(base_dir,data_dir, 'AtomAij.fits')
        
        atom = 's'
        ion = 'ii'
        s_ii_elj = atomneb.read_elj(atom_elj_file, atom, ion, level_num=5)
        s_ii_omij = atomneb.read_omij(atom_omij_file, atom, ion)
        s_ii_aij = atomneb.read_aij(atom_aij_file, atom, ion)
        
        upper_levels='1,2/'
        lower_levels='1,3/'
        temperature=np.float64(7000.0)#
        line_flux_ratio=np.float64(1.506)#
        density = pyequib.calc_density(line_flux_ratio=line_flux_ratio, temperature=temperature, 
                                       upper_levels=upper_levels, lower_levels=lower_levels, 
                                       elj_data=s_ii_elj, omij_data=s_ii_omij, aij_data=s_ii_aij)
        print("Electron Density:", density)

      which gives::
      
        Electron Density:       2312.6395

    - *Ionic Abundance*::

        import pyequib
        import atomneb
        import os
        base_dir = 'externals/atomneb'
        data_dir = os.path.join('atomic-data', 'chianti70')
        data_rc_dir = os.path.join('atomic-data-rc')
        atom_elj_file = os.path.join(base_dir,data_dir, 'AtomElj.fits')
        atom_omij_file = os.path.join(base_dir,data_dir, 'AtomOmij.fits')
        atom_aij_file = os.path.join(base_dir,data_dir, 'AtomAij.fits')
        atom_rc_sh95_file = os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
        
        atom = 'h'
        ion = 'ii' # H I Rec
        hi_rc_data = atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)
        
        atom = 'o'
        ion = 'iii' # [O III]
        o_iii_elj = atomneb.read_elj(atom_elj_file, atom, ion, level_num=5) # read Energy Levels (Ej)
        o_iii_omij = atomneb.read_omij(atom_omij_file, atom, ion) # read Collision Strengths (Omegaij)
        o_iii_aij = atomneb.read_aij(atom_aij_file, atom, ion) # read Transition Probabilities (Aij)

        levels5007='3,4/'
        temperature=np.float64(10000.0)
        density=np.float64(5000.0)
        iobs5007=np.float64(1200.0)
        abb5007 = pyequib.calc_abundance(temperature=temperature, density=density, 
                                         line_flux=iobs5007, atomic_levels=levels5007,
                                         elj_data=o_iii_elj, omij_data=o_iii_omij, aij_data=o_iii_aij, 
                                         h_i_aeff_data=hi_rc_data['aeff'][0])
        print('N(O^2+)/N(H+):', abb5007)

      which gives::
      
        N(O^2+)/N(H+):   0.00041256231 
        
    - *Emissivity*::
    
        import pyequib
        import atomneb
        import os
        base_dir = 'externals/atomneb'
        data_dir = os.path.join('atomic-data', 'chianti70')
        data_rc_dir = os.path.join('atomic-data-rc')
        atom_elj_file = os.path.join(base_dir,data_dir, 'AtomElj.fits')
        atom_omij_file = os.path.join(base_dir,data_dir, 'AtomOmij.fits')
        atom_aij_file = os.path.join(base_dir,data_dir, 'AtomAij.fits')
        atom_rc_sh95_file = os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
        
        atom = 'h'
        ion = 'ii' # H I Rec
        hi_rc_data = atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)
        
        atom = 'o'
        ion = 'iii' # [O III]
        o_iii_elj = atomneb.read_elj(atom_elj_file, atom, ion, level_num=5) # read Energy Levels (Ej)
        o_iii_omij = atomneb.read_omij(atom_omij_file, atom, ion) # read Collision Strengths (Omegaij)
        o_iii_aij = atomneb.read_aij(atom_aij_file, atom, ion) # read Transition Probabilities (Aij)
        
        levels5007='3,4/'
        temperature=np.float64(10000.0)
        density=np.float64(5000.0)
        iobs5007=np.float64(1200.0)
        emis = pyequib.calc_emissivity(temperature=temperature, density=density, atomic_levels=levels5007, 
                                       elj_data=o_iii_elj, omij_data=o_iii_omij, aij_data=o_iii_aij)
        print('Emissivity(O III 5007):', emis)

      which gives::
      
        Emissivity(O III 5007):   3.6041012e-21
        

    - *Atomic Level Population*::

        import pyequib
        import atomneb
        import os
        base_dir = 'externals/atomneb'
        data_dir = os.path.join('atomic-data', 'chianti70')
        atom_elj_file = os.path.join(base_dir,data_dir, 'AtomElj.fits')
        atom_omij_file = os.path.join(base_dir,data_dir, 'AtomOmij.fits')
        atom_aij_file = os.path.join(base_dir,data_dir, 'AtomAij.fits')
        
        atom = 's'
        ion = 'ii'
        s_ii_elj = atomneb.read_elj(atom_elj_file, atom, ion, level_num=5)
        s_ii_omij = atomneb.read_omij(atom_omij_file, atom, ion)
        s_ii_aij = atomneb.read_aij(atom_aij_file, atom, ion)
        
        density = np.float64(1000)
        temperature=np.float64(10000.0)#
        nlj = pyequib.calc_populations(temperature=temperature, density=density, 
                                       elj_data=s_ii_elj, omij_data=s_ii_omij, aij_data=s_ii_aij)
        print('Populations:', nlj)

      which prints::
      
        Populations: 0.96992832 0.0070036315 0.023062261 2.6593671e-06 3.1277019e-06

    - *Critical Density*::
    
        import pyequib
        import atomneb
        import os
        base_dir = 'externals/atomneb'
        data_dir = os.path.join('atomic-data', 'chianti70')
        atom_elj_file = os.path.join(base_dir,data_dir, 'AtomElj.fits')
        atom_omij_file = os.path.join(base_dir,data_dir, 'AtomOmij.fits')
        atom_aij_file = os.path.join(base_dir,data_dir, 'AtomAij.fits')
        
        atom = 's'
        ion = 'ii'
        s_ii_elj = atomneb.read_elj(atom_elj_file, atom, ion, level_num=5)
        s_ii_omij = atomneb.read_omij(atom_omij_file, atom, ion)
        s_ii_aij = atomneb.read_aij(atom_aij_file, atom, ion)
        
        temperature=np.float64(10000.0)
        n_crit = pyequib.calc_crit_density(temperature=temperature, 
                                           elj_data=s_ii_elj, omij_data=s_ii_omij, aij_data=s_ii_aij)
        print('Critical Densities:', n_crit)

      which gives::
      
        Critical Densities: 0.0000000 5007.8396 1732.8414 1072685.0 2220758.1

    - *All Ionic Level Information*::
    
        import pyequib
        import atomneb
        import os
        base_dir = 'externals/atomneb'
        data_dir = os.path.join('atomic-data', 'chianti70')
        data_rc_dir = os.path.join('atomic-data-rc')
        atom_elj_file = os.path.join(base_dir,data_dir, 'AtomElj.fits')
        atom_omij_file = os.path.join(base_dir,data_dir, 'AtomOmij.fits')
        atom_aij_file = os.path.join(base_dir,data_dir, 'AtomAij.fits')
        atom_rc_sh95_file = os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
        
        atom = 'h'
        ion = 'ii' # H I Rec
        hi_rc_data = atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)
        
        atom = 'o'
        ion = 'iii' # [O III]
        o_iii_elj = atomneb.read_elj(atom_elj_file, atom, ion, level_num=5) # read Energy Levels (Ej)
        o_iii_omij = atomneb.read_omij(atom_omij_file, atom, ion) # read Collision Strengths (Omegaij)
        o_iii_aij = atomneb.read_aij(atom_aij_file, atom, ion) # read Transition Probabilities (Aij)
        
        temperature=np.float64(10000.0)
        density=np.float64(5000.0)
        pyequib.print_ionic(temperature=temperature, density=density,
                    elj_data=o_iii_elj, omij_data=o_iii_omij, aij_data=o_iii_aij,
                    h_i_aeff_data=hi_rc_data['aeff'][0])

      which gives::
      
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


* **Recombination Unit** which have the API functions for plasma diagnostics and abundance analysis of recombination lines. Here are some examples of using *Recombination* Unit:

    - *He+ Ionic Abundance*::

        import pyequib
        import atomneb
        import os
        base_dir = 'externals/atomneb'
        data_rc_dir = os.path.join('atomic-data-rc')
        atom_rc_he_i_file = os.path.join(base_dir,data_rc_dir, 'rc_he_ii_PFSD12.fits')
        atom_rc_sh95_file = os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
        
        atom = 'h'
        ion = 'ii' # H I
        h_i_rc_data = atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)
        
        atom = 'he'
        ion = 'ii' # He I
        he_i_rc_data = atomneb.read_aeff_he_i_pfsd12(atom_rc_he_i_file, atom, ion)
        
        temperature=np.float64(10000.0)
        density=np.float64(5000.0)
        he_i_4471_flux= 2.104
        linenum=10# 4471.50
        abund_he_i = pyequib.calc_abund_he_i_rl(temperature=temperature, density=density,
                                        linenum=linenum, line_flux=he_i_4471_flux,
                                        he_i_aeff_data=he_i_aeff_data, h_i_aeff_data=h_i_aeff_data)
        print('N(He^+)/N(H^+):', abund_he_i)

      which gives::
      
        N(He^+)/N(H^+):     0.040848393

    - *He++ Ionic Abundance*::

        import pyequib
        import atomneb
        import os
        base_dir = 'externals/atomneb'
        data_rc_dir = os.path.join('atomic-data-rc')
        atom_rc_sh95_file = os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
        
        atom = 'h'
        ion = 'ii' # H I
        h_i_rc_data = atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)
        
        atom = 'he'
        ion = 'iii' # He II
        he_ii_rc_data = atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)
        
        temperature=np.float64(10000.0)
        density=np.float64(5000.0)
        he_ii_4686_flux = 135.833
        abund_he_ii = pyequib.calc_abund_he_ii_rl(temperature=temperature, density=density,
                                          line_flux=he_ii_4686_flux,
                                          he_ii_aeff_data=he_ii_aeff_data, h_i_aeff_data=h_i_aeff_data)
        print('N(He^2+)/N(H^+):', abund_he_ii)

      which gives::
      
        N(He^2+)/N(H^+):      0.11228817

    - *C++ Ionic Abundance*::
    
        import pyequib
        import atomneb
        import os
        base_dir = 'externals/atomneb'
        data_rc_dir = os.path.join('atomic-data-rc')
        atom_rc_all_file = os.path.join(base_dir,data_rc_dir, 'rc_collection.fits')
        atom_rc_sh95_file = os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
        
        atom = 'c'
        ion = 'iii' # C II
        c_ii_rc_data = atomneb.read_aeff_collection(atom_rc_all_file, atom, ion)
        
        atom = 'h'
        ion = 'ii' # H I
        h_i_rc_data = atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)

        temperature=np.float64(10000.0)
        density=np.float64(5000.0)
        wavelength=6151.43
        c_ii_6151_flux = 0.028
        abund_c_ii = pyequib.calc_abund_c_ii_rl(temperature=temperature, density=density,
                                        wavelength=wavelength, line_flux=c_ii_6151_flux,
                                        c_ii_rc_data=c_ii_rc_data, h_i_aeff_data=h_i_aeff_data)
        print('N(C^2+)/N(H+):', abund_c_ii)

      which gives::
      
        N(C^2+)/N(H+):   0.00063404650 
      
    - *C3+ Ionic Abundance*::

        import pyequib
        import atomneb
        import os
        base_dir = 'externals/atomneb'
        data_rc_dir = os.path.join('atomic-data-rc')
        atom_rc_ppb91_file = os.path.join(base_dir,data_rc_dir, 'rc_PPB91.fits')
        atom_rc_sh95_file = os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
        
        atom = 'c'
        ion = 'iv' # C III
        c_iii_rc_data = atomneb.read_aeff_ppb91(atom_rc_ppb91_file, atom, ion)
        
        atom = 'h'
        ion = 'ii' # H I
        h_i_rc_data = atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)
        
        temperature=np.float64(10000.0)
        density=np.float64(5000.0)
        wavelength=4647.42
        c_iii_4647_flux = 0.107
        abund_c_iii = pyequib.calc_abund_c_iii_rl(temperature=temperature, density=density,
                                          wavelength=wavelength,
                                          line_flux=c_iii_4647_flux, c_iii_rc_data=c_iii_rc_data,
                                          h_i_aeff_data=h_i_aeff_data)
        print('N(C^3+)/N(H+):', abund_c_iii)

      which gives::
      
        N(C^3+)/N(H+):   0.00017502840

    - *N++ Ionic Abundance*::

        import pyequib
        import atomneb
        import os
        base_dir = 'externals/atomneb'
        data_rc_dir = os.path.join('atomic-data-rc')
        atom_rc_all_file = os.path.join(base_dir,data_rc_dir, 'rc_collection.fits')
        atom_rc_sh95_file = os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
        
        atom = 'n'
        ion = 'iii' # N II
        n_ii_rc_data = atomneb.read_aeff_collection(atom_rc_all_file, atom, ion)
        n_ii_rc_data_br = atomneb.read_aeff_collection(atom_rc_all_file, atom, ion, br=True)
        
        atom = 'h'
        ion = 'ii' # H I
        h_i_rc_data = atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)
        
        wavelength=4442.02
        n_ii_4442_flux = 0.017
        abund_n_ii = pyequib.calc_abund_n_ii_rl(temperature=temperature, density=density,
                                        wavelength=wavelength, line_flux=n_ii_4442_flux,
                                        n_ii_rc_br=n_ii_rc_data_br, n_ii_rc_data=n_ii_rc_data,
                                        h_i_aeff_data=h_i_aeff_data)
        print('N(N^2+)/N(H+):', abund_n_ii)

      which gives::
      
        N(N^2+)/N(H+):   0.00069297541

    - *N3+ Ionic Abundance*::
    
        import pyequib
        import atomneb
        import os
        base_dir = 'externals/atomneb'
        data_rc_dir = os.path.join('atomic-data-rc')
        atom_rc_ppb91_file = os.path.join(base_dir,data_rc_dir, 'rc_PPB91.fits')
        atom_rc_sh95_file = os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
        
        atom = 'n'
        ion = 'iv' # N III
        n_iii_rc_data = atomneb.read_aeff_ppb91(atom_rc_ppb91_file, atom, ion)
           
        atom = 'h'
        ion = 'ii' # H I
        h_i_rc_data = atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)
        
        wavelength=4640.64
        n_iii_4641_flux = 0.245
        abund_n_iii = pyequib.calc_abund_n_iii_rl(temperature=temperature, density=density,
                                          wavelength=wavelength, line_flux=n_iii_4641_flux,
                                          n_iii_rc_data=n_iii_rc_data, h_i_aeff_data=h_i_aeff_data)
        print('N(N^3+)/N(H+):', abund_n_iii)

      which gives::
      
        N(N^3+)/N(H+):   6.3366175e-05

    - *O++ Ionic Abundance*::

        import pyequib
        import atomneb
        import os
        base_dir = 'externals/atomneb'
        data_rc_dir = os.path.join('atomic-data-rc')
        atom_rc_all_file = os.path.join(base_dir,data_rc_dir, 'rc_collection.fits')
        atom_rc_sh95_file = os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
        
        atom = 'o'
        ion = 'iii' # O II
        o_ii_rc_data = atomneb.read_aeff_collection(atom_rc_all_file, atom, ion)
        o_ii_rc_data_br = atomneb.read_aeff_collection(atom_rc_all_file, atom, ion, br=True)
                   
        atom = 'h'
        ion = 'ii' # H I
        h_i_rc_data = atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)
        
        wavelength=4613.68
        o_ii_4614_flux = 0.009
        abund_o_ii = pyequib.calc_abund_o_ii_rl(temperature=temperature, density=density,
                                        wavelength=wavelength, line_flux=o_ii_4614_flux,
                                        o_ii_rc_br=o_ii_rc_data_br,
                                        o_ii_rc_data=o_ii_rc_data,
                                        h_i_aeff_data=h_i_aeff_data)              
        print('N(O^2+)/N(H+):', abund_o_ii)
        
      which gives::
      
        N(O^2+)/N(H+):    0.0018886330

    - *Ne++ Ionic Abundance*::

        import pyequib
        import atomneb
        import os
        base_dir = 'externals/atomneb'
        data_rc_dir = os.path.join('atomic-data-rc')
        atom_rc_all_file = os.path.join(base_dir,data_rc_dir, 'rc_collection.fits')
        atom_rc_sh95_file = os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
        
        atom = 'ne'
        ion = 'iii' # Ne II
        ne_ii_rc_data = atomneb.read_aeff_collection(atom_rc_all_file, atom, ion)
                   
        atom = 'h'
        ion = 'ii' # H I
        h_i_rc_data = atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)
        
        wavelength=3777.14
        ne_ii_3777_flux = 0.056
        abund_ne_ii = pyequib.calc_abund_ne_ii_rl(temperature=temperature, density=density,
                                          wavelength=wavelength, line_flux=ne_ii_3777_flux,
                                          ne_ii_rc_data=ne_ii_rc_data, h_i_aeff_data=h_i_aeff_data)
        print('N(Ne^2+)/N(H+):', Abund_ne_ii)

      which gives::
      
        N(Ne^2+)/N(H+):   0.00043376850


    - *He I Emissivity*::

        import pyequib
        import atomneb
        import os
        base_dir = 'externals/atomneb'
        data_rc_dir = os.path.join('atomic-data-rc')
        atom_rc_he_i_file = os.path.join(base_dir,data_rc_dir, 'rc_he_ii_PFSD12.fits')
        
        atom = 'he'
        ion = 'ii' # He I
        he_i_rc_data = atomneb.read_aeff_he_i_pfsd12(atom_rc_he_i_file, atom, ion)
        
        temperature=np.float64(10000.0)
        density=np.float64(5000.0)
        linenum=10# 4471.50
        emiss_he_i = pyequib.calc_emiss_he_i_rl(temperature=temperature, density=density,
                                        linenum=linenum, he_i_aeff_data=he_i_aeff_data)
        print('He I Emissivity:', emiss_he_i)

      which gives::
      
        He I Emissivity:   6.3822830e-26

    - *He II Emissivity*::
    
        import pyequib
        import atomneb
        import os
        base_dir = 'externals/atomneb'
        data_rc_dir = os.path.join('atomic-data-rc')
        atom_rc_sh95_file = os.path.join(base_dir,data_rc_dir, 'rc_SH95.fits')
         
        atom = 'he'
        ion = 'iii' # He II
        he_ii_rc_data = atomneb.read_aeff_sh95(atom_rc_sh95_file, atom, ion)

        temperature=np.float64(10000.0)
        density=np.float64(5000.0)
        emiss_he_ii = pyequib.calc_emiss_he_ii_rl(temperature=temperature, density=density,
                                          he_ii_aeff_data=he_ii_aeff_data)
        print('He II Emissivity:', emiss_he_ii)

      which gives::
      
        He II Emissivity:   1.4989134e-24

    - *C II Emissivity*::
    
        import pyequib
        import atomneb
        import os
        base_dir = 'externals/atomneb'
        data_rc_dir = os.path.join('atomic-data-rc')
        atom_rc_all_file = os.path.join(base_dir,data_rc_dir, 'rc_collection.fits')
        
        atom = 'c'
        ion = 'iii' # C II
        c_ii_rc_data = atomneb.read_aeff_collection(atom_rc_all_file, atom, ion)
        
        temperature=np.float64(10000.0)
        density=np.float64(5000.0)
        wavelength=6151.43
        emiss_c_ii = pyequib.calc_emiss_c_ii_rl(temperature=temperature, density=density,
                                        wavelength=wavelength, c_ii_rc_data=c_ii_rc_data)
        print('C II Emissivity:', emiss_c_ii)

      which gives::
      
        C II Emissivity:   5.4719511e-26
      
    - *C III Emissivity*::

        import pyequib
        import atomneb
        import os
        base_dir = 'externals/atomneb'
        data_rc_dir = os.path.join('atomic-data-rc')
        atom_rc_ppb91_file = os.path.join(base_dir,data_rc_dir, 'rc_PPB91.fits')
        
        atom = 'c'
        ion = 'iv' # C III
        c_iii_rc_data = atomneb.read_aeff_ppb91(atom_rc_ppb91_file, atom, ion)
        
        temperature=np.float64(10000.0)
        density=np.float64(5000.0)
        wavelength=4647.42
        emiss_c_iii = pyequib.calc_emiss_c_iii_rl(temperature=temperature, density=density,
                                          wavelength=wavelength,
                                          c_iii_rc_data=c_iii_rc_data)
        print('C III Emissivity:', emiss_c_iii)

      which gives::
      
        C III Emissivity:   7.5749632e-25

    - *N II Emissivity*::
    
        import pyequib
        import atomneb
        import os
        base_dir = 'externals/atomneb'
        data_rc_dir = os.path.join('atomic-data-rc')
        atom_rc_all_file = os.path.join(base_dir,data_rc_dir, 'rc_collection.fits')
        
        atom = 'n'
        ion = 'iii' # N II
        n_ii_rc_data = atomneb.read_aeff_collection(atom_rc_all_file, atom, ion)
        n_ii_rc_data_br = atomneb.read_aeff_collection(atom_rc_all_file, atom, ion, br=True)
        
        wavelength=4442.02
        emiss_n_ii = pyequib.calc_emiss_n_ii_rl(temperature=temperature, density=density,
                                        wavelength=wavelength,
                                        n_ii_rc_br=n_ii_rc_data_br, n_ii_rc_data=n_ii_rc_data)
        print('N II Emissivity:', emiss_n_ii)

      which gives::
      
        N II Emissivity:   3.0397397e-26

    - *N III Emissivity*::
    
        import pyequib
        import atomneb
        import os
        base_dir = 'externals/atomneb'
        data_rc_dir = os.path.join('atomic-data-rc')
        atom_rc_ppb91_file = os.path.join(base_dir,data_rc_dir, 'rc_PPB91.fits')
        
        atom = 'n'
        ion = 'iv' # N III
        n_iii_rc_data = atomneb.read_aeff_ppb91(atom_rc_ppb91_file, atom, ion)
        
        wavelength=4640.64
        emiss_n_iii = pyequib.calc_emiss_n_iii_rl(temperature=temperature, density=density,
                                          wavelength=wavelength, n_iii_rc_data=n_iii_rc_data)
        print('N III Emissivity:', emiss_n_iii)

      which gives::
      
        N III Emissivity:   4.7908644e-24

    - *O II Emissivity*::

        import pyequib
        import atomneb
        import os
        base_dir = 'externals/atomneb'
        data_rc_dir = os.path.join('atomic-data-rc')
        atom_rc_all_file = os.path.join(base_dir,data_rc_dir, 'rc_collection.fits')
        
        atom = 'o'
        ion = 'iii' # O II
        o_ii_rc_data = atomneb.read_aeff_collection(atom_rc_all_file, atom, ion)
        o_ii_rc_data_br = atomneb.read_aeff_collection(atom_rc_all_file, atom, ion, br=True)
        
        wavelength=4613.68
        emiss_o_ii = pyequib.calc_emiss_o_ii_rl(temperature=temperature, density=density,
                                        wavelength=wavelength,
                                        o_ii_rc_br=o_ii_rc_data_br, o_ii_rc_data=o_ii_rc_data)
        print('O II Emissivity:', emiss_o_ii)
        
      which gives::
      
        O II Emissivity:   5.9047319e-27

    - *Ne II Emissivity*::

        import pyequib
        import atomneb
        import os
        base_dir = 'externals/atomneb'
        data_rc_dir = os.path.join('atomic-data-rc')
        atom_rc_all_file = os.path.join(base_dir,data_rc_dir, 'rc_collection.fits')
        
        atom = 'ne'
        ion = 'iii' # Ne II
        ne_ii_rc_data = atomneb.read_aeff_collection(atom_rc_all_file, atom, ion)
        
        wavelength=3777.14
        emiss_ne_ii = pyequib.calc_emiss_ne_ii_rl(temperature=temperature, density=density,
                                          wavelength=wavelength, ne_ii_rc_data=ne_ii_rc_data)
        print('Ne II Emissivity:', emiss_ne_ii)

      which gives::
      
        Ne II Emissivity:   1.5996881e-25
        
* **Reddening Unit** which have the API functions for estimating logarithmic extinctions at H-beta and dereddening observed fluxes based on reddening laws and extinctions. Here are some examples of using *Reddening* Unit:

    - *Reddening Law Function*::

        import pyequib
        wavelength=6563.0
        r_v=3.1
        fl=pyequib.redlaw(wavelength, rv=r_v, ext_law='GAL')
        print('fl(6563):', fl)

      which gives::
      
        fl(6563):     -0.32013816

    - *Galactic Reddening Law Function based on Seaton (1979), Howarth (1983), & CCM (1983)*::

        import pyequib
        wavelength=6563.0
        r_v=3.1
        fl=pyequib.redlaw_gal(wavelength, rv=r_v)
        print('fl(6563):', fl)

      which gives::
      
        fl(6563):     -0.32013816

    - *Galactic Reddening Law Function based on Savage & Mathis (1979)*::

        import pyequib
        wavelength=6563.0
        fl=pyequib.redlaw_gal2(wavelength)
        print('fl(6563):', fl)

      which gives::
      
        fl(6563):     -0.30925984

    - *Reddening Law Function based on Cardelli, Clayton & Mathis (1989)*::
    
        import pyequib
        wavelength=6563.0
        r_v=3.1
        fl=pyequib.redlaw_ccm(wavelength, rv=r_v)
        prin('fl(6563):', fl)

      which gives::
      
        fl(6563):     -0.29756615

    - *Galactic Reddening Law Function based on Whitford (1958), Seaton (1977), & Kaler(1976)*::
    
        import pyequib
        wavelength=6563.0
        fl=pyequib.redlaw_jbk(wavelength)
        print('fl(6563):', fl)

      which gives::
      
        fl(6563):     -0.33113684

    - *Reddening Law Function based on Fitzpatrick & Massa (1990), Fitzpatrick (1999), Misselt (1999)*::
    
        import pyequib
        wavelength=6563.0
        r_v=3.1
        fmlaw='AVGLMC'
        fl=pyequib.redlaw_fm(wavelength, fmlaw=fmlaw, rv=r_v)
        print('fl(6563):', fl)

      which gives::
      
        fl(6563):     -0.35053032

    - *Reddening Law Function for the Small Magellanic Cloud*::
    
        import pyequib
        wavelength=6563.0
        fl=pyequib.redlaw_smc(wavelength)
        print('fl(6563):', fl)

      which gives::
      
        fl(6563):     -0.22659261

    - *Reddening Law Function for the Large Magellanic Cloud*::
    
        import pyequib
        wavelength=6563.0
        fl=pyequib.redlaw_lmc(wavelength)
        print('fl(6563):', fl)

      which gives::
      
        fl(6563):     -0.30871187

    - *Dereddening Absolute Flux*::

        import pyequib
        wavelength=6563.0
        m_ext=1.0
        flux=1.0
        ext_law='GAL'
        r_v=3.1
        flux_deredden=pyequib.deredden_relflux(wavelength, flux, m_ext, ext_law=ext_law, rv=r_v)
        print('dereddened flux(6563)', flux_deredden)

      which gives::
      
        dereddened flux(6563)       4.7847785

    - *Dereddening Relative Flux*::

        import pyequib
        wavelength=6563.0
        m_ext=1.0
        flux=1.0
        ext_law='GAL'
        r_v=3.1
        flux_deredden=pyequib.deredden_flux(wavelength, flux, m_ext, ext_law=ext_law, rv=r_v)
        print('dereddened flux(6563)', flux_deredden)

      which gives::
      
        dereddened flux(6563)      0.47847785

Documentation
=============

For more information on how to use the API functions from the pyEQUIB libray, please read the `API Documentation  <https://equib.github.io/pyEQUIB/doc>`_ published on `equib.github.io/pyEQUIB <https://equib.github.io/pyEQUIB>`_.


References
==========
* Danehkar, A. (2020). pyEQUIB Python Package, an addendum to proEQUIB: IDL Library for Plasma Diagnostics and Abundance Analysis. *J. Open Source Softw.*, **5**, 2798. doi:`10.21105/joss.02798 <https://doi.org/10.21105/joss.02798>`_ ads:`2020JOSS....5.2798D <https://ui.adsabs.harvard.edu/abs/2020JOSS....5.2798D>`_.

* Danehkar, A. (2018). proEQUIB: IDL Library for Plasma Diagnostics and Abundance Analysis. *J. Open Source Softw.*, **3**, 899. doi:`10.21105/joss.00899 <https://doi.org/10.21105/joss.00899>`_ ads:`2018JOSS....3..899D <https://ui.adsabs.harvard.edu/abs/2018JOSS....3..899D>`_.


Citation
========

Using **pyEQUIB** in a scholarly publication? Please cite thess papers:

.. code-block:: bibtex

   @article{Danehkar2020,
     author = {{Danehkar}, Ashkbiz},
     title = {pyEQUIB Python Package, an addendum to proEQUIB: IDL Library for Plasma Diagnostics and Abundance Analysis},
     journal = {Journal of Open Source Software},
     volume = {5},
     number = {55},
     pages = {2798},
     year = {2020},
     doi = {10.21105/joss.02798}
   }

   @article{Danehkar2018,
     author = {{Danehkar}, Ashkbiz},
     title = {proEQUIB: IDL Library for Plasma Diagnostics and Abundance Analysis},
     journal = {Journal of Open Source Software},
     volume = {3},
     number = {32},
     pages = {899},
     year = {2018},
     doi = {10.21105/joss.00899}
   }

Learn More
==========

==================  =============================================
**Documentation**   https://pyequib.readthedocs.io/
**Repository**      https://github.com/equib/pyEQUIB
**Issues & Ideas**  https://github.com/equib/pyEQUIB/issues
**Conda-Forge**     https://anaconda.org/conda-forge/pyequib
**PyPI**            https://pypi.org/project/pyequib/
**DOI**             `10.21105/joss.02798 <https://doi.org/10.21105/joss.02798>`_
**Archive**         `10.5281/zenodo.4287576 <https://doi.org/10.5281/zenodo.4287576>`_
==================  =============================================
