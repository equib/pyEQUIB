---
title: "pyEQUIB Python Package, an addendum to proEQUIB: IDL Library for Plasma Diagnostics and Abundance Analysis"
tags:
  - python
  - astrophysics
  - gaseous nebulae
  - plasma diagnostics
  - abundance analysis
authors:
  - name: Ashkbiz Danehkar
    orcid: 0000-0003-4552-5997
    affiliation: "1, 2, 3"
affiliations:
 - name: Department of Physics and Astronomy, Macquarie University, Sydney, NSW 2109, Australia
   index: 1
 - name: Harvard-Smithsonian Center for Astrophysics, 60 Garden Street, Cambridge, MA 02138, USA 
   index: 2
 - name: Department of Astronomy, University of Michigan, 1085 S. University Avenue, Ann Arbor, MI 48109, USA 
   index: 3
date: 20 October 2020
bibliography: paper.bib
---

# Addendum

`pyEQUIB` is a pure Python open-source package containing several application programming interface (API) functions that can be employed for plasma diagnostics and abundance analysis of nebular emission lines. This package is a Python implementation of the IDL library `proEQUIB` [@Danehkar:2018b] that is coupled to the IDL library `AtomNeb` [@Danehkar:2019]. The collisional excitation and recombination units of this package need to have the energy levels, collision strengths, transition probabilities, and recombination coefficients, which can be retrieved from the Python package `AtomNeb` for _Atomic Data of Ionized Nebulae_ [@Danehkar:2020]. The API functions of this package can be used to deduce the electron temperature, electron concentration, chemical elements from CELs and Rls, and the interstellar extinction from the Balmer decrements emitted from ionized gaseous nebulae. This package can simply be used by astronomers, who are familiar with the high-level, general-purpose programming language Python.

This package requires the Python packages `NumPy` [@Walt:2011; @Harris:2020], `SciPy` [@Virtanen:2020], and `AtomNeb` [@Danehkar:2020]. This package is released under the GNU General Public License. The source code is publicly available on the GitHub platform. The latest version of this package can be installed directly from its repository on the GitHub, and its stable version from the Python Package Index (PyPi) via ``pip install pyequib`` or alternatively from the Anaconda Python package distributor via ``conda install -c conda-forge pyequib``. The online documentation, tutorials and examples are available on its GitHub page (https://equib.github.io/pyEQUIB/) and its Read the Docs documentation page (https://pyequib.readthedocs.io/).

# Acknowledgements

The author acknowledges the support of Research Excellence Scholarship from Macquarie University.

# References
