---
title: 'BroadLineRegions.jl: A fast and flexible toolkit for modeling the broad-line region (BLR) in Julia.'
tags:
  - Julia
  - astronomy
  - BLR
  - active galactic nuclei
authors:
  - name: Kirk Long
    orcid: 0009-0003-8293-2185
    equal-contrib: true
    affiliation: 1
affiliations:
 - name: JILA, University of Colorado Boulder, USA
   index: 1
date: 01 July 2025
bibliography: paper.bib

---

# Summary

Quasars and active galactic nuclei (AGN) show a remarkable degree of line-broadening, an observational hallmark whose origin is assumed to be due to the fast motions of line-emitting gas near the central supermassive black hole. The region from which this line emission originates is aptly named the broad-line region (BLR). Measuring and constraining properties of the BLR thus provides one of the only opportunities to directly constrain the masses of supermassive black holes outside our local universe as well as better understand the environments directly surrounding them. For example, assuming the gas is virialized we in principle need to measure just a few things about the BLR to get the black hole mass [@PetersonReview]: 

$$M_{\rm{BH}} = f\frac{R_{\rm{BLR}}\left(\Delta V\right)^2}{G}$$

Here $\Delta V$ is the characteristic velocity of the gas, which is often estimated from the width of the broad-line and can be measured from even just a single spectrum of the source. $R_{\rm{BLR}}$ is the characteristic distance from the supermassive black hole to the BLR, which has historically been constrained via reverberation mapping (RM) [@Peterson_RMReview] and more recently has been constrained with interferometric observations with GRAVITY [@GRAVITY19,@GRAVITY24_z2.4]. $f$ is the so-called "virial factor", which encodes extra information about the geometry and kinematics of the BLR, and $G$ is Newton's gravitational constant. While the black hole mass can be very roughly obtained from simple estimates of each of these quantities, to accurately measure the masses of supermassive black holes across cosmic time we must interpret data through a model of the BLR. There are thus both significant measurement and model-dependent uncertainties encoded in $\Delta V$, $R_{\rm{BLR}}$ and $f$. Sources with supermassive black hole masses measured from their broad-line regions underlie the measurements of many further reported black hole measurements, as they are used as calibrators in scaling relations [@GRAVITY_R-L_2024]. 

Beyond measuring black hole masses, if we want to understand the environment surrounding supermassive black holes we must characterize the kinematics and geometry of the BLR. There is a wealth of exciting new data from velocity-resolved reverberation mapping campaigns as well as interferometric observations with GRAVITY that present new opportunities to unravel the fundamental nature of the BLR, but to do so requires a degree of modeling to explain the data.

# Statement of need

[`BroadLineRegions.jl`](https://github.com/kirklong/BroadLineRegions.jl) enables fast and flexible modeling of the BLR in Julia, allowing for the quick creation of theoretical model BLRs for comparison to data or for making theoretical predictions. There are many possible models of the BLR, and in reality the BLR may be more complicated than any single component model can describe [@Long_2025]. `BroadLineRegions.jl` enables researchers to compare and contrast models as well as combine multiple models with easy syntax, allowing researchers to easily test their own bespoke models against and in concert with others. To most accurately measure the masses of supermassive black holes as well as to better understand the fundamental nature of the BLR we must understand what classes of models best fit observations. [`BroadLineRegions.jl`](https://github.com/kirklong/BroadLineRegions.jl) has already been used to this end in the literature, both in @Long_2023:2023 and @Long_2025:2025.

# Acknowledgements

We acknowledge support from Prof. Jason Dexter, who has supervised and advised the author throughout his Ph.D. and whose NSF grants AST-1909711 and AST-2307983 have supported the author.

# References