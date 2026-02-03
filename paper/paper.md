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

Here $\Delta V$ is the characteristic velocity of the gas, which is often estimated from the width of the broad-line and can be measured from even just a single spectrum of the source. $R_{\rm{BLR}}$ is the characteristic distance from the supermassive black hole to the BLR, which has historically been constrained via reverberation mapping (RM) [@Peterson_RMReview] and more recently has been constrained with interferometric observations with GRAVITY [@GRAVITY19;@GRAVITY24_z2.4]. $f$ is the so-called "virial factor", which encodes extra information about the geometry and kinematics of the BLR, and $G$ is Newton's gravitational constant. While the black hole mass can be very roughly obtained from simple estimates of each of these quantities, to accurately measure the masses of supermassive black holes across cosmic time we must interpret data through a model of the BLR. There are thus both significant measurement and model-dependent uncertainties encoded in $\Delta V$, $R_{\rm{BLR}}$ and $f$. Sources with supermassive black hole masses measured from their broad-line regions underlie the measurements of many further reported black hole measurements, as they are used as calibrators in scaling relations [@GRAVITY_R-L_2024]. 

Beyond measuring black hole masses, if we want to understand the environment surrounding supermassive black holes we must characterize the kinematics and geometry of the BLR. There is a wealth of exciting new data from velocity-resolved reverberation mapping campaigns as well as interferometric observations with GRAVITY that present new opportunities to unravel the fundamental nature of the BLR, but to do so requires a degree of modeling to explain the data.

# Statement of need

[`BroadLineRegions.jl`](https://github.com/kirklong/BroadLineRegions.jl) enables fast and (more importantly) flexible modeling of the BLR in Julia, allowing for the quick creation of theoretical model BLRs for comparison to data or for making theoretical predictions. There are many possible models of the BLR, and in reality the BLR may be more complicated than any single component model can describe [@Long_2025]. `BroadLineRegions.jl` enables researchers to compare and contrast models as well as combine multiple models with easy syntax, allowing researchers to easily test their own bespoke models against and in concert with others. To most accurately measure the masses of supermassive black holes as well as to better understand the fundamental nature of the BLR we must understand what classes of models best fit observations. Additionally, many existing BLR modeling tools are closed source and thus `BroadLineRegions.jl`'s open-source philosophy will enable anyone to use it as a tool to uncover the true nature of the gas surrounding supermassive black holes.

# State of the field
Currently researchers working to understand the BLR through modeling in general take one of two approaches: either they use simple, largely kinematic models to quickly fit the observations, or they attempt more in-depth and physically realistic simulations of the fluids and radiative transfer to try to produce similar results to observations from the ground up. While the second approach is admirable, it is computationally infeasible to apply such detailed modeling to actually fit the data observed in sources, thus this code follows the first approach &mdash; being largely a kinematic, simple modeling tool designed to flexibly fit and match data with a more limited set of physics. 

Currently the most popular kinematic modeling code used to fit the BLR is the `CARAMEL` code as described in @Pancoast2014 and @Pancoast2011. While a powerful tool to flexibly fit BLR line profiles, as pointed out in @Long_2025 this model of the BLR cannot fully explain more recent velocity-resolved reverberation mapping measurements. This fact, in addition to the fact that all existing BLR modeling codes are relatively narrow and fixed in scope/physics, motivated the development of `BroadLineRegions.jl` to more flexibly model the BLR and enable the combination and testing of different classes of models against each other. Additionally `CARAMEL` is closed-source, thus contributing directly is not possible, and many different research groups in the field have their own custom-tailored versions they use in their fitting. This further motivates `BroadLineRegion.jl`'s open-source design philosophy, which allows for researchers to independently contribute their models and ideas to this important tool that we must use to better characterize the geometry and kinematics of the BLR as well as most accurately measure supermassive black hole masses. Finally, `BroadLineRegions.jl` allows for easier interfacing with GRAVITY and RM data products, whereas `CARAMEL` and most other BLR modeling codes are focused solely on the reverberation problem.

# Software design 
`BroadLineRegions.jl` was built from the ground-up to be highly modular and flexible, such that end-users can replace certain parts/steps with custom/novel alternatives if they wish to. We use several mutable structs to accomplish this, the most important being the `ring` struct, which represents a ring on a "camera" the that represents the observer's view of the BLR.
A single `ring` is created by defining the relevant physical quantities (intensity, velocity, distance from the source, etc.) either by hand or through helper functions, which can be user supplied and thus enabling great flexibility for the end-user. A `model` is then a combination of camera `rings`, and the code has been designed such that one can easily combine models together to build up the BLR piece by piece. By default we supply constructors for "cloud" and "disk-wind" BLR models &mdash; two popular classes of model in the literature, with the "cloud" model similar in implementation to the `CARAMEL` BLR implementation &mdash; but users are not restricted to using these constructors, and can easily create their own or modify pieces of these defaults. For example, if a user wanted to model the BLR as a hybrid model of a disk-wind and cloud setup, this can be accomplished in just three lines of code: one to initialize the cloud model, one to initialize the disk-wind model, and a third to create the combined model, an approach demonstrated in @Long_2025. Or perhaps an end-user may like the cloud model of the BLR as described in @Pancoast2014 but wishes to alter the intensity from the default setting &mdash; this is again easy to accomplish as a result of the modular nature of the code, as they can just pass a custom intensity function to the cloud model constructor and `BroadLineRegions.jl` will call this custom intensity function instead of the default when generating the model. 

While this approach is very flexible, it does come with the design tradeoff of some lost performance due to requiring that the structs be mutable instead of immutable. This mutability enables the easy combination of models and allows for the user to overwrite quantities after model creation. The documentation attempts to make this distinction very clearly to warn users that models are not permanent, but of course this may also introduce some user confusion if they accidentally overwrite a variable and expected it to stay constant. Despite this modest tradeoff, we believe that this approach is required scientifically because this software was motivated by the idea that current BLR modeling codes are too rigid to explain all the data, and thus we designed `BroadLineRegions.jl` to intentionally use this flexibility as a distinct feature advantage. 

# Research impact statement
[`BroadLineRegions.jl`](https://github.com/kirklong/BroadLineRegions.jl) has already been used to better model the BLR and expose the systematic uncertainties in our choice of models [@Long_2023; @Long_2025; @Long_2026]. As demonstrated in @Long_2025, the BLR is likely more complicated than current existing single-component kinematic models can describe, thus motivating the newly possible approaches that can be taken with this code. To truly understand the environments around supermassive black holes as well as most accurately measure their masses we must have models that are capable of matching all of the data products, thus `BroadLineRegions.jl` has the potential to be a next-generation open-source tool for accomplishing this goal. 

# AI usage disclosure
Most of the software as described in this work was completed before the ubiquity of current LLM coding tools. Development was completed largely using the VSCode IDE with co-pilot functionality after mid-2024, but even after enabling co-pilot we never used generative AI for more than simple extended (~few lines at most) autocompletions (i.e. filling in a simple for loop, etc.) and never for design choices or extended code samples. No generative AI was used in the preparation of this manuscript or the `BroadLineRegions.jl` documentation.

# Acknowledgements
We acknowledge support from Prof. Jason Dexter, who has supervised and advised the author throughout his Ph.D. and whose NSF grants AST-1909711 and AST-2307983 have supported the author.

# References