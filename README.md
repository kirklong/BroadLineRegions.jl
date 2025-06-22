# BLR.jl

<!-- Tidyverse lifecycle badges, see https://www.tidyverse.org/lifecycle/ Uncomment or delete as needed. -->
![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
[![build](https://github.com/kirklong/BLR.jl/workflows/CI/badge.svg)](https://github.com/kirklong/BLR.jl/actions?query=workflow%3ACI)
![Static Badge](https://img.shields.io/badge/docs-dodgerblue?link=https%3A%2F%2Fwww.kirklong.space%2FBLR.jl)
<!-- travis-ci.com badge, uncomment or delete as needed, depending on whether you are using that service. -->
<!-- [![Build Status](https://travis-ci.com/kirklong/BLR.jl.svg?branch=master)](https://travis-ci.com/kirklong/BLR.jl) -->
<!-- NOTE: Codecov.io badge now depends on the token, copy from their site after setting up -->
<!-- Documentation -- uncomment or delete as needed -->
<!--
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://kirklong.github.io/BLR.jl/stable)
[![Documentation](https://img.shields.io/badge/docs-master-blue.svg)](https://kirklong.github.io/BLR.jl/dev)
-->
<!-- Aqua badge, see test/runtests.jl -->
<!-- [![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl) -->
# BLR.jl

A fast, flexible, and open-source toolkit for modeling the broad-line region (BLR) in Julia.

**Under heavy construction**: things are liable to change rapidly, there may be bugs, and things may not be fully documented.

## Installation 
### Julia
```julia
using Pkg
Pkg.add("BLR.jl")
```
Or install from the GitHub repo directly:
```julia
using Pkg
Pkg.add("https://github.com/kirklong/BLR.jl")
```
### Python
To access `BLR.jl` from within your Python installation, first you must install [`JuliaCall`](https://juliapy.github.io/PythonCall.jl/stable/juliacall/)
After successfully installing `JuliaCall` you can add `BLR.jl` to your new Julia installation in Python with: 
```python
from juliacall import Main as jl
from juliacall import Pkg as jlPkg
jlPkg.add("BLR.jl") #or use the github link
```

## Quickstart 
While this code is designed to be very flexible and modular such that you can implement your own bespoke models of the BLR easier, two popular models of the BLR are included as default models which one can play with. 

To generate a "cloud" model similar to that of Pancoast+ 2011 and 2014, use syntax like:
```julia
mClouds = BLR.cloudModel(100_000,μ=500.,β=1.,F=0.5,θₒ=30/180*π,γ=1.,ξ=1.,i=0.,
        I=BLR.IsotropicIntensity,v=BLR.vCircularCloud,rescale=1e-5,τ=0.0)
```

To generate a "disk-wind" model similar to that of Chiang and Murray 1996 and 1997 following the prescription laid out in Long+ 2023 and 2025 use syntax like:
```julia
mDisk = BLR.DiskWindModel(500.,5.,1.,30/180*π,nr=24,nϕ=48,scale=:log,
        f1=1.0,f2=1.0,f3=1.0,f4=1.0,reflect=false,τ=5.)
```

You can of course define fully custom models by passing your own "rings" (on a camera) and optionally any default profiles, a camera struct, and submodel start indices:

```julia
mCustom = BLR.model(myCustomRings,nothing,myCustomCamera,[1]) #a custom model with no profiles, a user-defined camera for raytracing/visualization, and with no submodels
```

Models can be combined simply by writing `mCombined = m1 + m2`. 

Generate profiles (i.e. line, phase, delay, or whatever else your heart desires) for models with syntax like:
```julia
p = BLR.getProfile(m,:line) #generate line profile with default parameters
BLR.setProfile!(m,p) #optinally store the profile in model data structure 
```

A few default visualization recipes exist as well:
```julia
profiles = BLR.profile(m) #plot all of the model's stored profiles 
img = BLR.image(m,:I) #make an image of the model intensity, can pass any other parameter as well to "image" them
geometry = BLR.plot3d(m) #visualize the geometry of the system in a 3d plot, can also color points according to any parameter
```

For more detailed examples, see the [Usage and Examples](https://www.kirklong.space/BLR.jl/dev/usage_examples/) page.

Full documentation is available in the [API](https://www.kirklong.space/BLR.jl/dev/api/) section.

## Referencing
If you find this code useful in your work, please cite it as:
```bibtex
lorem ipsum
```

## Contributing
If you would like to contribute to the package, please open a pull request on the [GitHub](https://github.com/kirklong/BLR.jl). For bug reports and feature requests, please open an issue on the GitHub. 
