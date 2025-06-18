# BLR.jl

***A fast and flexible toolkit for modeling the broad-line region (BLR) in Julia***

## Quickstart 
Define a cloud model with syntax like: 
```julia
mClouds = BLR.cloudModel(100_000,μ=500.,β=1.,F=0.5,θₒ=30/180*π,γ=1.,ξ=1.,i=0.,I=BLR.IsotropicIntensity,v=BLR.vCircularCloud,rescale=1e-5,τ=0.0)
```

## Referencing
If you find this code useful in your work, please cite it as:
```bibtex
lorem ipsum
```

## Contributing
If you would like to contribute to the package, please open a pull request on the GitHub. For bug reports and feature requests, please open an issue on the GitHub. 