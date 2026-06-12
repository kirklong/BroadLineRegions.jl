# API
Reference for `BroadLineRegions.jl`'s public interface. 

!!! note
    No methods are exported by default into the global namespace to prevent overlap with other modules, and you must prepend the module name to all methods to access them. `BroadLineRegions.jl` exports itself as both `BroadLineRegions` and `BLR`, so both of these prefixes are equivalent, i.e. `BroadLineRegions.model == BLR.model`. If prepending this to function calls annoys you you can always manually import whatever you desire into the global space with syntax like: ```using BLR: DiskWindModel, cloudModel```

## Performance notes

Several internal optimizations affect how you should interact with models:

- **Result caching.** Each `model` memoizes the arrays that `getVariable` gathers from its rings
  (intensities, velocities, delays, …) in `model.cache`, so repeated profile/transfer-function
  calculations on the same model are much faster than the first. The package's own mutating
  functions (`removeNaN!`, `raytrace!`, `zeroDiskObscuredClouds!`, `removeDiskObscuredClouds!`,
  model combination with `+`) manage this automatically — but **if you mutate ring fields
  directly** (e.g. `m.rings[1].I .= 0.0`) **you must call `reset!(m)` afterwards** so subsequent
  calculations see the new values. Set `m.cache = nothing` to disable caching for a model.

- **Cached 3D coordinates.** Every point's final system coordinates `(x, y, z)` are computed once at
  construction and stored on the ring (fields `x`, `y`, `z`; accessor `getXYZ`). 
  If you mutate a ring's geometry (`r`, `ϕ₀`, `i`, `rot`, `θₒ`, `reflect`),
  set `ring.x = nothing; ring.y = nothing; ring.z = nothing` to force recomputation.
  In performance-critical custom code prefer `rotate3D_scalar` (allocation-free tuple return) over
  `rotate3D` (allocates a `Vector`).

- **Multithreading.** Disk-wind model construction parallelizes across rings/pixels when Julia is
  started with multiple threads (`julia -t N`); results are bit-identical at any thread count.
  Custom `I`/`v` functions passed to the constructors
  must be thread-safe to benefit (the built-in ones are). Cloud model generation are currently deliberately
  *not* threaded so that seeded models (`rng=...`) reproduce identically.

- **Delays in combined models.** For combined (multi-submodel) models all time delays are computed
  with the general geometric formula `t = η(r − x)` (see `tCloud`) for every point to ensure the 
  most general formula is applied to be correct with all custom implementations. Thus it is critical in any combined 
  models for performance that `ring.x` is set properly for all submodels. 

## Full documentation

```@autodocs
Modules = [BroadLineRegions]
Order = [:function,:type,:macro]
Pages = ["structs.jl","clouds.jl","operators.jl","intensity.jl","velocity.jl","profiles.jl","transfer.jl","raytrace.jl","util.jl"]
```