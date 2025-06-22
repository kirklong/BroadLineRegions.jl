# API
Reference for `BroadLineRegions.jl`'s public interface. 

!!! note
    No methods are exported by default into the global namespace to prevent overlap with other modules, and you must prepend the module name to all methods to access them. `BroadLineRegions.jl` exports itself as both `BroadLineRegions` and `BLR`, so both of these prefixes are equivalent, i.e. `BroadLineRegions.model == BLR.model`. If prepending this to function calls annoys you you can always manually import whatever you desire into the global space with syntax like: ```using BLR: DiskWindModel, cloudModel```

## Full documentation

```@autodocs
Modules = [BroadLineRegions]
Order = [:function,:type,:macro]
Pages = ["structs.jl","clouds.jl","operators.jl","intensity.jl","velocity.jl","profiles.jl","transfer.jl","raytrace.jl","util.jl"]
```