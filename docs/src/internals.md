```@meta
CurrentModule = Dynema
```

# Internals

These types are not exported, but are documented here for contributors and
for users who want to work with `expand_geno`'s return value directly.

```@docs
Dynema.ExpandedGeno
Dynema.ExpandedGenoView
```

!!! note
    `Dynema.DynemaModel`, the struct returned by [`map_locus`](@ref), is
    intentionally accessed only through its `get_*`/`set_*` accessors (see
    the [API Reference](functions.md)) rather than by touching its fields
    directly, so that the internal layout can change without breaking user
    code.
