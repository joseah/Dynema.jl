
"""
    scorebootstrap(R, r; resp, scores, betas, A, clustid, small = false, rng, ptype, B)

Adaptive score bootstrapping via `WildBootTests.wildboottest`. WildBootTests
is an *optional* dependency: the real method lives in the
DynemaWildBootTestsExt extension, loaded automatically when the WildBootTests
package is installed and loaded alongside Dynema. This stub errors with
install instructions otherwise.
"""
function scorebootstrap(R, r; kwargs...)
    error("Score bootstrapping (boot = true) needs the optional WildBootTests package. " *
          "Install and load it with:\n" *
          "    using Pkg; Pkg.add(\"WildBootTests\"); using WildBootTests\n" *
          "alongside `using Dynema` (the dynema-map CLI does this automatically when --boot is used).")
end
