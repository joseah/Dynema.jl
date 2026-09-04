# ---------------------------------------------------------------------------- #
#                            DynemaWildBootTestsExt                            #
# ---------------------------------------------------------------------------- #
#
# Package extension providing the WildBootTests-backed functionality --
# adaptive score bootstrapping (`scorebootstrap`) and the generic analytical
# CRVE score test (`crvetest`, used as a fallback and for the fast path's
# one-time self-check). Loaded automatically when the optional WildBootTests
# package is installed and `using WildBootTests` has run alongside
# `using Dynema`; without it, Dynema's own direct score test covers all
# standard analytical mapping and the stubs in src/ explain how to enable
# these paths.

module DynemaWildBootTestsExt

using Dynema
using WildBootTests
using Dynema: pass_counts, combine_pass, compute_pvalue, beta_approx_pvalue
using Random: AbstractRNG

function __init__()
    Dynema.HAS_WILDBOOTTESTS[] = true
end

# ------------------------- Analytical CRVE score test ------------------------ #

# getci/getplot/getdist are disabled explicitly: they default to true in
# WildBootTests and trigger a full confidence-curve grid search per call
# whose output Dynema never reads.
function Dynema.crvetest(R::AbstractMatrix, r::AbstractVector;
                    resp::AbstractVector, scores::AbstractMatrix, betas::AbstractVector,
                    A::AbstractMatrix, clustid::AbstractMatrix, small::Bool = false)

    scoretest(R, r;
                resp = resp,
                scores = scores,
                beta = betas,
                A = A,
                clustid = clustid,
                ml = true,
                scorebs = true,
                small = small,
                getci = false,
                getplot = false,
                getdist = false)

end

# --------------------------- Adaptive score bootstrap ------------------------ #


function Dynema.scorebootstrap(R::AbstractMatrix, r::AbstractVector;
                    resp::AbstractVector, scores::AbstractMatrix, betas::AbstractVector,
                    A::AbstractMatrix, clustid::AbstractMatrix, small::Bool = false,
                    rng::AbstractRNG, ptype::Symbol, B::Vector{Int64})


            # getci/getplot are disabled explicitly: they default to true in
            # WildBootTests and trigger a confidence-curve grid search per
            # call (many extra bootstrap evaluations) whose output Dynema
            # never reads. getdist stays on: the distribution drives the
            # adaptive escalation below.
            test = wildboottest(R, r;
                                resp = resp,
                                scores = scores,
                                beta = betas,
                                A = A,
                                clustid = clustid,
                                ml = true,
                                scorebs = true,
                                imposenull = true,
                                small = small,
                                getci = false, getplot = false, getdist = true,
                                rng = rng, ptype = ptype, reps = B[1] - 1)

            # Extract results. Tail counts are accumulated incrementally per
            # round instead of re-concatenating and re-scanning the whole
            # distribution each round (which was quadratic in rounds).
            statistic = teststat(test)
            stattype = test.stattype
            dists = [dist(test)[1, :]]
            c1, c2 = pass_counts(statistic, dists[1], stattype)
            total = length(dists[1])
            counts = combine_pass(c1, c2, stattype)

            # ---------------- Remaining rounds of bootstrapping if needed --------------- #

            if counts <= 20

                for j in 2:length(B)

                    test = wildboottest(R, r;
                                resp = resp,
                                scores = scores,
                                beta = betas,
                                A = A,
                                clustid = clustid,
                                ml = true,
                                scorebs = true,
                                imposenull = true,
                                small = false,
                                getci = false, getplot = false, getdist = true,
                                rng = rng, ptype = ptype,
                                reps = B[j])


                    bootdist_i = dist(test)[1, :]
                    push!(dists, bootdist_i)
                    d1, d2 = pass_counts(statistic, bootdist_i, stattype)
                    c1 += d1; c2 += d2
                    total += length(bootdist_i)
                    counts = combine_pass(c1, c2, stattype)

                    if counts > 20
                        break
                    end

                end

            end

            # Compute final p-values; concatenate the distribution once.
            pval = compute_pvalue(counts, total, stattype)
            bootdist = length(dists) == 1 ? dists[1] : reduce(vcat, dists)

            # FastQTL-style beta approximation over the bootstrap
            # distribution: smooth, and not capped at the empirical
            # resolution floor of ~1/B (see beta_approx_pvalue).
            pbeta = beta_approx_pvalue(statistic, bootdist, stattype, size(R, 1))

            return (pval, pbeta, bootdist)


end

end # module
