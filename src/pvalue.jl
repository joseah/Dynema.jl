
function pass(stat, boot, stattype)

    if stattype == "z"
        minimum([sum(stat .> boot), sum(stat .< boot)])
    elseif stattype == "χ²"
        sum(boot .>= stat)
    end

end

function compute_pvalue(stat, boot, stattype)

    pass_obs = pass(stat, boot, stattype)

    if stattype == "z"

        pass_obs == 0 ? 2 / (length(boot) + 1) : 2 * pass_obs / length(boot)

    elseif stattype == "χ²"

        pass_obs == 0 ? 1 / (length(boot) + 1) : pass_obs / length(boot)

    end

end


"""
    crvetest(R, r; resp, scores, betas, A, clustid, small = false)

Analytical CRVE score (Lagrange multiplier) test via
`WildBootTests.scoretest`. `getci`/`getplot`/`getdist` are disabled
explicitly: they default to `true` in WildBootTests and trigger a full
confidence-curve grid search per call whose output Dynema never reads.
"""
function crvetest(R, r; resp::AbstractVector, scores::AbstractMatrix, betas::AbstractVector,
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


"""
    crvetest_direct(R, A, Sg, clustshare)

Direct implementation of the analytical CRVE score test for the shared-null
fast path, replicating `WildBootTests.scoretest`'s ml-mode statistic
(`reps = 0`, `imposenull = true`, `scorebs = true`, `small = false`,
one-way clustering) without the library's per-call setup:

  - `numer = (A R')' S`, with `S` the total score sum;
  - cluster-level `K = Sg * (A R')`, recentered by each cluster's share of
    observations, so `denom = K'K` is the CRVE meat transformed by `A R'`;
  - q = 1: `z = numer / sqrt(denom)`, `p = ccdf(Chisq(1), z²)` (two-tailed);
    q > 1: `χ² = numer' denom⁻¹ numer`, `p = ccdf(Chisq(q), χ²)`.

`Sg` is the G × k matrix of per-cluster score sums; `clustshare[g]` is
cluster g's share of observations. Returns `(stat, p, stattype)` in the same
shape `crvetest` results are consumed.
"""
function crvetest_direct(R::AbstractMatrix, A::AbstractMatrix,
                         Sg::AbstractMatrix, clustshare::AbstractVector)

    ARt = A * R'                          # k × q
    S = vec(sum(Sg, dims = 1))            # total score sum (k)
    numer = ARt' * S                      # q

    K = Sg * ARt                          # G × q
    K .-= clustshare * (S' * ARt)         # recenter, as WildBootTests does for scorebs
    denom = K' * K                        # q × q CRVE denominator

    q = size(R, 1)
    if q == 1
        z = numer[1] / sqrt(denom[1, 1])
        return (stat = z, p = ccdf(Chisq(1), z^2), stattype = "z")
    else
        χ² = dot(numer, Symmetric(denom) \ numer)
        return (stat = χ², p = ccdf(Chisq(q), χ²), stattype = "χ²")
    end

end
