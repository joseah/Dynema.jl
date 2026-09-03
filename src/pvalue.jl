
"""
    pass_counts(stat, boot, stattype) -> (c1, c2)

Allocation-free tail counts of `stat` against a bootstrap distribution chunk:
for "z", the (greater, less) counts (whose running sums stay additive across
chunks, unlike their min); for "χ²", `(count(boot .>= stat), 0)`.
"""
function pass_counts(stat, boot, stattype)
    if stattype == "z"
        (count(b -> stat > b, boot), count(b -> stat < b, boot))
    else
        (count(b -> b >= stat, boot), 0)
    end
end

# Final tail count from accumulated (c1, c2) totals
combine_pass(c1, c2, stattype) = stattype == "z" ? min(c1, c2) : c1

function pass(stat, boot, stattype)
    c1, c2 = pass_counts(stat, boot, stattype)
    combine_pass(c1, c2, stattype)
end

# From a precomputed tail count and total replication count
function compute_pvalue(pass_obs::Integer, nboot::Integer, stattype)

    if stattype == "z"

        pass_obs == 0 ? 2 / (nboot + 1) : 2 * pass_obs / nboot

    elseif stattype == "χ²"

        pass_obs == 0 ? 1 / (nboot + 1) : pass_obs / nboot

    end

end

compute_pvalue(stat, boot, stattype) =
    compute_pvalue(pass(stat, boot, stattype), length(boot), stattype)


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
