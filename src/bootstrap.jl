
function scorebootstrap(R, r; resp::AbstractVector, scores::AbstractMatrix, betas::AbstractVector,
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

            # Extract results
            statistic = teststat(test)
            bootdist = dist(test)[1, :]
            stattype = test.stattype
            counts = pass(statistic, bootdist, stattype)

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
                    bootdist = vcat(bootdist, bootdist_i)
                    counts = pass(statistic, bootdist, stattype)

                    if counts > 20
                        break
                    end

                end

            end

            # Compute final p-value
            pval = compute_pvalue(statistic, bootdist, stattype)


            return (pval, bootdist)


end
