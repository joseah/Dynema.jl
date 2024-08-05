function sample_index(rng, n_cells)

    index = collect(1:n_cells)
    index = sample(rng, index, length(index); replace=true)
    return index

end
