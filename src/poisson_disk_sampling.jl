function fast_neighbors(data, n_neighbors; metric=Euclidean())
    """
    Calculates the nearest neighboring points for each point in `data`
    
    Parameters
    ----------
    data : a coordinate matrix of shape <num coords> x <num dimensions>
    n_neighbors : the number of nearest neighbors
    metric : distance metric for calculating neighbor distances.
             See `Distances`` package for metric options.
    
    Returns
    ------
    matrix of shape <num coords> x <n_neighbors>,
      with row i containing the neighbors of data[i]
    """
    

    
    graph = nndescent(transpose(data), n_neighbors - 1, metric)
    indices, _ = knn_matrices(graph)

    indices = vcat(transpose(1:size(indices, 2)), indices)
    indices = [row for row in eachcol(indices)]
    
    return indices

end


function graph_poisson_disk(rng, neighbors, n_pseudocells, n_candidates=100)
    """
    Calculates Poisson disk samples on a graph specified by `neighbors`
    
    Parameters
    ----------
    rng : Random number generator
    neighbors : a 2D list specifying the neighbors of each node in the graph,
                neighbors[i] is a list of the neighbors of node i
    n_pseudocells : the number of Poisson disk samples to generate
    n_candidates : the number of candidates to choose between when drawing a new sample,
                   by default is 100    
    Returns
    ------
    a vector with each element being a sublist containing the nodes in a Poisson disk sample
    """
    
    # Select cell at random
    sample = rand(rng, 1:length(neighbors))  # randomly choose a node
    samples = [sample]
    
    # Remove it from the available samples
    available_samples = setdiff(1:length(neighbors), sample)  # indices except `sample`
    
    # Get neighbors of cell `sample`
    included_samples = neighbors[sample]  # neighbors of the sample node
    
    # Create sets for all neighbors
    neighbor_sets = [Set(i) for i in neighbors]
    
    while length(samples) < n_pseudocells && length(available_samples) > 0
        # Choose `n_candidates` out of the available cells (indices)
        sample_candidates = rand(rng, available_samples, n_candidates)  # select candidate nodes
        sample_candidates = unique(sample_candidates)  # get unique candidates
        
        # Add to whole list of candidates
        sample_candidates_all = sample_candidates
        
        # For each cell in `sample_candidates`, check if any of their neighbors are in the `included_samples` set
        sample_candidates = filter(i -> isempty(Set(neighbor_sets[i]) ∩ included_samples), sample_candidates)
        
        # If no candidates, reset to all sample candidates
        if length(sample_candidates) == 0
            sample_candidates = sample_candidates_all
        end
        
        # Retrieve neighbors for all `sample_candidates`
        super_neighbors = []
        for i in sample_candidates
            push!(super_neighbors, neighbors[i])
        end
        
        # Get intersection of neighbors between `included_samples` and all `sample_candidates` neighbors
        sample_density = [length(intersect(Set(i), included_samples)) / length(i) for i in super_neighbors]
        
        # Select the candidate with the minimum density
        new_sample = sample_candidates[argmin(sample_density)]
        
        # Add the cell with the lowest overlap density to the `samples` list
        push!(samples, new_sample)
        
        # Update the available samples by removing the neighbors of `new_sample`
        available_samples = setdiff(available_samples, neighbors[new_sample])
        
        # Update `included_samples` with the neighbors of the new sample
        included_samples = Set(included_samples) ∪ Set(neighbors[new_sample])
        

    end
    

    return samples

end
