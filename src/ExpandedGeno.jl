
struct ExpandedGeno{T<:AbstractMatrix, R<:AbstractVector{<:Integer},
                    RN<:Union{Nothing,AbstractVector}, CN<:Union{Nothing,AbstractVector}} <: AbstractMatrix{eltype(T)}
    mat::T
    rows::R
    rownames::RN
    colnames::CN
end


"""
    expand_geno(mat::AbstractMatrix, rows::AbstractVector;
                 rownames=nothing, colnames=nothing)

Memory-efficient wrapper for genotype matrices that allows lazy row and column subsetting.
Instead of building a full matrix at the single-cell level, it provides a lazy genotype 
matrix representation compatible with downstream operations. Genotype data at the single-cell
level is only materialized when subsetting a single column or row.

# Arguments
- `mat::AbstractMatrix`: donor × SNP genotype matrix (e.g., `Matrix{Float64}` or `SparseMatrixCSC`).
- `rows::AbstractVector`: mapping from single-cell to donor rows.
- `rownames::Union{Nothing,String,AbstractVector}=nothing`: optional cell ids
- `colnames::Union{Nothing,String,AbstractVector}=nothing`: optional SNP names.

# Indexing
- **2D indexing only**: `E[i, j]`, `E[i, :]`, `E[:, j]`, `E[I, J]`.
- **Single-index access is forbidden**: `E[i]` throws an error.
- Materializes **single rows** or **single columns** only.
- Submatrices return `ExpandedGenoView` (lazy, no copying).

# Examples
```julia
# Construct an ExpandedGeno
E = expand_geno(geno, indices;
                 rownames = cell_id,
                 colnames = snp_names)

# Access a single element
val = E[10, 5]        # Float64

# Access a single row
row = E[10, :]         # Vector{Float64}

# Access a single column
col = E[:, 5]          # Vector{Float64}

# Lazy submatrix (returns ExpandedGenoView)
view = E[1:100, 5:10]

# Access row and column names
rnames = names(E, 1)
cnames = names(E, 2)
```
"""
function expand_geno(mat, rows; rownames=nothing, colnames=nothing)
    rn = rownames isa String ? [rownames] : rownames === nothing ? nothing : rownames
    cn = colnames isa String ? [colnames] : colnames === nothing ? nothing : colnames
    return ExpandedGeno(mat, rows, rn, cn)
end

Base.size(E::ExpandedGeno) = (length(E.rows), size(E.mat, 2))
Base.IndexStyle(::Type{<:ExpandedGeno}) = IndexCartesian()
Base.eltype(E::ExpandedGeno) = eltype(E.mat)

# Single element
@inline Base.getindex(E::ExpandedGeno, i::Int, j::Int) = @inbounds E.mat[E.rows[i], j]

# Single row, all columns
Base.getindex(E::ExpandedGeno, i::Int, ::Colon) = @inbounds E.mat[E.rows[i], :]

# All rows, single column
Base.getindex(E::ExpandedGeno, ::Colon, j::Int) = @inbounds E.mat[E.rows, j]

# Vector of rows, single column
Base.getindex(E::ExpandedGeno, I::AbstractVector{<:Integer}, j::Int) =
    @inbounds E.mat[E.rows[I], j]

# Single row, vector of columns
Base.getindex(E::ExpandedGeno, i::Int, J::AbstractVector{<:Integer}) =
    @inbounds E.mat[E.rows[i], J]

# All rows, vector of columns → lazy view
Base.getindex(E::ExpandedGeno, ::Colon, J::AbstractVector{<:Integer}) =
    ExpandedGenoView(E.mat, E.rows, J,
                     E.rownames,
                     E.colnames === nothing ? nothing : E.colnames[J])

# Vector of rows, vector of columns → lazy view
Base.getindex(E::ExpandedGeno, I::AbstractVector{<:Integer}, J::AbstractVector{<:Integer}) =
    ExpandedGenoView(E.mat, E.rows[I], J,
                     E.rownames === nothing ? nothing : E.rownames[I],
                     E.colnames === nothing ? nothing : E.colnames[J])

# forbid single-index linear access
function Base.getindex(E::ExpandedGeno, I::Union{Int, AbstractVector})
    throw(ArgumentError("Single-index access not allowed; use row,column indexing, e.g. E[1:10, :]"))
end

Base.names(E::ExpandedGeno, dim::Int) =
    dim == 1 ? E.rownames : dim == 2 ? E.colnames :
    throw(ArgumentError("dimension must be 1 or 2"))


"""
    ExpandedGenoView(mat::AbstractMatrix, rows::AbstractVector,
                     cols::AbstractVector;
                     rownames=nothing, colnames=nothing)

A lazy view of an `ExpandedGeno`, representing a submatrix without copying data.

# Arguments
- `mat::AbstractMatrix`: underlying donor × SNP matrix.
- `rows::AbstractVector`: row indices of the view.
- `cols::AbstractVector`: column indices of the view.
- `rownames::Union{Nothing,AbstractVector}=nothing`: optional subset of row names.
- `colnames::Union{Nothing,AbstractVector}=nothing`: optional subset of column names.

# Indexing
- **2D indexing only**: `V[i, j]`, `V[i, :]`, `V[:, j]`, `V[I, J]`.
- Materializes **single rows** or **single columns** only.
- Supports `names(view, 1)` and `names(view, 2)`.

# Examples
```julia
# Create a lazy submatrix
view = E[1:100, 5:10]    # returns ExpandedGenoView

# Materialize a single row from the view
row = view[1, :]          # Vector{Float64}

# Materialize a single column from the view
col = view[:, 1]          # Vector{Float64}

# Access row and column names from the view
rnames = names(view, 1)
cnames = names(view, 2)
```
"""
struct ExpandedGenoView{T<:AbstractMatrix, R<:AbstractVector{<:Integer},
                        C<:AbstractVector{<:Integer},
                        RN<:Union{Nothing,AbstractVector}, CN<:Union{Nothing,AbstractVector}} <: AbstractMatrix{eltype(T)}
    mat::T
    rows::R
    cols::C
    rownames::RN
    colnames::CN
end

Base.size(E::ExpandedGenoView) = (length(E.rows), length(E.cols))
Base.IndexStyle(::Type{<:ExpandedGenoView}) = IndexCartesian()
Base.eltype(E::ExpandedGenoView) = eltype(E.mat)

@inline Base.getindex(E::ExpandedGenoView, i::Int, j::Int) = @inbounds E.mat[E.rows[i], E.cols[j]]

Base.names(E::ExpandedGenoView, dim::Int) =
    dim == 1 ? E.rownames : dim == 2 ? E.colnames :
    throw(ArgumentError("dimension must be 1 or 2"))
