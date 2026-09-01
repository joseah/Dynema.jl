# ---------------------------------------------------------------------------- #
#                                  vcf_geno.jl                                 #
# ---------------------------------------------------------------------------- #
#
# Shared VCF -> donor x variant dosage matrix extraction logic, `include()`d
# by both dynema_extract_geno.jl (extract-to-file) and dynema_map.jl
# (extract-on-the-fly before mapping). Not a registered package -- both
# scripts activate the same bin/Project.toml environment, so this file just
# needs to be included after that environment is active.
#
# Assumes the VCF already carries per-sample genotype dosages (FORMAT field
# `DS`) and/or genotype probabilities (FORMAT field `GP`) -- as produced by
# standard imputation pipelines (Minimac, IMPUTE2/5, Beagle, etc.). Hard-call
# genotypes (`GT`) are not read; a variant with neither `DS` nor `GP` is
# skipped. Region seeking is delegated to the `tabix` executable bundled by
# `htslib_jll` (a Julia package artifact -- no system htslib install needed);
# all GP/DS-to-dosage parsing happens here in Julia.

using DataFrames
using CSV
using Statistics
using htslib_jll

const VCF_SAMPLE_START = 10 # VCF fixed columns: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT, then samples

"""
    resolve_region(; chr, tss, tss_file, gene, window) -> (chr, tss, start_pos, end_pos)

Resolves the chromosome/TSS either directly from `chr`/`tss`, or by looking
`gene` up in `tss_file` (a table with columns gene_id, chr, tss), then
applies `window` to get the [start_pos, end_pos] region to query.
"""
function resolve_region(; chr, tss, tss_file, gene, window::Int)

    if gene !== nothing
        tss_file === nothing && error("--gene requires --tss-file to look its TSS up")
        ann = CSV.read(tss_file, DataFrame)
        for c in ("gene_id", "chr", "tss")
            c in names(ann) || error("--tss-file is missing required column '$c'")
        end
        rows = findall(==(gene), ann.gene_id)
        isempty(rows) && error("Gene '$gene' not found in --tss-file")
        length(rows) > 1 && error("Gene '$gene' matches $(length(rows)) rows in --tss-file; expected exactly one")
        chr = string(ann.chr[rows[1]])
        tss = Int(ann.tss[rows[1]])
    else
        (chr === nothing || tss === nothing) &&
            error("Either --gene (with --tss-file), or both --chr and --tss, must be provided to locate the cis-window")
    end

    start_pos = max(1, tss - window)
    end_pos = tss + window

    return chr, tss, start_pos, end_pos

end

"""
    run_tabix(vcf, chr, start_pos, end_pos) -> (header_fields, data_lines)

Uses the `tabix` executable bundled by `htslib_jll` to fetch the VCF header
(for sample names) and the data lines overlapping the region in a single
indexed query.
"""
function run_tabix(vcf::AbstractString, chr::AbstractString, start_pos::Int, end_pos::Int)

    isfile(vcf) || error("VCF file not found: $vcf")
    (isfile(vcf * ".tbi") || isfile(vcf * ".csi")) ||
        error("No tabix index found for $vcf. Run `tabix -p vcf $vcf` first (see README.md for how to do this with the htslib_jll-bundled tabix, with no system install needed).")

    region = "$(chr):$(start_pos)-$(end_pos)"
    bullet("querying region $region with tabix (via htslib_jll)...")

    out = try
        htslib_jll.tabix() do tabix_exe
            read(`$tabix_exe -h $vcf $region`, String)
        end
    catch e
        error("tabix failed to query $vcf for region $region: $e")
    end

    lines = split(out, '\n'; keepempty = false)
    header_idx = findfirst(l -> startswith(l, "#CHROM"), lines)
    header_idx === nothing && error("Could not find a #CHROM header line in tabix output for $vcf")
    header_fields = split(lines[header_idx], '\t')
    data_lines = filter(l -> !startswith(l, "#"), lines)

    return header_fields, data_lines

end

"""
    variant_dosage(fields, sample_start, field_pref) -> (variant_id, dosage_or_nothing, reason_or_nothing)

Parses one VCF data line into a variant id and a per-sample dosage vector,
preferring GP or DS per `field_pref` ("auto", "GP", "DS"). Returns
`dosage = nothing` with a skip `reason` if the variant should be dropped
(multiallelic, or missing the requested field).
"""
function variant_dosage(fields::Vector{<:AbstractString}, sample_start::Int, field_pref::AbstractString)

    chrom, pos, id, ref, alt = fields[1], fields[2], fields[3], fields[4], fields[5]
    variant_id = id == "." ? "$(chrom)_$(pos)_$(ref)_$(alt)" : id

    occursin(",", alt) && return variant_id, nothing, "multiallelic"

    format_keys = split(fields[9], ':')
    idx_gp = findfirst(==("GP"), format_keys)
    idx_ds = findfirst(==("DS"), format_keys)

    use_field = if field_pref == "GP"
        idx_gp === nothing ? nothing : :GP
    elseif field_pref == "DS"
        idx_ds === nothing ? nothing : :DS
    else # auto: prefer GP, fall back to DS
        idx_gp !== nothing ? :GP : (idx_ds !== nothing ? :DS : nothing)
    end

    use_field === nothing && return variant_id, nothing, "neither GP nor DS present (FORMAT=$(fields[9]))"

    idx = use_field == :GP ? idx_gp : idx_ds
    sample_fields = @view fields[sample_start:end]
    n = length(sample_fields)
    dosage = Vector{Union{Missing,Float64}}(undef, n)

    for (i, sf) in enumerate(sample_fields)
        subfields = split(sf, ':')
        raw = idx <= length(subfields) ? subfields[idx] : "."
        dosage[i] = if raw == "." || isempty(raw)
            missing
        elseif use_field == :DS
            parse(Float64, raw)
        else # GP: three genotype probabilities (0/0, 0/1, 1/1) -> expected dosage
            p = parse.(Float64, split(raw, ','))
            length(p) == 3 ? (p[2] + 2 * p[3]) : missing
        end
    end

    return variant_id, dosage, nothing

end

"""
    extract_geno_dataframe(; vcf, chr=nothing, tss=nothing, tss_file=nothing, gene=nothing,
                            window=1_000_000, field="auto", samples_file=nothing,
                            donor_col="donor_id", maf=0.0, max_missing=0.1) -> DataFrame

High-level, reusable VCF cis-window extraction: resolves the region, queries
the VCF with tabix, converts GP/DS to dosage, applies sample remapping and
MAF/missingness filtering, and returns a DataFrame with one donor-id column
(named `donor_col`) and one column per retained variant -- the same shape
both `dynema_map.jl --geno` and `dynema_extract_geno.jl --out` use, so the
result can either be written straight to a file or handed directly to
`map_locus` in-process.
"""
function extract_geno_dataframe(; vcf::AbstractString,
                                 chr::Union{Nothing,AbstractString} = nothing,
                                 tss::Union{Nothing,Int} = nothing,
                                 tss_file::Union{Nothing,AbstractString} = nothing,
                                 gene::Union{Nothing,AbstractString} = nothing,
                                 window::Int = 1_000_000,
                                 field::AbstractString = "auto",
                                 samples_file::Union{Nothing,AbstractString} = nothing,
                                 donor_col::AbstractString = "donor_id",
                                 maf::Float64 = 0.0,
                                 max_missing::Float64 = 0.1)

    section("Extracting genotypes from VCF")
    bullet("file: $vcf")

    chr, tss, start_pos, end_pos = resolve_region(chr = chr, tss = tss, tss_file = tss_file, gene = gene, window = window)
    bullet("cis-window: $(chr):$(start_pos)-$(end_pos) (TSS $(tss) +/- $(window) bp)")

    header_fields, data_lines = run_tabix(vcf, chr, start_pos, end_pos)
    vcf_samples = String.(header_fields[VCF_SAMPLE_START:end])

    # ------------------------------ Sample mapping ----------------------------- #

    donor_ids = vcf_samples
    keep_idx = collect(1:length(vcf_samples))
    if samples_file !== nothing
        smap = CSV.read(samples_file, DataFrame)
        for c in ("vcf_id", "donor_id")
            c in names(smap) || error("--samples is missing required column '$c'")
        end
        id_to_donor = Dict(zip(smap.vcf_id, smap.donor_id))
        keep_idx = findall(s -> haskey(id_to_donor, s), vcf_samples)
        isempty(keep_idx) && error("None of the VCF sample names matched 'vcf_id' in --samples")
        n_unmatched = length(vcf_samples) - length(keep_idx)
        n_unmatched > 0 && @warn "$n_unmatched VCF sample(s) not found in --samples were dropped"
        donor_ids = [id_to_donor[vcf_samples[i]] for i in keep_idx]
    end
    n_samples = length(keep_idx)
    bullet("samples: $(length(vcf_samples)) in VCF, $n_samples retained after sample matching")

    out_df = DataFrame()
    out_df[!, donor_col] = donor_ids

    if isempty(data_lines)
        @warn "No variants found in $(chr):$(start_pos)-$(end_pos); returning donor-id-only table"
        return out_df
    end

    # ------------------------------ Parse variants ------------------------------ #

    variant_ids = String[]
    columns = Vector{Float64}[]
    n_multiallelic = 0
    n_no_field = 0
    n_high_missing = 0

    for line in data_lines
        fields = split(line, '\t')
        vid, dosage, reason = variant_dosage(fields, VCF_SAMPLE_START, field)

        if dosage === nothing
            reason == "multiallelic" ? (n_multiallelic += 1) : (n_no_field += 1)
            continue
        end

        dosage = dosage[keep_idx]
        n_missing = count(ismissing, dosage)
        if n_missing / n_samples > max_missing
            n_high_missing += 1
            continue
        end
        if n_missing > 0
            m = mean(skipmissing(dosage))
            dosage = coalesce.(dosage, m)
        end
        dosage = Float64.(dosage)

        maf_i = mean(dosage) / 2
        maf_i = min(maf_i, 1 - maf_i)
        maf_i < maf && continue

        push!(variant_ids, vid)
        push!(columns, dosage)
    end

    bullet("variants in region: $(length(data_lines))")
    bullet("skipped (multiallelic): $n_multiallelic", indent = 2)
    bullet("skipped (missing GP/DS field): $n_no_field", indent = 2)
    bullet("skipped (missingness > $(max_missing)): $n_high_missing", indent = 2)
    bullet("retained after MAF >= $(maf) filter: $(length(variant_ids))", indent = 2)

    for (vid, col) in zip(variant_ids, columns)
        out_df[!, vid] = col
    end

    return out_df

end
