module NetworkBridge

using Graphs, DataFrames, StatsBase, Statistics, Printf

# --- 1️⃣ Detect Communities ---
function detect_communities(g::SimpleGraph)
    labs = nothing
    try
        @eval using CommunityDetection
        labs = CommunityDetection.labels(CommunityDetection.louvain(g))
    catch
        @warn "CommunityDetection.jl not found; using Graphs.label_propagation instead."
        labs = label_propagation(g)
    end
    labs isa Tuple ? first(labs) : labs
end

# --- 2️⃣ Participation Coefficients ---
function participation_coefficients(g::SimpleGraph, comms::AbstractVector{<:Integer})
    n = nv(g)
    uniq = unique(comms)
    remap = Dict(uniq[i] => i for i in eachindex(uniq))
    C = length(uniq)
    comm = [remap[c] for c in comms]
    deg = degree(g)
    k_is = zeros(Float64, n, C)

    for e in edges(g)
        u, v = src(e), dst(e)
        k_is[u, comm[v]] += 1
        k_is[v, comm[u]] += 1
    end

    pc = similar(deg, Float64)
    @inbounds for i in 1:n
        k = max(deg[i], 1)
        pc[i] = 1.0 - sum((k_is[i, :] ./ k).^2)
    end
    pc
end

# --- 3️⃣ Normalize helper ---
normalize(x::AbstractVector) = (x .- minimum(x)) ./ (maximum(x) - minimum(x) + eps())

# --- 4️⃣ Summarize Bridges ---
function summarize_bridges(
    g::SimpleGraph;
    top_n::Int=10,
    targets_df::Union{Nothing,DataFrame}=nothing,
    id2idx::Union{Nothing,Dict}=nothing,
    labels::Union{Nothing,AbstractVector}=nothing
)
    println("🔹 Computing community structure and bridge scores...")

    # --- Compute metrics ---
    comms = detect_communities(g)
    pc    = participation_coefficients(g, comms)
    btw   = betweenness_centrality(g)
    score = normalize(pc) .+ normalize(btw)

    # --- Build Base DataFrame ---
    df = DataFrame(
        node = 1:nv(g),
        bridge_score = score,
        participation = pc,
        betweenness = btw
    )

    # --- Add original Facebook IDs ---
    if id2idx !== nothing
        idx2id = Dict(v => k for (k, v) in id2idx)
        df.original_id = [idx2id[i] for i in 1:nv(g)]
    else
        df.original_id = missing
    end

    # --- Attach Page Type ---
    if labels !== nothing
        # Directly aligned from preprocessing
        df.page_type = labels
    elseif targets_df !== nothing && id2idx !== nothing
        idx2id = Dict(v => k for (k, v) in id2idx)
        page_lookup = Dict(targets_df.id .=> targets_df.page_type)
        df.page_type = [get(page_lookup, idx2id[i], "Unknown") for i in 1:nv(g)]
    else
        df.page_type = ["Unknown" for _ in 1:nv(g)]
    end

    # --- Sort and Display ---
    sort!(df, :bridge_score, rev=true)
    top = first(df, min(top_n, nrow(df)))

    println("\n====== TOP STRUCTURAL BRIDGES ======")
    show(top, allrows=true, allcols=true, truncate=80)

    return (top=top, comms=comms, btw=btw, pc=pc, df=df)
end

export summarize_bridges, detect_communities, participation_coefficients

end # module
