module NetworkBridge

using Graphs, DataFrames, StatsBase

# --- 1) Community detection (Louvain/Leiden if available, fallback to label propagation)
"""
    detect_communities(g::SimpleGraph) -> Vector{Int}

Returns a vector of community labels (1..K). Tries CommunityDetection.jl (Louvain/Leiden),
falls back to Graphs.jl's `label_propagation`.
"""
function detect_communities(g::SimpleGraph)
    labs = nothing
    try
        @eval using CommunityDetection
        # Prefer Louvain; if you want Leiden instead, swap the next line.
        labs = CommunityDetection.labels(CommunityDetection.louvain(g))
    catch
        @warn "CommunityDetection.jl not found; using Graphs.label_propagation"
        labs = label_propagation(g)
    end
    labs isa Tuple ? first(labs) : labs
end

# --- 2) Participation coefficient
# PC_i = 1 - Σ_s (k_is / k_i)^2, where k_is is degree of i to comm s
function participation_coefficients(g::SimpleGraph, comms::AbstractVector{<:Integer})
    n = nv(g)
    # normalize labels to 1..K for safe indexing
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
        k = max(deg[i], 1)                        # avoid divide-by-zero
        pc[i] = 1.0 - sum((k_is[i, :] ./ k).^2)   # participation coefficient
    end
    pc
end

# --- 3) Bridge ranking = z(betweenness) + z(participation)
"""
    bridge_table(g::SimpleGraph; top_n=10)

Returns a DataFrame with node, community, betweenness, participation, and bridge_score,
sorted by bridge_score (descending).
"""
function bridge_table(g::SimpleGraph; top_n::Int=10)
    comms = detect_communities(g)
    btw   = betweenness_centrality(g)
    pc    = participation_coefficients(g, comms)

    # z-score both, then sum
    score = zscore(btw) .+ zscore(pc)

    df = DataFrame(
        node          = 1:nv(g),
        community     = comms,
        betweenness   = btw,
        participation = pc,
        bridge_score  = score,
    )
    sort!(df, :bridge_score, rev=true)
    first(df, min(top_n, nrow(df)))
end

"""
    summarize_bridges(g::SimpleGraph; top_n=10, meta::Union{Nothing,DataFrame}=nothing)

Prints a tiny summary and returns the ranked table. If you pass `meta` with an `:id`
column (matching node ids) and optional `:page_name`/`:page_type`, they’ll be joined.
"""
function summarize_bridges(g::SimpleGraph; top_n::Int=10, meta::Union{Nothing,DataFrame}=nothing)
    comms = detect_communities(g)
    btw   = betweenness_centrality(g)
    pc    = participation_coefficients(g, comms)
    score = zscore(btw) .+ zscore(pc)

    println("Nodes: $(nv(g))   Edges: $(ne(g))   Communities: $(length(unique(comms)))")
    println("Betweenness: mean=$(round(mean(btw), digits=4))  max=$(round(maximum(btw), digits=4))")

    df = DataFrame(node=1:nv(g), community=comms, betweenness=btw,
                   participation=pc, bridge_score=score)
    sort!(df, :bridge_score, rev=true)
    top = first(df, min(top_n, nrow(df)))

    if meta !== nothing && :id ∈ names(meta)
        keep = intersect([:id, :page_name, :page_type], names(meta))
        ren  = Dict(:id=>:node)
        top  = leftjoin(top, rename(select(meta, keep), ren; ignore=true), on=:node)
    end

    println("\nTop $(nrow(top)) structural bridges (high betweenness + high participation):")
    show(select(top, names(top)); allrows=true, truncate=100)
    println("\nInterpretation:")
    println("• Participation captures how evenly a node’s ties spread across communities.")
    println("• High betweenness + high participation ⇒ strong cross-community connector.")
    return top
end

export detect_communities, participation_coefficients, bridge_table, summarize_bridges

end # module
