module NetworkBridge

using Graphs, DataFrames, Statistics, StatsBase, Printf

# ---------------------------
# STRUCTURAL BRIDGES
# ---------------------------
# functions: edge_betweenness, bridge_edges, bridging_nodes, etc.

# --- community labels
function community_labels(g::SimpleGraph)
    try
        @eval using CommunityDetection
        return CommunityDetection.labels(CommunityDetection.louvain(g))
    catch
        @warn "CommunityDetection.jl not found; using label_propagation from Graphs.jl"
        return label_propagation(g)
    end
end


# --- compute edge betweenness centrality for all edges
function edge_betweenness(g::SimpleGraph)
    bc = edge_betweenness_centrality(g)
    edges = collect(edges(g))
    return DataFrame(src = [src(e) for e in edges],
                     dst = [dst(e) for e in edges],
                     betweenness = bc)
end

# --- return edges whose removal increases number of connected componenets
function bridge_edges(g::SimpleGraph)
    bridges = bridges(g)
    return DataFrame(src = [src(e) for e in bridges],
                     dst = [dst(e) for e in bridges])
end


# --- identify nodes with high betweenness and diverse community connections
function bridging_nodes(g::SimpleGraph, communities::Vector{Int}; top_n::Int=10)
    n = nv(g)
    bc = betweenness_centrality(g)
    
    # Participation coefficient (diversity of community connections)
    community_set = unique(communities)
    degree_per_comm = zeros(Float64, n, length(community_set))
    for e in edges(g)
        u, v = src(e), dst(e)
        degree_per_comm[u, communities[v]] += 1
        degree_per_comm[v, communities[u]] += 1
    end
    
    total_deg = degree(g)
    participation = [1 - sum((degree_per_comm[i, :] ./ max(total_deg[i], 1)).^2)
                     for i in 1:n]
    
    df = DataFrame(node = 1:n,
                   betweenness = bc,
                   participation = participation)
    df[!, :bridge_score] = zscore(df.betweenness) .+ zscore(df.participation)
    
    top_nodes = sort(df, :bridge_score, rev=true)[1:top_n, :]
    return top_nodes
end


# --- summary
function summarize_bridges(g::SimpleGraph, communities::Vector{Int}, targets_df::DataFrame; top_n::Int=10)
    println("\n================ STRUCTURAL BRIDGES SUMMARY ================\n")

    # --- Edge-level metrics ---
    edge_btw = edge_betweenness(g)
    bridge_df = bridge_edges(g)
    println("Edges in network:                 $(ne(g))")
    println("Bridge edges (disconnect graph):  $(nrow(bridge_df))")
    println("Avg. edge betweenness:            $(round(mean(edge_btw.betweenness), digits=4))")
    println("Max edge betweenness:             $(round(max(edge_btw.betweenness), digits=4))")

    # --- Node-level bridging metrics ---
    top_nodes = bridging_nodes(g, communities; top_n=top_n)

    # Join with metadata (page_name, page_type) if available
    if all(∈(["id", "page_name", "page_type"]), names(targets_df))
        meta_cols = [:id => :node, :page_name, :page_type]
        top_nodes = leftjoin(top_nodes, select(targets_df, meta_cols), on = :node => :id)
    elseif :id ∈ names(targets_df)
        top_nodes = leftjoin(top_nodes, select(targets_df, [:id, :page_name => :page_name, :page_type => :page_type], renamecols=false), on = :node => :id)
    end

    println("\nTop $(top_n) Bridging Nodes (by bridge score):")
    show(select(top_nodes, [:node, :page_name, :page_type, :betweenness, :participation, :bridge_score]); allrows=true, truncate=100)
    
    println("\nInterpretation:")
    println("• Bridge edges = links whose removal would fragment the graph.")
    println("• High betweenness + high participation ⇒ cross-community connectors.")
    println("• Pages with high bridge scores often include media outlets, organizations, or influencers connecting different interest groups.\n")

    return (edge_betweenness=edge_btw, bridge_edges=bridge_df, bridging_nodes=top_nodes)
end

export edge_betweenness, bridge_edges, bridging_nodes, summarize_bridges


end # module