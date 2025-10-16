module NetworkMetrics

using Graphs, StatsBase, Statistics, LinearAlgebra

# ---------------------------
# HOMOPHILY
# ---------------------------
# functions: categorical_mixing_matrix, categorical_assortativity, node_homophily, summarize_homophily

# --- Mixing matrix (normalized, sums to 1) ---
function categorical_mixing_matrix(g::Graphs.SimpleGraph, label_code::Vector{Int}, K::Int)
    M = zeros(Float64, K, K)
    for ed in Graphs.edges(g)                 # <-- fully qualified
        u = Graphs.src(ed); v = Graphs.dst(ed)  # <-- fully qualified
        lu = label_code[u]; lv = label_code[v]
        if lu == lv
            M[lu, lv] += 2.0
        else
            M[lu, lv] += 1.0
            M[lv, lu] += 1.0
        end
    end
    z = sum(M)
    return z > 0 ? M ./ z : M
end

# --- Newman assortativity (categorical attributes) ---
function categorical_assortativity(g::Graphs.SimpleGraph, label_code::Vector{Int}, K::Int)
    e = categorical_mixing_matrix(g, label_code, K)
    a = vec(sum(e, dims=2))
    tr_e = tr(e)
    a2 = sum(a .^ 2)
    denom = 1.0 - a2
    return iszero(denom) ? 0.0 : (tr_e - a2) / denom
end

# --- Node-level homophily (share of neighbors with same label) ---
function node_homophily(g::Graphs.SimpleGraph, label_code::Vector{Int})
    H = fill(NaN, Graphs.nv(g))
    for v in 1:Graphs.nv(g)
        nbrs = Graphs.neighbors(g, v)
        d = length(nbrs)
        if d > 0
            same = count(u -> label_code[u] == label_code[v], nbrs)
            H[v] = same / d
        end
    end
    return H
end

# --- Return summary of homophily ---
function summarize_homophily(g::SimpleGraph, label_code::Vector{Int}, K::Int; label_levels=nothing)
    e = categorical_mixing_matrix(g, label_code, K)
    edge_h = tr(e)                                   # observed same-category edge share
    freqs = counts(label_code, 1:K) ./ length(label_code)   # category frequencies p_k
    base = sum(freqs .^ 2)                                  # random-mixing baseline
    r = categorical_assortativity(g, label_code, K)         # assortativity coefficient
    H = node_homophily(g, label_code)
    Hnz = H[.!isnan.(H)]
    percat = [sum(e[i, i]) / sum(e[i, :]) for i in 1:K]     # per-category internal share (Vector{Float64} of length K)
    return (; edge_h, base, ratio = edge_h / base, r,
            H_mean = mean(Hnz), H_median = median(Hnz),     # node-level homophily stats
            percat, freqs)
end



export summarize_homophily
end # module
