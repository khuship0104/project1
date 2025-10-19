module Homophily

using Graphs
using LinearAlgebra, Statistics

export compute_homophily, chance_homophily, observed_homophily, category_homophily

# ------------------------------------------------------------
# 1. Chance homophily = Σ p_k²
#    (Baseline probability that two random nodes share a label)
# ------------------------------------------------------------
function chance_homophily(labels::Vector{Int}, K::Int)
    freqs = [count(==(k), labels) / length(labels) for k in 1:K]
    return sum(freqs .^ 2)
end

# ------------------------------------------------------------
# 2. Observed homophily
#    = fraction of edges connecting nodes with the same label
# ------------------------------------------------------------
function observed_homophily(g::SimpleGraph, labels::Vector{Int})
    total_edges = 0
    matching_edges = 0

    for e in edges(g)
        u, v = src(e), dst(e)
        total_edges += 1
        matching_edges += (labels[u] == labels[v]) ? 1 : 0
    end

    return total_edges == 0 ? 0.0 : matching_edges / total_edges
end

# ------------------------------------------------------------
# 3. Category-wise observed homophily
#    = fraction of edges involving a category that connect
#      to the same category
# ------------------------------------------------------------
function category_homophily(g::SimpleGraph, labels::Vector{Int}, category::Int)
    total_edges = 0
    matching_edges = 0

    for e in edges(g)
        u, v = src(e), dst(e)
        if labels[u] == category || labels[v] == category
            total_edges += 1
            matching_edges += (labels[u] == labels[v] && labels[u] == category) ? 1 : 0
        end
    end

    return total_edges == 0 ? 0.0 : matching_edges / total_edges
end

# ------------------------------------------------------------
# 4. Compute and display overall + category-wise results
# ------------------------------------------------------------
function compute_homophily(g::SimpleGraph, labels::Vector{Int}, category_map::Dict)
    K = length(category_map)

    # Overall metrics
    chance = chance_homophily(labels, K)
    observed = observed_homophily(g, labels)
    ratio = observed / chance

    println("\n🌐 Overall Homophily Results")
    println("-------------------------------------")
    println("Chance Homophily   : ", round(chance, digits=4))
    println("Observed Homophily : ", round(observed, digits=4))
    println("Ratio (Obs/Chance) : ", round(ratio, digits=2))
    println("-------------------------------------\n")

    # Category-wise metrics
    println("📊 Category-wise Homophily:")
    println("Category        | Observed | Chance | Ratio")
    println("----------------|----------|--------|-------")

    for (name, cid) in category_map
        obs = category_homophily(g, labels, cid)
        chance_c = (count(==(cid), labels) / nv(g))^2
        ratio_c = chance_c == 0 ? 0.0 : obs / chance_c

        println(rpad(name, 15), "| ",
                round(obs, digits=4), "  | ",
                round(chance_c, digits=4), " | ",
                round(ratio_c, digits=2))
    end

    return (; observed, chance, ratio)
end

end # module
