module Describe

using CSV
using DataFrames
using Graphs

export analyze_network

"""
    analyze_network(edge_csv_path::String)

Reads a CSV file of edges and builds an undirected graph.
Prints number of nodes, edges, global clustering coefficient, and bridge statistics.
Returns the SimpleGraph object.
"""
function analyze_network(edge_csv_path::String)
    # Read the edge CSV
    edge = CSV.read(edge_csv_path, DataFrame)
    
    # Build graph
    total_nodes = unique(vcat(edge.id_1, edge.id_2))
    n = maximum(total_nodes) + 1
    g = SimpleGraph(n)
    for row in eachrow(edge)
        add_edge!(g, row.id_1 + 1, row.id_2 + 1)
    end
    
    # Print stats
    println("No. of Nodes in the Graph: ", nv(g))
    println("No. of Edges in the Graph: ", ne(g))
    println("Fraction of triangles that are closed: ", global_clustering_coefficient(g))
    
    # Bridges
    b = bridges(g)
    n_bridges = length(b)
    te = ne(g)
    percentage = n_bridges / te * 100
    println("Percent of edges that are bridges: ", percentage)
    println("No. of edges that form a bridge: ", n_bridges)
    
    return g
end

end
