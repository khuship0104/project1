module Describe

using Graphs

export describe_graph

"""
    describe_graph(g::SimpleGraph)

Prints number of nodes, edges, clustering coefficient, and bridge statistics
for a given graph `g`.
"""
function describe_graph(g::SimpleGraph)
    println("No. of Nodes in the Graph: ", nv(g))
    println("No. of Edges in the Graph: ", ne(g))
    println("Fraction of triangles that are closed: ", global_clustering_coefficient(g))
    
    b = bridges(g)
    n_bridges = length(b)
    te = ne(g)
    percentage = n_bridges / te * 100
    println("Percent of edges that are bridges: ", percentage)
    println("No. of edges that form a bridge: ", n_bridges)
end

end