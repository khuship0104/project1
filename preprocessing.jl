module NetworkPreprocessing

using CSV, DataFrames, JSON3, Graphs
using SparseArrays, CategoricalArrays

# ---------------------------
# Load raw files
# ---------------------------
function load_raw(edges_path::AbstractString,
                  targets_path::AbstractString,
                  features_path::AbstractString)
    edges_df   = CSV.read(edges_path, DataFrame)
    targets_df = CSV.read(targets_path, DataFrame)
    features   = JSON3.read(read(features_path, String))
    return edges_df, targets_df, features
end

# ---------------------------
# Column Resolution
# ---------------------------

# Pick standardized source/target column name
function resolve_edge_columns(edges_df::DataFrame)
    src_col = hasproperty(edges_df, :id_1)  ? :id_1  :
              hasproperty(edges_df, :source) ? :source : names(edges_df)[1]
    dst_col = hasproperty(edges_df, :id_2)  ? :id_2  :
              hasproperty(edges_df, :target) ? :target : names(edges_df)[2]
    return src_col, dst_col
end

# Pick standardized id/label column names.
function resolve_target_columns(targets_df::DataFrame)
    id_col    = hasproperty(targets_df, :id) ? :id : names(targets_df)[1]
    label_col = hasproperty(targets_df, :page_type) ? :page_type :
                (hasproperty(targets_df, :target) ? :target : names(targets_df)[2])
    return id_col, label_col
end


# ---------------------------
# Node indexing and graph
# ---------------------------

# Creates a compact 1..N index for node ids found in edges.
function build_index(edges_df::DataFrame, src_col::Symbol, dst_col::Symbol)
    all_ids = unique(vcat(Vector{Int}(edges_df[:, src_col]), Vector{Int}(edges_df[:, dst_col])))
    sort!(all_ids)
    id2idx = Dict(id => i for (i, id) in enumerate(all_ids))
    return all_ids, id2idx
end

# Constructs an undirected SimpleGraph using the 1..N index.
function build_graph(edges_df::DataFrame, id2idx::Dict{Int,Int},
                     src_col::Symbol, dst_col::Symbol)
    N = length(id2idx)
    g = SimpleGraph(N)
    for r in eachrow(edges_df)
        u = id2idx[Int(r[src_col])]
        v = id2idx[Int(r[dst_col])]
        u != v && add_edge!(g, u, v)
    end
    return g
end

# ---------------------------
# Labels
# ---------------------------
function build_labels(targets_df::DataFrame, all_ids::Vector{Int}, id_col::Symbol, label_col::Symbol; unknown="Unknown")
    lab_map = Dict{Int,Union{Missing,String}}()
    for r in eachrow(targets_df)
        lab_map[Int(r[id_col])] = String(r[label_col])
    end
    labels_raw = [get(lab_map, id, missing) for id in all_ids]

    # Replace missings with "Unknown" 
    labels_str = [ismissing(x) ? unknown : x for x in labels_raw]
    labels = CategoricalArray(labels_str)

    # Convenience: numeric codes 1..K
    label_levels = levels(labels)
    label_code = Vector{Int}(labels.refs)

    return labels, label_code, label_levels
end

# ---------------------------
# Feature Matrix
# ---------------------------
function build_feature_matrix(features::Dict, all_ids::Vector{Int}, id2idx::Dict{Int,Int})
    node_keys = collect(keys(features))
    node_ids_in_feat = parse.(Int, string.(node_keys))  # normalize keys -> Int
    feat_lists = [features[k] for k in node_keys]

    max_feat = maximum(Iterators.flatten(feat_lists))
    is_zero_based = any(f -> f == 0, Iterators.flatten(feat_lists))
    feat_index = f -> is_zero_based ? (f .+ 1) : f
    F = (is_zero_based ? max_feat + 1 : max_feat)

    rowI = Int[]; colJ = Int[]; valV = Int[]
    id_set = Set(all_ids)

    for (nid, feats) in zip(node_ids_in_feat, feat_lists)
        if nid in id_set
            i = id2idx[nid]
            for f in feat_index(feats)
                push!(rowI, i); push!(colJ, f); push!(valV, 1)
            end
        end
    end

    X = sparse(rowI, colJ, valV, length(all_ids), F)

    return X, F, is_zero_based
end

# ---------------------------
# One-Shot convenience
# ---------------------------
function preprocess(edges_path::AbstractString,
                    targets_path::AbstractString,
                    features_path::AbstractString)
    edges_df, targets_df, features_raw = load_raw(edges_path, targets_path, features_path)
    features = Dict(string(k) => v for (k, v) in pairs(features_raw))

    src_col, dst_col = resolve_edge_columns(edges_df)
    id_col, label_col = resolve_target_columns(targets_df)

    all_ids, id2idx = build_index(edges_df, src_col, dst_col)
    g = build_graph(edges_df, id2idx, src_col, dst_col)

    labels, label_code, label_levels = build_labels(targets_df, all_ids, id_col, label_col)
    X, F, zero_based = build_feature_matrix(features, all_ids, id2idx)

    return (; g, labels, label_code, label_levels, X, id2idx, all_ids, F, zero_based)
end

export load_raw, resolve_edge_columns, resolve_target_columns,
       build_index, build_graph, build_labels, build_feature_matrix,
       preprocess

end # module