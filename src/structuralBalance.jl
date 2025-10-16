module NetworkBalance

using Graphs, StatsBase, Statistics, Random

# ---------------------------
# STRUCTURAL BALANCE
# ---------------------------
# functions: triad_signs, balance_ratio, friend_of_friend_positive_closure, structural_balance_summary

# --- enumerate all triangles and returns their edge-sign triplet ---
@inline edge_sign(u::Int, v::Int, label_code::Vector{Int}) =
    (label_code[u] == label_code[v]) ? 1 : -1

function triad_signs(g::SimpleGraph, label_code::Vector{Int})
    ts = Tuple{Int,Int,Int}[]
    nverts = nv(g)                     # <- renamed so we don't shadow Graphs.nv
    @inbounds for u in 1:nverts
        nu = neighbors(g, u)
        for v in nu
            v <= u && continue         # enforce u < v
            nv_ = neighbors(g, v)

            # two-pointer intersection (neighbors are sorted)
            i = 1; j = 1
            while i <= length(nu) && j <= length(nv_)
                wu = nu[i]; wv = nv_[j]
                if wu == v; i += 1; continue; end
                if wv == u; j += 1; continue; end
                if wu == wv
                    w = wu
                    w <= v && (i += 1; j += 1; continue)  # enforce v < w
                    push!(ts, (
                        (label_code[u] == label_code[v]) ? 1 : -1,
                        (label_code[v] == label_code[w]) ? 1 : -1,
                        (label_code[u] == label_code[w]) ? 1 : -1
                    ))
                    i += 1; j += 1
                elseif wu < wv
                    i += 1
                else
                    j += 1
                end
            end
        end
    end
    return ts
end




# --- Computes proportion of balanced triads (balanced if product of signs is positive) ---
balance_ratio(triads::Vector{Tuple{Int,Int,Int}}) =
    isempty(triads) ? NaN :
    count(t -> (t[1] * t[2] * t[3]) > 0, triads) / length(triads)


# Returns friend-of-friend closure
function friend_of_friend_positive_closure(g::SimpleGraph, label_code::Vector{Int})
    wedges = 0
    closed_pos = 0
    for a in vertices(g)
        nbrs = neighbors(g, a)
        for i in 1:length(nbrs)-1, j in i+1:length(nbrs)
            b, c = nbrs[i], nbrs[j]
            if label_code[a] == label_code[b] && label_code[a] == label_code[c]
                wedges += 1
                if has_edge(g, b, c) && (label_code[b] == label_code[c])
                    closed_pos += 1
                end
            end
        end
    end
    return wedges == 0 ? (NaN, 0, 0) : (closed_pos / wedges, closed_pos, wedges)
end

function friend_of_friend_positive_closure_baseline(
    g::SimpleGraph, label_code::Vector{Int}; R::Int=100, rng=Random.default_rng()
)
    rates = Float64[]
    for _ in 1:R
        shuffled = copy(label_code)
        randperm!(rng, shuffled)
        r, _, _ = friend_of_friend_positive_closure(g, shuffled)
        if !isnan(r); push!(rates, r); end
    end
    return isempty(rates) ? (NaN, NaN) : (mean(rates), std(rates))
end

# --- Runs full structural balance analysis
function structural_balance_summary(
    g::SimpleGraph, label_code::Vector{Int}; R::Int=100, rng=Random.default_rng()
)
    # Global structural balance via signed triangles
    ts = triad_signs(g, label_code)
    b_ratio = balance_ratio(ts)

    # Targeted diagnostic: “friends of a friend are friends”
    fof_rate, closed_pos, wedges = friend_of_friend_positive_closure(g, label_code)
    base_mean, base_std = friend_of_friend_positive_closure_baseline(g, label_code; R=R, rng=rng)

    lift = (isnan(fof_rate) || isnan(base_mean) || base_mean == 0) ? NaN : fof_rate / base_mean
    z = (isnan(fof_rate) || isnan(base_mean) || isnan(base_std) || base_std == 0) ? NaN : (fof_rate - base_mean) / base_std

    println("====== STRUCTURAL BALANCE SUMMARY ======")
    println("Triangles (closed triads):            $(length(ts))")
    println("Balanced triads ratio:                $(isnan(b_ratio) ? "NaN" : string(round(b_ratio, digits=4)))")
    println()
    println("Friend-of-friend positive closure:")
    println("  Qualifying wedges (A–B, A–C +pos):  $wedges")
    println("  Closed positive wedges (B–C +pos):  $closed_pos")
    println("  Closure rate (observed):            $(isnan(fof_rate) ? "NaN" : string(round(fof_rate, digits=4)))")
    println()

    return (
        n_triads = length(ts),
        balance_ratio = b_ratio,
        fof_pos_closure = fof_rate,
        fof_baseline_mean = base_mean,
        fof_baseline_std = base_std,
        fof_lift = lift,
        fof_zscore = z
    )
end

export structural_balance_summary
end    # module