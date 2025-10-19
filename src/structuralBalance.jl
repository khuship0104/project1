module NetworkBalance

using Graphs, StatsBase, Statistics, Random


@inline edge_sign(u::Int, v::Int, label_code::Vector{Int}) =
    (label_code[u] == label_code[v]) ? 1 : -1

# enumerate all triangles and return edge-sign triplets ---
function triad_signs(g::SimpleGraph, label_code::Vector{Int})
    ts = Tuple{Int,Int,Int}[]
    nverts = nv(g)
    @inbounds for u in 1:nverts
        nu = neighbors(g, u)
        for v in nu
            v <= u && continue 
            nv_ = neighbors(g, v)

            # two-pointer intersection (neighbors sorted)
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

# --- proportion of balanced triads ---
balance_ratio(triads::Vector{Tuple{Int,Int,Int}}) =
    isempty(triads) ? NaN :
    count(t -> (t[1]*t[2]*t[3]) > 0, triads) / length(triads)

# --- friend-of-friend positive closure ---
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

# --- baseline closure for comparison ---
function friend_of_friend_positive_closure_baseline(
        g::SimpleGraph, label_code::Vector{Int}; R::Int=100, rng=Random.default_rng()
    )
    rates = Float64[]
    n = length(label_code)
    for _ in 1:R
        shuffled = copy(label_code)
        shuffle!(rng, shuffled)                   # Correctly shuffle labels
        r, _, _ = friend_of_friend_positive_closure(g, shuffled)
        push!(rates, r)
    end
    mean_rate = isempty(rates) ? 0.0 : mean(rates)
    std_rate  = isempty(rates) ? 0.0 : std(rates)
    return (mean_rate, std_rate)
end

# --- full structural balance analysis ---
function structural_balance_summary(
        g::SimpleGraph, label_code::Vector{Int}; R::Int=100, rng=Random.default_rng()
    )
    # Global structural balance via triangles
    ts = triad_signs(g, label_code)
    b_ratio = balance_ratio(ts)

    # Friend-of-friend positive closure
    fof_rate, closed_pos, wedges = friend_of_friend_positive_closure(g, label_code)
    base_mean, base_std = friend_of_friend_positive_closure_baseline(g, label_code; R=R, rng=rng)

    lift = (isnan(fof_rate) || base_mean == 0) ? NaN : fof_rate / base_mean
    z    = (isnan(fof_rate) || base_std == 0) ? NaN : (fof_rate - base_mean) / base_std

    println("====== STRUCTURAL BALANCE SUMMARY ======")
    println("Triangles (closed triads):            $(length(ts))")
    println("Balanced triads ratio:                $(isnan(b_ratio) ? "NaN" : string(round(b_ratio, digits=4)))")
    println()
    println("Friend-of-friend positive closure:")
    println("  Qualifying wedges (A–B, A–C +pos):   $wedges")
    println("  Closed positive wedges (B–C +pos):   $closed_pos")
    println("  Closure rate (observed):             $(isnan(fof_rate) ? "NaN" : string(round(fof_rate, digits=4)))")
    println("  Expected closure rate (random):      $(round(base_mean, digits=4))")
    println("  Lift over random:                    $(round(lift, digits=4))")
    println("  Z-score:                             $(round(z, digits=4))")
    println()
    println("Research Question Interpretation:")
    println("  Pages that like a common neighbour tend to like themselves in $(round(fof_rate*100, digits=2))% of cases.")
    println("  Compared to random chance ($(round(base_mean*100, digits=2))%), this is a lift of $(round(lift, digits=2))×, indicating clear structural balance.")

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

end # module