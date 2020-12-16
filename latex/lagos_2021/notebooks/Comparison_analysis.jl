# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: jl:percent
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Julia 1.4.2
#     language: julia
#     name: julia-1.4
# ---

# %%
# Ensure your working directory is: https://github.com/henriquebecker91/phd/tree/master/latex/revised_PPG2KP
include("notebook_setup.jl")

# %%
# Read the data, show nothing.
raw_cdf_path = "./data/lagos2021.csv"
raw_cdf = DataFrame(CSV.File(raw_cdf_path))
#nothing
showtable(raw_cdf)

# %%
# Read the data, show nothing.
raw_idf_path = "./data/dataset_c_info.csv"
raw_idf = DataFrame(CSV.File(raw_idf_path))
showtable(raw_idf)

# %%
# Clean the data a little
idf = let raw_idf = deepcopy(raw_idf)
    raw_idf[!, :instance] = basename.(raw_idf[!, :instance])
    select!(raw_idf, :instance => :instance_name, Not(:instance))
    raw_ne2df_path = "./data/dataset_c_2NE_opt.csv"
    raw_ne2df = DataFrame(CSV.File(raw_ne2df_path))
    raw_idf = innerjoin(raw_idf, raw_ne2df; on = :instance_name)

    raw_idf
end
showtable(idf)

# %%
# Clean the data: replace sentinel values by missing; add a string id combining the parameter combinations;
# create some boolean columns from the data.
cdf = let cdf = deepcopy(raw_cdf)
    # Keep only the instance name (not the path).
    @with(cdf, :instance_name .= basename.(:instance_name))
    @with(cdf, :datafile .= basename.(:datafile))

    cdf[!, :finished] = .!cdf[!, :run_ended_by_exception]

    colnames = names(cdf)
    for name in colnames
        column = cdf[!, name]
        if eltype(column) <: AbstractFloat
            cdf[!, name] = ifelse.(isnan.(column), missing, column)
        elseif eltype(column) <: Integer
            cdf[!, name] = ifelse.(column .< 0, missing, column)
        end
    end
    
    # Create columns for those attributes that indicate their final state before the final solving.
    for attr in (:cmvars, :pevars, :plates)
        @assert all(cdf[!, "qt_$(attr)_after_preprocess"] .== cdf[!, "qt_$(attr)_after_purge"])
        rename!(cdf, "qt_$(attr)_after_preprocess" => attr)
        select!(cdf, Not("qt_$(attr)_after_purge"))
    end
    
    # this is brittle, but it works because every run finished with a solution (even if not optimal)
    cdf[!, :g2csp] = (!ismissing).(cdf[!, :solution_value])
    cdf[!, :g2kp] = (!ismissing).(cdf[!, :solution_profit])
    @assert all(cdf[!, :g2csp] .⊻ cdf[!, :g2kp])
    
    cdf[!, :opt] = @with(cdf, ifelse.(
        :finished, ifelse.(:g2kp, :solution_profit, :solution_value), missing
    ))

    cdf
end
showtable(cdf)
#nothing

# %%
# Remove all gcut instances for now
filter!(:instance_name => s -> isnothing(match(r"^gcut", s)), cdf)

# %%
# Remove enhanced without PSN and original with PSN.
filter!([:faithful, :round2disc] => (f, psn) -> (f && !psn) || (!f  && psn), cdf)

# %%
# Work only with the multi-seed runs there.
let cdf = deepcopy(cdf)
    # Get only the runs with seeds greater than one, and then add the respectives seed one runs.
    extra = filter(:seed => !isequal(1), cdf)
    pk_no_seed = [:instance_name, :round2disc, :faithful, :threads, :g2csp]
    extra = semijoin(cdf, extra; on = pk_no_seed)
    pk = push!(copy(pk_no_seed), :seed)
    @assert all(extra[!, :threads] .== 12)
    @assert all(extra[!, :g2csp])
    summary = sort!(extra, pk) |> df -> select!(df, pk..., :build_and_solve_time) |>
        df -> groupby(df, [:instance_name, :faithful]) |>
        df -> @based_on(df, min_time = minimum(:build_and_solve_time), max_time = maximum(:build_and_solve_time))
    summary[!, :abs_gap] = @with(summary, :max_time .- :min_time)
end

# %%
# Remove all seeds different from one
filter!(:seed => isequal(1), cdf)

# %%
cidf = let cdf = deepcopy(cdf), idf = deepcopy(idf)
    function something_if_possible(itr)
        nmiss = skipmissing(itr)
        return isempty(nmiss) ? missing : first(nmiss)
    end
    
    select!(idf, :instance_name, :N, :n_, :L, :W, :min_l => :l_, :min_w => :w_,
        :kp_ub => :kp_tub, :csp_lb => :csp_tub, :opt2ne => :csp_opt2ne)
    csp_cdf = @linq filter(r -> r.g2csp, cdf) |> groupby(:instance_name) |>
        based_on(
            csp_opt = something_if_possible(:opt),
            csp_rel = mean(skipmissing(:final_root_relaxation_value))
        )
    kp_cdf = @linq filter(r -> r.g2kp, cdf) |> groupby(:instance_name) |>
        based_on(
            kp_opt = something_if_possible(:opt),
            kp_rel = mean(skipmissing(:final_root_relaxation_value))
        )
    
    cidf = innerjoin(idf, csp_cdf, kp_cdf; on = :instance_name)
    select!(cidf,
        :instance_name, :N, :n_, :L, :W, :l_, :w_,
        :kp_tub, :kp_rel, :kp_opt, :csp_tub, :csp_rel, :csp_opt, :csp_opt2ne
    )
    sort!(cidf, :instance_name)
end
showtable(cidf)

# %%
cidf[!, :kp_gap_tub_rel] = @with(cidf, :kp_tub ./ :kp_rel) .* 100 .- 100
sort(cidf, :kp_gap_tub_rel; rev = true) |> df -> select(df, :instance_name, :kp_gap_tub_rel) |> showtable
#mean(cidf[!, :kp_gap_tub_rel])
cidf[!, :kp_gap_rel_opt] = @with(cidf, :kp_rel ./ :kp_opt) .* 100 .- 100
#@show mean(cidf[!, :kp_gap_rel_opt])
sort(cidf, :kp_gap_rel_opt; rev = true) |> df -> select(df, :instance_name, :kp_gap_rel_opt) |> showtable

cidf[!, :csp_gap_tub_rel] = @with(cidf, ceil.(:csp_rel) .- ceil.(:csp_tub))
@show sum(iszero.(cidf[!, :csp_gap_tub_rel]))
sort(cidf, :csp_gap_tub_rel) |> df -> select(df, :instance_name, :csp_gap_tub_rel) |> showtable

#=filter(r -> r.instance_name in ("W", "A1", "3"), cidf) |>
    df -> select(df, :instance_name, :csp_gap_tub_rel, :csp_opt) |>
    df -> (df[!, :per_csp_gap_tub_rel] = df[!, :csp_gap_tub_rel] ./ df[!, :csp_opt]; df) |>
    showtable
=#
cidf[!, :csp_gap_rel_opt] = @with(cidf, :csp_opt .- ceil.(:csp_rel))
sort(cidf, :csp_gap_rel_opt) |> df -> select(df, :instance_name, :csp_gap_rel_opt) |> showtable

#cidf[!, :csp_gap_tub_opt] = @with(cidf, :csp_opt .- ceil.(:csp_tub))
#@show sum(iszero.(skipmissing(cidf[!, :csp_gap_tub_opt])))

# %%
latex_insts = let cidf = deepcopy(cidf)
    #:instance_name, :N, :n_, :L, :W, :kp_tub, :kp_rel, :kp_opt, :csp_tub, :csp_rel, :csp_opt
    cidf[!, :kp_tub] = trunc.(Int, cidf[!, :kp_tub]) # the decimal places are not relevant here
    cidf[!, :kp_rel] = trunc.(Int, cidf[!, :kp_rel]) # the decimal places are not relevant here
    
    for col in names(cidf)
        col == "instance_name" && continue
        cidf[!, col] = number2latex.(cidf[!, col]; enclose = false)
    end

    cidf
end

pretty_table(
    latex_insts; backend = :latex, nosubheader = true,
    alignment = vcat([:l], repeat([:r], ncol(latex_insts) - 1))
)

# %%
# Create a code for identifying the distinct runs for the same instance.

#= Version considering all combinations of PSN = #
function build_variant(threads, faithful, psn, problem)
    e, n = convert.(Int, (!faithful, psn))
    p = problem ? "k" : "c"
    return "$(threads)t_$(e)e_$(n)n_$(p)"
end

cdf[!, :variant] = build_variant.(
    getindex.((cdf,), (!,), [:threads, :faithful, :round2disc, :g2kp])...
)
# =#

#= Version considering only enhanced+PSN e faithful-PSN =#
function build_variant(threads, faithful, psn, problem)
    e, n = convert.(Int, (!faithful, psn))
    @assert e == n
    p = problem ? "k" : "c"
    return "$(threads)t_$(e)e_$(p)"
end

cdf[!, :variant] = build_variant.(
    getindex.((cdf,), (!,), [:threads, :faithful, :round2disc, :g2kp])...
)
# =#

# %%
select(cdf, :instance_name, :faithful, :threads, :g2kp, :g2csp, :total_instance_time, :rnr_solved_by_simplex) |>
    df -> filter(r -> r.rnr_solved_by_simplex && r.threads == 12, df) |>
    df -> select(df, :instance_name, :faithful, :g2kp, :g2csp, :total_instance_time) |>
    df -> sort(df, :instance_name) |>
    df -> groupby(df, :instance_name) |>
    df -> @based_on(df, faithful = sum(:faithful), g2kp = sum(:g2kp), g2csp = sum(:g2csp)) |>
    showtable
    #df -> filter(r -> r.rnr_solved_by_simplex, df)# |> 
    #df -> groupby(df, :instance_name) |>
    #df -> based_on(df, )

# %%
@assert all(ifelse.(cdf[!, :had_timeout], .!cdf[!, :finished], true)) 

# %%
# Check how many did not finish in each group.
@linq cdf |> groupby([:round2disc, :faithful, :threads, :g2kp]) |>
    based_on(; qt_total = length(:instance_name), qt_not_finished = sum(.!:finished)) |>
    sort([:threads, :faithful, :round2disc, :g2kp]) |>
    select!(:threads, :faithful, :round2disc, :g2kp, :qt_total, :qt_not_finished)

# %%
# Check which instances did not finish the "normal" way.
let cdf = deepcopy(cdf)
    atypical = @linq cdf |> where((.!:finished) .| :had_timeout) |>
        select(:instance_name, :threads, :faithful, :round2disc, :g2kp, :build_and_solve_time, :datafile) |>
        sort([:instance_name, :threads, :faithful, :round2disc, :g2kp])
    showtable(atypical)
end

# %%
# Just check if the optimals and the relaxations agree.
let cdf = deepcopy(cdf)
    problem = :g2csp # :g2kp #
    my_field = :solution_value # :solution_profit #
    relaxations = filter(row -> row[:finished] .& row[problem], cdf) |>
        df -> select(df, :instance_name, :variant, :final_root_relaxation_value) |>
        df -> unstack(df, :instance_name, :variant, :final_root_relaxation_value)
    optimals = filter(row -> row[:finished] .& row[problem], cdf) |>
        df -> select(df, :instance_name, :variant, my_field) |>
        df -> unstack(df, :instance_name, :variant, my_field)
    showtable(optimals)
    #showtable(relaxations)
    opt = map(eachrow(optimals)) do row
        instance = row[:instance_name]
        v = filter(!ismissing, collect(values(row[Not(:instance_name)])))
        #println(instance, " -- ", maximum(v) - minimum(v))
        @assert all(iszero(maximum(v) - minimum(v)))
        instance => first(v)
    end
    
    #=differences = map(eachrow(relaxations)) do row
        instance = row[:instance_name]
        v = filter(!ismissing, collect(values(row[Not(:instance_name)])))
        #println(instance, " -- ", maximum(v) - minimum(v))
        instance => maximum(v) - minimum(v)
    end
    rel_diff_by_opt = sort(first.(opt)) .=> (last.(sort(differences, by = first)) ./
        last.(sort(opt, by = first)))
    sort!(rel_diff_by_opt, by = last, rev = true)=#
end

# %%
latex_times = let cdf = deepcopy(cdf)
    cdf[!, :build_and_solve_time] = @with(cdf,
        ifelse.(:finished, :build_and_solve_time, missing)
    )
    times = filter(row -> row[:finished], cdf) |>
        df -> select(df, :instance_name, :variant, :build_and_solve_time) |>
        df -> unstack(df, :instance_name, :variant, :build_and_solve_time)
    rnr_simplex = cdf |>
        df -> select(df, :instance_name, :variant, :rnr_solved_by_simplex) |>
        df -> unstack(df, :instance_name, :variant, :rnr_solved_by_simplex)
    showtable(rnr_simplex)
    # Put columns in the desired order. rnr_simplex is refereed by column and
    # does not need that.
    col_order = [
        "$(t)t_$(e)e_$(p)" for p in ("k", "c") for e in (0, 1) for t in (1, 12)
    ]
    @show col_order
    select!(times, :instance_name, col_order...)

    for col in names(times)
        col == "instance_name" && continue
        times[!, col] = number2latex.(times[!, col]; enclose = false)
        times[!, col] = ifelse.(rnr_simplex[!, col],
            wrap_in_textit.(times[!, col]),
            times[!, col]
        )
    end
    showtable(times)

    highlight_best_values!(times; columns = 2:5, cleaner = dirt2number)
    highlight_best_values!(times; columns = 6:9, cleaner = dirt2number)
    
    sort!(times, :instance_name)
end

pretty_table(
    latex_times; backend = :latex, nosubheader = true, alignment = vcat([:l], repeat([:r], ncol(latex_times) - 1))
)

