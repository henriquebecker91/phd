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


# %%

# %%
# Create the comparison time table
ctt = let cdf = deepcopy(cdf)
    # build_and_solve_time: measures only what we want to measure, but does not exist if unfinished
    # run_total_time: measures things like instance reading too, but exists even if unfinished
    used_columns = [
        :instance_name, :, :build_and_solve_time, :run_total_time, :finished,
        :pevars, :cmvars, :plates, :build_stop_reason
    ]
    select!(cdf, used_columns)

    # We cannot remove the unfinished runs from the table, because things like
    # number of plates and variables of non-priced variants need to be summed even if
    # the run ended by timeout. This needs some jugglery with missing values.

    # Create a Bool column indicating if that row has the best time for all table.
    cdf[!, :was_best] = let cdf = deepcopy(cdf)
        replace!((v -> ismissing(v) ? Inf : v), cdf[!, :build_and_solve_time])
        b = @linq cdf |> select(:instance_name, :build_and_solve_time) |> groupby(:instance_name) |>
            based_on(; best_time = minimum(:build_and_solve_time))
        d = Dict{String, Float64}(map(=>, unique(cdf[!, :instance_name]), Iterators.cycle(Inf)))
        for row in eachrow(b)
            d[row.instance_name] = row.best_time
        end
        best_time_col = getindex.((d,), cdf[!, :instance_name])
        was_best_col = best_time_col .≈ cdf[!, :build_and_solve_time]
    end

    cdf[!, :build_and_solve_time] .= ifelse.(cdf[!, :finished], cdf[!, :build_and_solve_time], 0.0)
    cdf[!, :variables] = cdf[!, :cmvars] .+ cdf[!, :pevars]
    # Create the summary table. Transforms in zero the missing values that come from FOUND_OPTIMUM
    # (i.e., when the pricing procedure finds an optimal solution before generating a final model).
    # This depends on non-priced runs always enumerating the plates and variables
    cdf[!, :was_built] = @with cdf (.!occursin.(("priced",), :model_variant) .|
        (:build_stop_reason .== "BUILT_MODEL"))
    cdf[!, :variables] .= @with cdf ifelse.(:was_built, :variables, 0)
    cdf[!, :plates] .= @with cdf ifelse.(:was_built, :plates, 0)
    @assert iszero.(cdf[!, :variables]) == iszero.(cdf[!, :plates]) 
    cdf = @linq cdf |> groupby(:model_variant) |> based_on(;
        qt_solved = sum(:finished),
        solved_time = sum(:build_and_solve_time),
        qt_best = sum(:was_best),
        qt_model_built = sum(:was_built),
        qt_early = sum(:build_stop_reason .== "FOUND_OPTIMUM"),
        variables = sum(:variables),
        plates = sum(:plates)
    )
    cdf[!, :solved_time] = trunc.(Int, cdf[!, :solved_time])
    cdf[!, :total_time] = @with cdf (:solved_time .+ ((59 .- :qt_solved) * 3 * 60 * 60))
    never_finish_early = .!occursin.(("priced",), cdf[!, :model_variant])
    # If it never finishes early displays a dash instead of zero.
    cdf[!, :qt_early] = ifelse.(never_finish_early, missing, cdf[!, :qt_early])

    numeric_columns = [
        :total_time, :qt_solved, :qt_early, :solved_time, :qt_best, :qt_model_built, :variables, :plates
    ]
    for col in numeric_columns
        cdf[!, col] = number2latex.(cdf[!, col]; enclose = false)
    end

    # Rename the columns to the names in the table
    select!(cdf,
        :model_variant => "Variant",
        :total_time => "T. T.",
        :qt_early => "#e",
        :qt_model_built => "#m",
        :qt_solved => "#s",
        :qt_best => "#b",
        :solved_time => "S. T. T.",
        :variables => "#variables",
        :plates => "#plates",
    )
    select!(cdf, names(cdf) .=> esc_latex.(names(cdf)))
    # Rename the variants to the names in the table (also, define their order of appearance)
    pretty_variant_names = Dict{String, NamedTuple{(:pretty_name,:order), Tuple{String, Int}}}(
        "faithful" => (pretty_name = "Faithful", order = 1),
        "simple_revised" => (pretty_name = "Enhanced", order = 2),
        "rounded_faithful" => (pretty_name = "F. +Rounding", order = 3),
        "rounded_revised" => (pretty_name = "E. +Rounding", order = 4),
        "warmed_rounded_faithful" => (pretty_name = "F. +R. +Warming", order = 5),
        "warmed_rounded_revised" => (pretty_name = "E. +R. +Warming", order = 6),
        "priced_faithful" => (pretty_name = "Priced F. +R. +W.", order = 7),
        "priced_revised" => (pretty_name = "Priced E. +R. +W.", order = 8),
        "no_purge_priced_faithful" => (pretty_name = "P. F. +R. +W. -Purge", order = 9),
        "no_purge_priced_revised" => (pretty_name = "P. E. +R. +W. -Purge", order = 10)
    )
    sort!(cdf, "Variant"; by = (name -> pretty_variant_names[name].order))
    cdf[!, "Variant"] = getindex.(getindex.((pretty_variant_names,), cdf[!, "Variant"]), :pretty_name)

    cdf
end

pretty_table(
    ctt; backend = :latex, nosubheader = true, alignment = vcat([:l], repeat([:r], ncol(ctt) - 1))
)

#showtable(ctt)


# %%
# Create the comparison time table
ctt2 = let cdf = deepcopy(cdf)
    used_columns = [
        :model_variant, :enumeration_time, :heuristic_lb_time, :restricted_pricing_time,
        :iterated_pricing_time, :final_pricing_time, :final_solving_time, :build_and_solve_time,
        :final_root_relaxation_time, :finished, :build_stop_reason, :instance_name
    ]
    select!(cdf, used_columns)
    cdf = @linq cdf |> where(.!isnothing.(match.(r".*priced.*", :model_variant)))
    # We only consider runs that finished executing all phases (i.e., they cannot have hit
    # the time limit nor found an optimum early) in ALL variants.
    filter!([:finished, :build_stop_reason] => (f,bsr) -> f && bsr == "BUILT_MODEL", cdf)
    instances_to_consider = unique(cdf[!, :instance_name])
    for variant in unique(cdf[!, :model_variant])
        variant_set = @linq select(cdf, ["model_variant", "instance_name"]) |>
            where(:model_variant .== variant)
        instances_to_consider = intersect(
            instances_to_consider,
            unique(variant_set[!, :instance_name])
        )
    end
    filter!(:instance_name => (i -> i in instances_to_consider), cdf)

    select!(cdf, Not([:finished, :build_stop_reason, :instance_name]))
    
    # length(unique!(sort(raw_cdf[!, "instance_name"])))
    cdf = @linq cdf |> groupby(:model_variant) |> based_on(;
        qt_finished_all_phases = length(:model_variant),
        total_time = sum(:build_and_solve_time),
        enumeration_time = sum(:enumeration_time),
        restricted_heuristic_time = sum(:heuristic_lb_time),
        restricted_lp_and_mip_time = sum(:restricted_pricing_time) - sum(:heuristic_lb_time),
        iterated_pricing_time = sum(:iterated_pricing_time),
        final_pricing_time = sum(:final_pricing_time),
        final_solve_lp = sum(:final_root_relaxation_time),
        final_solve_beb = sum(:final_solving_time) - sum(:final_root_relaxation_time)
    )
    cdf[!, :other_time] = deepcopy(cdf[!, :total_time])
    for column in names(cdf) # compute other_time
        column in ("model_variant", "qt_finished_all_phases", "total_time", "other_time") && continue
        cdf[!, :other_time] .-= cdf[!, column]
    end
    for column in names(cdf) # transforms in percentages
        column in ("model_variant", "qt_finished_all_phases", "total_time") && continue
        cdf[!, column] ./= (cdf[!, :total_time] ./ 100.0) # gets percents, trust me
    end
    cdf[!, :total_time] = trunc.(Int, cdf[!, :total_time])

    select!(cdf,
        :model_variant => "Variant",
        :qt_finished_all_phases => "#",
        :total_time => "Time",
        :enumeration_time => "E",
        :restricted_heuristic_time => "RH",
        :restricted_lp_and_mip_time => "RP",
        :iterated_pricing_time => "IP",
        :final_pricing_time => "FP",
        :final_solve_lp => "LP",
        :final_solve_beb => "BB",
        #:other_time => "Ot."
    )
    # All rows should work over the same subset of instances.
    @assert isone(length(unique(cdf[!, "#"])))
    select!(cdf, Not("#"))

    #@show names(cdf)
    #@show unique(cdf[!, "Variant"])
    pretty_variant_names = Dict{String, NamedTuple{(:pretty_name,:order), Tuple{String, Int}}}(
        "priced_faithful" => (pretty_name = "Priced Faithful +R. +W.", order = 1),
        "priced_revised" => (pretty_name = "Priced Enhanced +R. +W.", order = 2),
        "no_purge_priced_faithful" => (pretty_name = "P. F. +R. +W. -Purge", order = 3),
        "no_purge_priced_revised" => (pretty_name = "P. E. +R. +W. -Purge", order = 4),
    )
    sort!(cdf, "Variant"; by = (name -> pretty_variant_names[name].order))
    cdf[!, "Variant"] = getindex.(getindex.((pretty_variant_names,), cdf[!, "Variant"]), :pretty_name)

    for column in names(cdf)
        column == "Variant" && continue
        cdf[!, column] = number2latex.(cdf[!, column]; enclose = false)
    end
    select!(cdf, names(cdf) .=> esc_latex.(names(cdf)))

    cdf
end

pretty_table(
    ctt2; backend = :latex, nosubheader = true, alignment = vcat([:l], repeat([:r], ncol(ctt2) - 1))
)

showtable(ctt2)

# %%
