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
#     display_name: Julia 1.5.4
#     language: julia
#     name: julia-1.5
# ---

# %%
# Ensure your working directory is: https://github.com/henriquebecker91/phd/tree/master/latex/revised_PPG2KP
include("notebook_setup.jl")

# %%
# Read the data, show nothing.
raw_cdf_path = "./data/comparison.csv"
raw_cdf = DataFrame(CSV.File(raw_cdf_path))
#nothing
showtable(raw_cdf)

# %%
# Clean the data: replace sentinel values by missing; replace bool parameters by a single string id;
# create columns with the relevant pevars/cmvars/plates for a given configuration.
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
        final_column = ifelse.(cdf[!, :pricing_method] .== "none",
            cdf[!, "qt_$(attr)_after_preprocess"],
            ifelse.(cdf[!, :purge_disabled],
                cdf[!, "qt_$(attr)_after_pricing"],
                cdf[!, "qt_$(attr)_after_purge"]
            )
        )
        cdf[!, attr] = final_column
    end

    primary_keys = ["pricing_method", "round2disc", "mip_start", "purge_disabled", "faithful"]
    @assert all(name -> name ∈ colnames, primary_keys)
    args2id = Dict{NTuple{5, Bool}, String}(
        (false, false, false, false, false) => "simple_revised",
        (false, true, false, false, false) => "rounded_revised",
        (false, true, true, false, false) => "warmed_rounded_revised",
        (true, true, true, false, false) => "priced_revised",
        (true, true, true, true, false) => "no_purge_priced_revised",
        (true, true, true, false, true) => "priced_faithful",
        (true, true, true, true, true) => "no_purge_priced_faithful",
        (false, false, false, false, true) => "faithful",
        (false, true, false, false, true) => "rounded_faithful",
        (false, true, true, false, true) => "warmed_rounded_faithful"
    )
    variant_column = getindex.((args2id,), tuple.(
        cdf[!, "pricing_method"] .== "furini",
        cdf[!, "round2disc"],
        (|).((cdf[!, "mip_start"] .== "guaranteed"),
            (&).(cdf[!, "mip_start"] .== "expected", cdf[!, "pricing_method"] .!= "none")),
        cdf[!, "purge_disabled"],
        cdf[!, "faithful"]
    ))
    #@show variant_column
    cdf[!, "model_variant"] = variant_column
    select!(cdf, Not(Symbol.(primary_keys)))
    
    cdf
end
showtable(cdf)
#nothing

# %%
# Should return an empty table. If the table have any lines they are runs that finished with
# the incorrect optimal value (so there is something wrong with the code).
let
    finished = select(cdf, :instance_name, :model_variant, :solution_profit, :finished)
    #finished = @linq filter!([:solution_profit, :finished] => ((p, f) -> f && !ismissing(p)), finished) |>
    #    select!(Not(:finished))
    finished = @linq filter!(:finished => identity, finished) |> select!(Not(:finished))
    finished = @linq finished |> groupby(:instance_name) |> based_on(;
        solution_profit = maximum(:solution_profit), qt_distinct = length(unique(:solution_profit))
    )
    @assert all(isone, finished[!, "qt_distinct"])
    finished = select!(finished, Not(:qt_distinct))
    opt59 = DataFrame(CSV.File("./data/opt59_by_thomo.csv"))
    joined = innerjoin(finished, opt59; on = :instance_name)
    wrong = filter(:solution_profit => isequal(:best), joined)
end

# %%
# Check how many did not finish in each group.
@linq cdf |> groupby(:model_variant) |>
    based_on(; qt_total = length(:model_variant), qt_not_finished = sum(.!:finished))

# %%
# Check which instances did not finish the "normal" way.
let cdf = deepcopy(cdf)
    atypical = @linq cdf |> where((.!:finished) .| (:build_stop_reason .!= "BUILT_MODEL") .| :had_timeout) |>
        select(:instance_name, :model_variant, :total_pricing_time, :build_stop_reason, :datafile) |>
        sort([:instance_name, :model_variant]) |>
        where(:model_variant .== "no_purge_priced_faithful")
    showtable(atypical)
end

# %%
# Check if the runs with missing pevars/cmvars/plates info are always a finished or FOUND_OPTIMUM run.
@linq cdf |> select(:datafile, :finished, :pevars, :cmvars, :plates, :build_stop_reason) |>
    where(ismissing.(:pevars) .| ismissing.(:cmvars) .| ismissing.(:plates))

# %%
# The data extracted from DOI 10.6092/unibo/amsdottorato/7399 by hand. The times of Priced PP-G2KP
# for the 59 instances used by Thomopulos in their 2016 thesis.
# The column 'total_time' has the time in seconds for the whole run, the column 'iterative_percentage'
# has which percentage of that time was spent in the iterative pricing, the column 'had_timeout'
# (boolean) indicates if some of the three time limits imposed to the run was broken.
th59_time_csv_path = "./data/thomo_ppg2kp_timefractions.csv"
th59_time = DataFrame(CSV.File(th59_time_csv_path))
th59_time[!, :had_timeout] = convert.(Bool, th59_time[!, :had_timeout])
th59_time[!, :iterative_time] = @with(th59_time, :total_time .* (:iterative_percent ./ 100.0))
@show names(th59_time)
@show typeof.(eachcol(th59_time))
showtable(th59_time)
#nothing

# %%
# Get the percentage of time spent inside iterative pricing by thomo and our method,
# using only the instances that had no timeout in any of them.
let th59_time = deepcopy(th59_time), cdf = deepcopy(cdf)
    select!(th59_time, :instance_name, :total_time, :iterative_time, :had_timeout)
    select!(cdf, :instance_name, :total_instance_time => :total_time,
        :iterated_pricing_time => :iterative_time, :had_timeout, :build_stop_reason, :model_variant
    )
    filter!(:model_variant => isequal("priced_faithful"), cdf)
    @assert nrow(cdf) == nrow(th59_time)
    for df in (th59_time, cdf)
        filter!(:had_timeout => !, df)
    end
    filter!(:build_stop_reason => isequal("BUILT_MODEL"), cdf)
    solved_on_both = intersect(th59_time[!, :instance_name], cdf[!, :instance_name])
    for df in (th59_time, cdf)
        filter!(:instance_name => (s -> s in solved_on_both), df)
    end
    @show solved_on_both
    @show length(solved_on_both)
    @assert sort(cdf[!, :instance_name]) == sort(th59_time[!, :instance_name])
    @show @with(cdf, sum(:iterative_time)/sum(:total_time))
    @show @with(th59_time, sum(:iterative_time)/sum(:total_time))
    @show sum(cdf[!, :total_time])/sum(th59_time[!, :total_time])
end

# %%
let cdf = deepcopy(cdf)
    temp = @linq cdf |> select(:datafile, :instance_name, :model_variant, :build_and_solve_time) |>
        where((:model_variant .== "rounded_revised") .| (:model_variant .== "warmed_rounded_revised")) |>
        sort([:instance_name, :model_variant])
    showtable(temp)
end


# %%
# Create the comparison time table
ctt = let cdf = deepcopy(cdf)
    # build_and_solve_time: measures only what we want to measure, but does not exist if unfinished
    # run_total_time: measures things like instance reading too, but exists even if unfinished
    used_columns = [
        :instance_name, :model_variant, :build_and_solve_time, :run_total_time, :finished,
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
