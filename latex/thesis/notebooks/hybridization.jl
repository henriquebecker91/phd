# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: jl:percent
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.5.4
#     language: julia
#     name: julia-1.5
# ---

# %%
# Ensure your working directory is: https://github.com/henriquebecker91/phd/tree/master/latex/revised_PPG2KP
include("notebook_setup.jl")

# %%
function only_value(xs)
    x = first(xs)
    @assert all(==(x), xs)
    x
end

# %%
# Read the data, show nothing.
raw_cdf_path = "./data/G2KP_hybridization.csv"
raw_cdf = DataFrame(CSV.File(raw_cdf_path))
#nothing
showtable(raw_cdf)

# %%
# No run was done with the faithful algorithm.
@assert all(iszero, raw_cdf.faithful)
# No pricing was used.
@assert all(==("none"), raw_cdf.pricing)
# All runs used round2disc
@assert all(isone, raw_cdf.round2disc)
# No rotation runs.
@assert all(iszero, raw_cdf.rotation .| raw_cdf.mirror_plates)
# Every instance appear exactly three times.
@assert all(==(3), count.(.==(unique(raw_cdf.instance_path)), (raw_cdf.instance_path,)))
# No interrupted instances.
@assert all(!isnan, raw_cdf.build_and_solve_time)
@assert all(==("BUILT_MODEL"), raw_cdf.build_stop_reason)
# The value of the objective function and of the extracted solution match. 
@assert all(isapprox(replace(raw_cdf.obj_value, NaN => -1.0), raw_cdf.solution_value))
@assert all(
    (((h, qt_h),) -> iszero(h) ? iszero(qt_h) : true),
    zip(raw_cdf.hybridization, raw_cdf.hybridizations)
)
@assert all(
    (((ah, qt_extra_vars),) -> iszero(ah) ? iszero(qt_extra_vars) : true),
    zip(raw_cdf.aggressive_hybridization, raw_cdf.aggressive_hybridization_extra_vars)
)
@assert all(
    (((h, fb_killed),) -> iszero(h) ? iszero(fb_killed) : true),
    zip(raw_cdf.hybridization, raw_cdf.first_borns_killed_by_hybridization)
)
@show (sort!(filter(:hybridization => iszero, raw_cdf), :instance_path).num_vars .+
    sort!(filter(:aggressive_hybridization => isone, raw_cdf), :instance_path).aggressive_hybridization_extra_vars) .>=
    sort!(filter(:aggressive_hybridization => isone, raw_cdf), :instance_path).num_vars


# %%
function flags2variant!(df)
    primary_keys = ["hybridization", "aggressive_hybridization"]

    @assert all(name -> name in names(df), primary_keys)
    args2id = Dict{NTuple{2, Bool}, String}(
        (false, false) => "no_hybridization",
        (true,  false) => "hybridization",
        (true,   true) => "agg_hybridization",
    )
    variant_column = getindex.((args2id,), tuple.(
        isone.(df[!, :hybridization]), isone.(df[!, :aggressive_hybridization])
    ))
    #@show variant_column
    df[!, "model_variant"] = variant_column
    select!(df, Not(Symbol.(primary_keys)))
    df
end

function get_clean_data(df)
    cdf = deepcopy(df)
    # Keep only the instance name (not the path).
    cdf[!, :instance_name] = basename.(cdf[!, :instance_path])
    cdf[!, :datafile] = basename.(cdf[!, :this_data_file])

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
    @assert all(==("none"), cdf[!, :pricing])
    for attr in (:cmvars, :pevars, :plates)
        cdf[!, attr] = cdf[!, "qt_$(attr)_after_preprocess"]
    end

    flags2variant!(cdf)

    cdf
end

# Clean the data: replace sentinel values by missing; replace bool parameters by a single string id;
# create columns with the relevant pevars/cmvars/plates for a given configuration.
cdf = get_clean_data(raw_cdf)
showtable(cdf)
#nothing

# %%
# Just check if all variants agree on the objective
let df = deepcopy(cdf)
    combine(groupby(cdf, :instance_name), :obj_value => only_value => :obj_value)
end

# %%
# Read the data, show nothing.
raw_fmt59_stats_path = "./data/G2KP_FMT59_instance_stats.csv"
raw_fmt59_stats = DataFrame(CSV.File(raw_fmt59_stats_path))
#nothing
@show eltype.(identity.(eachcol(raw_fmt59_stats)))
showtable(raw_fmt59_stats)

# Clean the data: replace sentinel values by missing; replace bool parameters by a single string id;
# create columns with the relevant pevars/cmvars/plates for a given configuration.
raw_fmt59 = let df = deepcopy(raw_fmt59_stats)
    for col in [:rrl, :rrw]
        df[!, col] = parse.((Rational{BigInt},), df[!, col])
    end
    df
end

let df = deepcopy(raw_fmt59)
    @show mean(df[!, :rrw])
    @show mean(df[!, :rrl])
    @show Float64(mean(df[!, :rrw]))
    @show Float64(mean(df[!, :rrl]))
end

# %%
rnr_gap_df = let df = deepcopy(cdf)
    select!(df, :instance_name, :model_variant, :root_node_value)
    udf = unstack(df, :model_variant, :root_node_value)
    fudf = filter(
        [:no_hybridization, :hybridization, :agg_hybridization] => (n, c, a) -> abs(c/n - 1.0) >= 10e-5 || abs(a/n - 1.0) > 10e-5,
        udf
    )
    fudf[!, :hybridization_percent] = ((fudf[!, :hybridization] ./ fudf[!, :no_hybridization]) .- 1.0) * 100.0
    fudf[!, :agg_hybridization_percent] = ((fudf[!, :agg_hybridization] ./ fudf[!, :no_hybridization]) .- 1.0) * 100.0
    # Order columns
    select!(
        fudf, [:instance_name, :no_hybridization, :hybridization, :hybridization_percent,
            :agg_hybridization, :agg_hybridization_percent]
    )
    #showtable(fudf)
end

showtable(rnr_gap_df)

# %%
let df = deepcopy(rnr_gap_df)
    for col in names(df)
        col == "instance_name" && continue
        df[!, col] = number2latex.(df[!, col]; enclose = false)
    end
    pretty_table(
        df; backend = :latex, nosubheader = true,
        alignment = vcat([:l], repeat([:r], ncol(df) - 1))
    )
end

# %%
let df = deepcopy(cdf)
    df[!, :extra_vars_per_hyb] = df[!, :aggressive_hybridization_extra_vars] ./ df[!, :hybridizations]
    showtable(sort!(
        df, :extra_vars_per_hyb
    ))
end

# %%
# Should return an empty table. If the table have any lines they are runs that finished with
# the incorrect optimal value (so there is something wrong with the code).
let
    finished = select(cdf, :instance_name, :model_variant, :solution_value, :finished)
    #finished = @linq filter!([:solution_value, :finished] => ((p, f) -> f && !ismissing(p)), finished) |>
    #    select!(Not(:finished))
    finished = @linq filter!(:finished => identity, finished) |> select!(Not(:finished))
    finished = @linq finished |> groupby(:instance_name) |> based_on(;
        solution_value = maximum(:solution_value), qt_distinct = length(unique(:solution_value))
    )
    @assert all(isone, finished[!, "qt_distinct"])
    finished = select!(finished, Not(:qt_distinct))
    opt59 = DataFrame(CSV.File("./data/opt59_by_thomo.csv"))
    joined = innerjoin(finished, opt59; on = :instance_name)
    wrong = filter(:solution_value => isequal(:best), joined)
end

# %%
# Check how many did not finish in each group.
@linq cdf |> groupby(:model_variant) |>
    based_on(; qt_total = length(:model_variant), qt_not_finished = sum(.!:finished))

# %%
# Check if the runs with missing pevars/cmvars/plates info are always a finished or FOUND_OPTIMUM run.
@linq cdf |> select(:datafile, :finished, :pevars, :cmvars, :plates, :build_stop_reason) |>
    where(ismissing.(:pevars) .| ismissing.(:cmvars) .| ismissing.(:plates))

# %%
# Create the comparison time table
ctt = let cdf = deepcopy(cdf)
    # build_and_solve_time: measures only what we want to measure, but does not exist if unfinished
    # run_total_time: measures things like instance reading too, but exists even if unfinished
    used_columns = [
        :instance_name, :model_variant, :build_and_solve_time, :run_total_time, :finished,
        :pevars, :cmvars, :plates, :build_stop_reason, :hybridizations,
        :first_borns_killed_by_hybridization
    ]
    select!(cdf, used_columns)

    # We cannot remove the unfinished runs from the table, because things like
    # number of plates and variables of non-priced variants need to be summed even if
    # the run ended by timeout. This needs some jugglery with missing values.

    # Create a Float64 column indicating the best time for that instance in the table and
    # a Bool column indicating if the row has the best time for all table.
    cdf[!, :best_time_for_inst] = let cdf = deepcopy(cdf)
        replace!((v -> ismissing(v) ? Inf : v), cdf[!, :build_and_solve_time])
        b = @linq cdf |> select(:instance_name, :build_and_solve_time) |> groupby(:instance_name) |>
            based_on(; best_time = minimum(:build_and_solve_time))
        d = Dict{String, Float64}(map(=>, unique(cdf[!, :instance_name]), Iterators.cycle(Inf)))
        for row in eachrow(b)
            d[row.instance_name] = row.best_time
        end
        best_time_col = getindex.((d,), cdf[!, :instance_name])
    end

    cdf[!, :was_best] = cdf[!, :best_time_for_inst] .≈ cdf[!, :build_and_solve_time]
    cdf[!, :loss_on_non_best] = cdf[!, :build_and_solve_time] .- cdf[!, :best_time_for_inst]
    #sort!(select(cdf, :instance_name, :model_variant, :loss_on_non_best), [:instance_name, :model_variant])

    cdf[!, :build_and_solve_time] .= ifelse.(cdf[!, :finished], cdf[!, :build_and_solve_time], 0.0)
    cdf[!, :variables] = cdf[!, :cmvars] .+ cdf[!, :pevars]

    @assert iszero.(cdf[!, :variables]) == iszero.(cdf[!, :plates]) 
    cdf = @linq cdf |> groupby(:model_variant) |> based_on(;
        solved_time = sum(:build_and_solve_time),
        qt_best = sum(:was_best),
        time_loss = sum(:loss_on_non_best),
        variables = sum(:variables),
        plates = sum(:plates),
        qt_hybrid = sum(:hybridizations),
        qt_killed_fc = sum(:first_borns_killed_by_hybridization),
    )
    cdf[!, :per_hybrid] = (cdf[!, :qt_hybrid] ./ cdf[!, :variables]) .* 100
    cdf[!, :per_killed_fc] = (cdf[!, :qt_killed_fc] ./ cdf[!, :variables]) .* 100

    numeric_columns = [
        :solved_time, :time_loss, :qt_best, :variables, :plates,
        :per_hybrid, :per_killed_fc
    ]
    for col in numeric_columns
        cdf[!, col] = number2latex.(round.((Int,), cdf[!, col]); enclose = false)
    end
    
    # Rename the columns to the names in the table
    select!(cdf,
        :model_variant => "Variant",
        :solved_time => "T. T.",
        :time_loss => "T. L.",
        :qt_best => "#b",
        :variables => "#variables",
        :per_hybrid => "h %",
        :per_killed_fc => "k %",
        :plates => "#plates",
    )
    select!(cdf, names(cdf) .=> esc_latex.(names(cdf)))
    # Rename the variants to the names in the table (also, define their order of appearance)
    pretty_variant_names = Dict{String, NamedTuple{(:pretty_name,:order), Tuple{String, Int}}}(
        "no_hybridization" => (pretty_name = "N. H.", order = 1),
        "hybridization" => (pretty_name = "C. H.", order = 2),
        "agg_hybridization" => (pretty_name = "A. H.", order = 3),
    )
    sort!(cdf, "Variant"; by = (name -> pretty_variant_names[name].order))
    cdf[!, "Variant"] = getindex.(getindex.((pretty_variant_names,), cdf[!, "Variant"]), :pretty_name)

    cdf
end

# #=
pretty_table(
    ctt; backend = :latex, nosubheader = true, alignment = vcat([:l], repeat([:r], ncol(ctt) - 1))
)
# =#

#showtable(ctt)


# %%
# Create plot of the three variants, y axis is time, x axis is the instance (ordered by?).
hyb_bar_plot = let df = deepcopy(cdf)
    select!(df, :instance_name, :model_variant, :build_and_solve_time, :hybridizations, :num_vars)
    df[!, :worst_time_for_inst] = let cdf = deepcopy(df)
        b = @linq cdf |> select(:instance_name, :build_and_solve_time) |> groupby(:instance_name) |>
            based_on(; worst_time = maximum(:build_and_solve_time))
        d = Dict{String, Float64}(map(=>, unique(cdf[!, :instance_name]), Iterators.cycle(Inf)))
        for row in eachrow(b)
            d[row.instance_name] = row.worst_time
        end
        getindex.((d,), cdf[!, :instance_name])
    end

    df[!, :per_hybrid] = number2latex.(
        trunc.((Int,), df[!, :hybridizations] ./ df[!, :num_vars] .* 100); enclose = false
    )
    suffix_df = combine(
        groupby(sort!(
            select(df, :instance_name, :model_variant, :per_hybrid),
        [:instance_name, :model_variant]; rev = true), :instance_name),
    :per_hybrid => (s -> " (" * join(s[2:3], ";") * ")") => :suffix)

    df[!, :ylabel] = df[!, :instance_name] .* suffix_df[
        findfirst.(.==(df[!, :instance_name]), (suffix_df[!, :instance_name],)),
    :suffix]

    filter!(:worst_time_for_inst => >(10.0), df)
    sort!(df, :worst_time_for_inst)
    df[!, "Variant"] = replace(
        df[!, "model_variant"], "agg_hybridization" => "A. H.",
        "hybridization" => "C. H.", "no_hybridization" => "N. H."
    )
    
    plot(
        df, x = "ylabel", y = "build_and_solve_time", color="Variant",
        Geom.bar(; position = :dodge),# Scale.y_log10,
        Guide.xlabel("Instance Name (h % for C. H. and A. H., resp.)"),
        Guide.ylabel("Time to solve (seconds)"),
        Guide.yticks(; ticks = [1000, 2000, 3000, 3600])
        #Gadfly.Coord.cartesian(; ymax = 3600)
        #Guide.title("Run times (only instances for which worst run took 10s or more)")
    )
end

#hyb_bar_plot |> PDF("../plots/hyb_bar_plot.pdf")

# %%
# Read the data, show nothing.
raw_edf_path = "./data/G2KP_hybridization_extra.csv"
raw_edf = DataFrame(CSV.File(raw_edf_path))
    #nothing
showtable(raw_edf)

# %%
# Combine extra runs with the runs with the first seed 
edf = let odf = deepcopy(raw_cdf), edf = deepcopy(raw_edf)
    odf[!, :solver_seed] = ones(nrow(odf))
    @assert all(in(2:10), edf[!, :solver_seed])
    filter!(:instance_path => in(edf[!, :instance_path]), odf)
    df = vcat(odf, edf)
    get_clean_data(df)
end

showtable(edf)

# %%
# Plot percentage of barrier time to time
per_barrier_time_extra_edf = let df = deepcopy(edf)
    df[!, :stop_reason] = last.(split.(df[!, :stop_reason], ('.',)))
    @assert all(==("OPTIMAL"), df[!, :stop_reason])
    @assert all(!isnan, df[!, :root_node_time])
    # Fix the precision problem with Gurobi output
    df[!, :root_node_time] = ifelse.(df[!, :root_node_time] .<= 0.01, df[!, :build_and_solve_time], df[!, :root_node_time])
    
    df[!, :percent_barrier] = (df[!, :root_node_time] ./ df[!, :build_and_solve_time]) .* 100
    
    #flatten(
    df = transform(groupby(df, [:instance_name, :model_variant]),
        (:build_and_solve_time => (t -> std(t)/mean(t)) => :cv))
    #, setdiff(names(df), [:instance_name, :instance_path, :model_variant]))
    
    @assert all(==("gcut1"), filter(:cv => >(1.0), df).instance_name)
    
    df
end
@show filter(:cv => >(1.0), select(per_barrier_time_extra_edf, :instance_name, :root_node_time, :build_and_solve_time, :percent_barrier, :cv))
plot(per_barrier_time_extra_edf, x = "percent_barrier", color="stop_reason", Geom.histogram(;bincount = 10))
per_barrier_time_extra_edf[!, :cv] = string.(per_barrier_time_extra_edf[!, :cv])
sort!(per_barrier_time_extra_edf, :cv)
plot(filter(:root_node_time => >(10.0), per_barrier_time_extra_edf), y = "percent_barrier", x = "cv", color = "instance_name", Geom.boxplot)

# %%
# Create plot of the three variants, y axis is time, x axis is the instance (ordered by?).
hyb_box_plot = let df = deepcopy(edf)
    select!(df, :instance_name, :model_variant, :build_and_solve_time, :hybridizations, :num_vars, :solver_seed)
    df[!, :worst_time_for_inst] = let cdf = deepcopy(df)
        b = @linq cdf |> select(:instance_name, :build_and_solve_time) |> groupby(:instance_name) |>
            based_on(; worst_time = maximum(:build_and_solve_time))
        d = Dict{String, Float64}(map(=>, unique(cdf[!, :instance_name]), Iterators.cycle(Inf)))
        for row in eachrow(b)
            d[row.instance_name] = row.worst_time
        end
        getindex.((d,), cdf[!, :instance_name])
    end

    # Gets the number of hybridizations for each combination of instance/variant
    df[!, :per_hybrid] = number2latex.(
        trunc.((Int,), df[!, :hybridizations] ./ df[!, :num_vars] .* 100); enclose = false
    )
    suffix_df = combine(
        groupby(sort!(
            select(filter(:solver_seed => isone, df), :instance_name, :model_variant, :per_hybrid),
        [:instance_name, :model_variant]; rev = true), :instance_name),
    :per_hybrid => (s -> " (" * join(s[2:3], ";") * ")") => :suffix)

    df[!, :ylabel] = df[!, :instance_name] .* suffix_df[
        findfirst.(.==(df[!, :instance_name]), (suffix_df[!, :instance_name],)),
    :suffix]

    filter!(:worst_time_for_inst => >(100.0), df)
    sort!(df, :worst_time_for_inst)
    df[!, "Variant"] = replace(
        df[!, "model_variant"], "agg_hybridization" => "A. H.",
        "hybridization" => "C. H.", "no_hybridization" => "N. H."
    )
    
    plot(
        df, x = "ylabel", y = "build_and_solve_time", color="Variant",
        Geom.boxplot(),
        #Geom.bar(; position = :dodge),
        #Scale.y_log10,
        Guide.xlabel("Instance Name (h % for C. H. and A. H., resp.)"),
        Guide.ylabel("Time to solve (seconds)"),
        #Guide.yticks(; ticks = [1000, 2000, 3000, 3600])
        #Gadfly.Coord.cartesian(; ymax = 3600)
        #Guide.title("Run times (only instances for which worst run took 10s or more)")
    )
end

hyb_box_plot |> PDF("../plots/hyb_box_plot.pdf")
hyb_box_plot

# %%
let df = deepcopy(edf)
    select!(df, :instance_name, :solver_seed, :model_variant, :root_node_value)
    udf = unstack(df, :model_variant, :root_node_value)
    fudf = filter(
        [:no_hybridization, :hybridization, :agg_hybridization] => (n, c, a) -> !isapprox(n, c) || !isapprox(c, a),
        udf
    )
    sort!(fudf, [:instance_name, :solver_seed])
    showtable(fudf)
end

# %%
let df = deepcopy(edf)
    replace!(
        df[!, :model_variant],
        "no_hybridization" => "Base", "hybridization" => "Hyb.",
        "agg_hybridization" => "Agg. H."
    )
    df[!, :order] = replace(df[!, :model_variant],
        "Base" => 1, "Hyb." => 2, "Agg. H." => 3
    )
    #df[!, :discrete_seed] = string.(df[!, :solver_seed])
    sort!(df, [:order, :solver_seed])
    
    # #=
    ps = map(instance -> plot(
        filter(:instance_name => ==(instance), df), Guide.title(instance),
        Guide.ylabel("Time (s)"; orientation = :vertical),
        Guide.xlabel("Variant"; orientation = :horizontal),
        Scale.y_continuous(; format = :plain),
        #color = "discrete_seed",
        x = "model_variant", y = "build_and_solve_time", Geom.point
    ), unique(df[!, :instance_name]))
    gp = gridstack(reshape(ps, (2, 2)))
    gp |> PS("../plots/hyb_hard4_seeded.ps")
    # =#
    #=
    plot(
        df, xgroup = "instance_name", x = "build_and_solve_time", y = "model_variant",
        Geom.subplot_grid(Geom.point)
    )
    =#
end

# %%
extra_ussdf = let df = deepcopy(edf)
    #=
    replace!(
        df[!, :model_variant],
        "no_hybridization" => "Base", "hybridization" => "Hyb.",
        "agg_hybridization" => "Agg. H."
    )
    =#
    
    sdf = combine(
        groupby(df, [:instance_name, :model_variant]),
        :build_and_solve_time => length => :nrow,
        :build_and_solve_time => mean => :mean,
        :build_and_solve_time => std  => :std,
        :build_and_solve_time => minimum => :min,
        :build_and_solve_time => maximum => :max,
        [:build_and_solve_time, :root_node_time] => ((b, r) -> (sum(r)/sum(b)) * 100) => :per_avg_barrier_time,
        [:build_and_solve_time, :root_node_time] => ((b, r) -> (1.0 - (sum(r)/sum(b))) * 100) => :per_non_root_time,
        [:hybridizations, :pevars, :cmvars] => ((h, p, e) -> 100.0 * (sum(h) / (sum(p) + sum(e)))) => :h_per,
        [:first_borns_killed_by_hybridization, :pevars, :cmvars] => ((k, p, e) -> 100.0 * (sum(k) / (sum(p) + sum(e)))) => :k_per,
    )
    @assert all(==(10), sdf[!, :nrow])
    select!(sdf, Not(:nrow))
    sdf[!, :variant_order] = replace(sdf[!, :model_variant],
        "no_hybridization" => 1, "hybridization" => 2, "agg_hybridization" => 3
    )
    stacked_cols = [:h_per, :k_per, :mean, :std, :min, :max, :per_avg_barrier_time, :per_non_root_time]
    ssdf = stack(
        sdf, stacked_cols;
        variable_eltype = Symbol
    )
    ssdf[!, :stat_order] = replace(ssdf[!, :variable],
        (stacked_cols .=> 1:length(stacked_cols))...
    )
    ssdf[!, :order] = @. (ssdf[!, :stat_order] * 10) + ssdf[!, :variant_order]
    @. ssdf[!, :variable] = string(ssdf[!, :order]) * " " * ssdf[!, :model_variant] *
        " " * string(ssdf[!, :variable])
    select!(ssdf, Not([:model_variant, :variant_order, :stat_order, :order]))
    
    ussdf = unstack(ssdf)
    @show names(ussdf)
    # COMPUTE PERCENT COLUMNS
    non_key_names = filter(!=("instance_name"), names(ussdf))
    rename!(ussdf, non_key_names .=> join.(getindex.(split.(non_key_names, " "), ([2, 3],)), " "))
    # #=
    for variant in ("hybridization", "agg_hybridization")
        for stat in ("mean", #=, "std"=#)
            ussdf[!, "percent $variant $stat"] = (
                ((ussdf[!, "$variant $stat"] ./ ussdf[!, "no_hybridization $stat"]) #=.- 1.0=#) .* 100
            )
        end
    end
    # Changes standard deviation to coefficient of deviation
    for variant in ("no_hybridization", "hybridization", "agg_hybridization")
        ussdf[!, "$variant cv"] = (
            ((ussdf[!, "$variant std"] ./ ussdf[!, "$variant mean"]) #=.- 1.0=#) .* 100
        )
    end
    # =#
    
    #showtable(ussdf)
    ussdf
end

# %%
let df = deepcopy(extra_ussdf)
    # Remove min and max columns
    #select!(df, filter(s -> !occursin("min", s) && !occursin("max", s), names(df)))
    times = collect(zip(
        eachcol(select(df, ["no_hybridization mean", "hybridization mean", "agg_hybridization mean"]))...
    ))
    @show names(eachcol(select(df, ["no_hybridization mean", "hybridization mean", "agg_hybridization mean"])))
    df[!, :best_variant] = last.(findmin.(times))

    select!(df,
        :instance_name, "hybridization h_per", "agg_hybridization h_per",
        "no_hybridization mean", "percent hybridization mean", "percent agg_hybridization mean",
        #="no_hybridization std", "percent hybridization std", "percent agg_hybridization std",=#
        "no_hybridization cv", "hybridization cv", "agg_hybridization cv",
        "no_hybridization per_non_root_time", "hybridization per_non_root_time", "agg_hybridization per_non_root_time",
        "best_variant"
    )

    selected_instances = [
        "Hchl4s", "okp2", "Hchl7s", "okp3", "Hchl8s", "Hchl3s", "Hchl2", "CHL6",
        "CHL7", "Hchl6s", "CHL1", "CHL1s", "okp5"
    ]
    filter!(:instance_name => in(selected_instances), df)
    df[!, :order] = findfirst.(.==(df.instance_name), (selected_instances,))
    sort!(df, :order)
    select!(df, Not(:order))
    for col in names(df)
        (col == "instance_name" || col == "best_variant") && continue
        df[!, col] = number2latex.(round.((Int,), df[!, col]); enclose = false)
        #replace!(df[!, col], "0" => ">1")
    end
    for col in ("no_hybridization per_non_root_time", "hybridization per_non_root_time", "agg_hybridization per_non_root_time")
        replace!(df[!, col], "100" => ">99")
    end
    # Apply the best time bold just before returning.
    final_time_cols = ["no_hybridization mean", "percent hybridization mean", "percent agg_hybridization mean"]
    for (i, t) in enumerate(df[!, :best_variant])
        best_col = final_time_cols[t]
        df[i, best_col] = "\bestcolumnemph{$(df[i, best_col])}"
    end
    select!(df, Not(:best_variant))
    pretty_table(
        df; backend = :latex, nosubheader = true,
        alignment = vcat([:l], repeat([:r], ncol(df) - 1))
    )
end

# %%
let df = deepcopy(edf)
    sdf = filter(:instance_name => in(["CHL1s", "CHL1"]), df)
    cols = [:model_variant, :instance_name, :solver_seed, :build_and_solve_time]
    ssdf = sort(sdf, cols)
    sssdf = select(ssdf, cols)
end

# %%
# Create the comparison time table with extra runs
ectt = let df = deepcopy(edf)
    # build_and_solve_time: measures only what we want to measure, but does not exist if unfinished
    # run_total_time: measures things like instance reading too, but exists even if unfinished
    used_columns = [
        :instance_name, :model_variant, :build_and_solve_time, :run_total_time, :finished,
        :pevars, :cmvars, :plates, :build_stop_reason, :hybridizations,
        :first_borns_killed_by_hybridization
    ]
    select!(df, used_columns)
    @assert all(isone, df[!, :finished])
    cdf = combine(
        groupby(df, [:instance_name, :model_variant]),
        :build_and_solve_time => mean => :build_and_solve_time,
        :run_total_time => mean => :run_total_time,
        :finished => only_value => :finished,
        :pevars => only_value => :pevars,
        :cmvars => only_value => :cmvars,
        :plates => only_value => :plates,
        :build_stop_reason => only_value => :build_stop_reason,
        :hybridizations => only_value => :hybridizations,
        :first_borns_killed_by_hybridization => only_value => :first_borns_killed_by_hybridization
    )

    # We cannot remove the unfinished runs from the table, because things like
    # number of plates and variables of non-priced variants need to be summed even if
    # the run ended by timeout. This needs some jugglery with missing values.

    # Create a Float64 column indicating the best time for that instance in the table and
    # a Bool column indicating if the row has the best time for all table.
    cdf[!, :best_time_for_inst] = let cdf = deepcopy(cdf)
        replace!((v -> ismissing(v) ? Inf : v), cdf[!, :build_and_solve_time])
        b = @linq cdf |> select(:instance_name, :build_and_solve_time) |> groupby(:instance_name) |>
            based_on(; best_time = minimum(:build_and_solve_time))
        d = Dict{String, Float64}(map(=>, unique(cdf[!, :instance_name]), Iterators.cycle(Inf)))
        for row in eachrow(b)
            d[row.instance_name] = row.best_time
        end
        best_time_col = getindex.((d,), cdf[!, :instance_name])
    end

    cdf[!, :was_best] = cdf[!, :best_time_for_inst] .≈ cdf[!, :build_and_solve_time]
    cdf[!, :loss_on_non_best] = cdf[!, :build_and_solve_time] .- cdf[!, :best_time_for_inst]
    #sort!(select(cdf, :instance_name, :model_variant, :loss_on_non_best), [:instance_name, :model_variant])

    cdf[!, :build_and_solve_time] .= ifelse.(cdf[!, :finished], cdf[!, :build_and_solve_time], 0.0)
    #cdf[!, :variables] = cdf[!, :cmvars] .+ cdf[!, :pevars]

    #@assert iszero.(cdf[!, :variables]) == iszero.(cdf[!, :plates]) 
    cdf = @linq cdf |> groupby(:model_variant) |> based_on(;
        solved_time = sum(:build_and_solve_time),
        qt_best = sum(:was_best),
        time_loss = sum(:loss_on_non_best),
        pevars = sum(:pevars),
        cmvars = sum(:cmvars),
        plates = sum(:plates),
        qt_hybrid = sum(:hybridizations),
        qt_killed_fc = sum(:first_borns_killed_by_hybridization),
    )
    cdf[!, :per_hybrid] = (cdf[!, :qt_hybrid] ./ cdf[!, :cmvars]) .* 100
    cdf[!, :per_killed_fc] = (cdf[!, :qt_killed_fc] ./ cdf[!, :cmvars]) .* 100

    numeric_columns = [
        :solved_time, :time_loss, :qt_best, :pevars, :cmvars, :plates,
        :per_hybrid, :per_killed_fc
    ]
    for col in numeric_columns
        cdf[!, col] = number2latex.(round.((Int,), cdf[!, col]); enclose = false)
    end
    # Rename the columns to the names in the table
    select!(cdf,
        :model_variant => "Variant",
        :solved_time => "T. T.",
        :time_loss => "T. L.",
        :qt_best => "#b",
        :pevars => "#extr.",
        :cmvars => "#cuts",
        :per_hybrid => "h %",
        :per_killed_fc => "k %",
        :plates => "#plates",
    )
    select!(cdf, names(cdf) .=> esc_latex.(names(cdf)))
    # Rename the variants to the names in the table (also, define their order of appearance)
    pretty_variant_names = Dict{String, NamedTuple{(:pretty_name,:order), Tuple{String, Int}}}(
        "no_hybridization" => (pretty_name = "N. H.", order = 1),
        "hybridization" => (pretty_name = "C. H.", order = 2),
        "agg_hybridization" => (pretty_name = "A. H.", order = 3),
    )
    sort!(cdf, "Variant"; by = (name -> pretty_variant_names[name].order))
    cdf[!, "Variant"] = getindex.(getindex.((pretty_variant_names,), cdf[!, "Variant"]), :pretty_name)

    cdf
end

# #=
pretty_table(
    ectt; backend = :latex, nosubheader = true, alignment = vcat([:l], repeat([:r], ncol(ectt) - 1))
)
# =#

#showtable(ctt)

# %%
# Read the data, show nothing.
raw_oppdf_path = "./data/G2OPP_hybridization.csv"
raw_oppdf = DataFrame(CSV.File(raw_oppdf_path))
#nothing
showtable(raw_oppdf)

# %%
# Clean the data: replace sentinel values by missing; replace bool parameters by a single string id;
# create columns with the relevant pevars/cmvars/plates for a given configuration.
oppdf = let cdf = deepcopy(raw_oppdf)
    # Keep only the instance name (not the path).
    cdf[!, :instance_name] = basename.(cdf[!, :instance_path])
    cdf[!, :datafile] = basename.(cdf[!, :this_data_file])

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
    @assert all(==("none"), cdf[!, :pricing])
    for attr in (:cmvars, :pevars, :plates)
        cdf[!, attr] = cdf[!, "qt_$(attr)_after_preprocess"]
    end

    flags2variant!(cdf)

    cdf
end

# %%
# Read the data, show nothing.
raw_c42_stats_path = "./data/G2OPP_Clautiaux42_instance_stats.csv"
raw_c42_stats = DataFrame(CSV.File(raw_c42_stats_path))
#nothing
@show eltype.(identity.(eachcol(raw_c42_stats)))
showtable(raw_c42_stats)

# %%
# Clean the data: replace sentinel values by missing; replace bool parameters by a single string id;
# create columns with the relevant pevars/cmvars/plates for a given configuration.
raw_c42 = let df = deepcopy(raw_c42_stats)
    for col in [:rrl, :rrw]
        df[!, col] = parse.((Rational{BigInt},), df[!, col])
    end
    df
end

# %%
# Create the comparison time table with extra runs
opp_ctt = let df = deepcopy(oppdf)
    # build_and_solve_time: measures only what we want to measure, but does not exist if unfinished
    # run_total_time: measures things like instance reading too, but exists even if unfinished
    used_columns = [
        :instance_name, :model_variant, :build_and_solve_time, :run_total_time, :finished,
        :pevars, :cmvars, :plates, :build_stop_reason, :hybridizations,
        :first_borns_killed_by_hybridization
    ]
    select!(df, used_columns)
    @assert all(isone, df[!, :finished])
    cdf = combine(
        groupby(df, [:instance_name, :model_variant]),
        :build_and_solve_time => mean => :build_and_solve_time,
        :run_total_time => mean => :run_total_time,
        :finished => only_value => :finished,
        :pevars => only_value => :pevars,
        :cmvars => only_value => :cmvars,
        :plates => only_value => :plates,
        :build_stop_reason => only_value => :build_stop_reason,
        :hybridizations => only_value => :hybridizations,
        :first_borns_killed_by_hybridization => only_value => :first_borns_killed_by_hybridization
    )

    # We cannot remove the unfinished runs from the table, because things like
    # number of plates and variables of non-priced variants need to be summed even if
    # the run ended by timeout. This needs some jugglery with missing values.

    # Create a Float64 column indicating the best time for that instance in the table and
    # a Bool column indicating if the row has the best time for all table.
    cdf[!, :best_time_for_inst] = let cdf = deepcopy(cdf)
        replace!((v -> ismissing(v) ? Inf : v), cdf[!, :build_and_solve_time])
        b = @linq cdf |> select(:instance_name, :build_and_solve_time) |> groupby(:instance_name) |>
            based_on(; best_time = minimum(:build_and_solve_time))
        d = Dict{String, Float64}(map(=>, unique(cdf[!, :instance_name]), Iterators.cycle(Inf)))
        for row in eachrow(b)
            d[row.instance_name] = row.best_time
        end
        best_time_col = getindex.((d,), cdf[!, :instance_name])
    end

    cdf[!, :was_best] = cdf[!, :best_time_for_inst] .≈ cdf[!, :build_and_solve_time]
    cdf[!, :loss_on_non_best] = cdf[!, :build_and_solve_time] .- cdf[!, :best_time_for_inst]
    #sort!(select(cdf, :instance_name, :model_variant, :loss_on_non_best), [:instance_name, :model_variant])

    cdf[!, :build_and_solve_time] .= ifelse.(cdf[!, :finished], cdf[!, :build_and_solve_time], 0.0)
    #cdf[!, :variables] = cdf[!, :cmvars] .+ cdf[!, :pevars]

    #@assert iszero.(cdf[!, :variables]) == iszero.(cdf[!, :plates]) 
    cdf = @linq cdf |> groupby(:model_variant) |> based_on(;
        solved_time = sum(:build_and_solve_time),
        qt_best = sum(:was_best),
        time_loss = sum(:loss_on_non_best),
        pevars = sum(:pevars),
        cmvars = sum(:cmvars),
        plates = sum(:plates),
        qt_hybrid = sum(:hybridizations),
        qt_killed_fc = sum(:first_borns_killed_by_hybridization),
    )
    cdf[!, :per_hybrid] = (cdf[!, :qt_hybrid] ./ cdf[!, :cmvars]) .* 100
    cdf[!, :per_killed_fc] = (cdf[!, :qt_killed_fc] ./ cdf[!, :cmvars]) .* 100

    numeric_columns = [
        :solved_time, :time_loss, :qt_best, :pevars, :cmvars, :plates,
        :per_hybrid, :per_killed_fc
    ]
    for col in numeric_columns
        cdf[!, col] = number2latex.(cdf[!, col]; enclose = false)
    end
    # Rename the columns to the names in the table
    select!(cdf,
        :model_variant => "Variant",
        :solved_time => "T. T.",
        :time_loss => "T. L.",
        :qt_best => "#b",
        :pevars => "#extr.",
        :cmvars => "#cuts",
        :per_hybrid => "h %",
        :per_killed_fc => "k %",
        :plates => "#plates",
    )
    select!(cdf, names(cdf) .=> esc_latex.(names(cdf)))
    # Rename the variants to the names in the table (also, define their order of appearance)
    pretty_variant_names = Dict{String, NamedTuple{(:pretty_name,:order), Tuple{String, Int}}}(
        "no_hybridization" => (pretty_name = "N. H.", order = 1),
        "hybridization" => (pretty_name = "C. H.", order = 2),
        "agg_hybridization" => (pretty_name = "A. H.", order = 3),
    )
    sort!(cdf, "Variant"; by = (name -> pretty_variant_names[name].order))
    cdf[!, "Variant"] = getindex.(getindex.((pretty_variant_names,), cdf[!, "Variant"]), :pretty_name)

    cdf
end

# #=
pretty_table(
    opp_ctt; backend = :latex, nosubheader = true, alignment = vcat([:l], repeat([:r], ncol(ectt) - 1))
)
# =#

#showtable(ctt)

# %%
let df = deepcopy(raw_c42)
    @show mean(df[!, :rrw])
    @show mean(df[!, :rrl])
    @show Float64(mean(df[!, :rrw]))
    @show Float64(mean(df[!, :rrl]))
end

# %%
opp_ussdf = let df = deepcopy(oppdf)
    #=
    replace!(
        df[!, :model_variant],
        "no_hybridization" => "Base", "hybridization" => "Hyb.",
        "agg_hybridization" => "Agg. H."
    )
    =#
    
    sdf = combine(
        groupby(df, [:instance_name, :model_variant]),
        :build_and_solve_time => length => :nrow,
        :build_and_solve_time => mean => :mean,
        :build_and_solve_time => std  => :std,
        :build_and_solve_time => minimum => :min,
        :build_and_solve_time => maximum => :max,
        [:build_and_solve_time, :root_node_time] => ((b, r) -> (sum(r)/sum(b)) * 100) => :per_avg_barrier_time,
        [:build_and_solve_time, :root_node_time] => ((b, r) -> (1.0 - (sum(r)/sum(b))) * 100) => :per_non_root_time,
        [:hybridizations, :pevars, :cmvars] => ((h, p, e) -> 100.0 * (sum(h) / (sum(p) + sum(e)))) => :h_per,
        [:first_borns_killed_by_hybridization, :pevars, :cmvars] => ((k, p, e) -> 100.0 * (sum(k) / (sum(p) + sum(e)))) => :k_per,
    )
    @assert all(==(10), sdf[!, :nrow])
    select!(sdf, Not(:nrow))
    sdf[!, :variant_order] = replace(sdf[!, :model_variant],
        "no_hybridization" => 1, "hybridization" => 2, "agg_hybridization" => 3
    )
    stacked_cols = [:h_per, :k_per, :mean, :std, :min, :max, :per_avg_barrier_time, :per_non_root_time]
    ssdf = stack(
        sdf, stacked_cols;
        variable_eltype = Symbol
    )
    ssdf[!, :stat_order] = replace(ssdf[!, :variable],
        (stacked_cols .=> 1:length(stacked_cols))...
    )
    ssdf[!, :order] = @. (ssdf[!, :stat_order] * 10) + ssdf[!, :variant_order]
    @. ssdf[!, :variable] = string(ssdf[!, :order]) * " " * ssdf[!, :model_variant] *
        " " * string(ssdf[!, :variable])
    select!(ssdf, Not([:model_variant, :variant_order, :stat_order, :order]))
    
    ussdf = unstack(ssdf)
    @show names(ussdf)
    # COMPUTE PERCENT COLUMNS
    non_key_names = filter(!=("instance_name"), names(ussdf))
    rename!(ussdf, non_key_names .=> join.(getindex.(split.(non_key_names, " "), ([2, 3],)), " "))
    # #=
    for variant in ("hybridization", "agg_hybridization")
        for stat in ("mean", #=, "std"=#)
            ussdf[!, "percent $variant $stat"] = (
                ((ussdf[!, "$variant $stat"] ./ ussdf[!, "no_hybridization $stat"]) #=.- 1.0=#) .* 100
            )
        end
    end
    # Changes standard deviation to coefficient of deviation
    for variant in ("no_hybridization", "hybridization", "agg_hybridization")
        ussdf[!, "$variant cv"] = (
            ((ussdf[!, "$variant std"] ./ ussdf[!, "$variant mean"]) #=.- 1.0=#) .* 100
        )
    end
    # =#
    
    #showtable(ussdf)
    ussdf
end

# %%
let df = deepcopy(opp_ussdf)
    # Remove min and max columns
    #select!(df, filter(s -> !occursin("min", s) && !occursin("max", s), names(df)))

    times = collect(zip(
        eachcol(select(df, ["no_hybridization mean", "hybridization mean", "agg_hybridization mean"]))...
    ))
    @show names(eachcol(select(df, ["no_hybridization mean", "hybridization mean", "agg_hybridization mean"])))
    df[!, :best_variant] = last.(findmin.(times))

    select!(df,
        :instance_name, "hybridization h_per", "agg_hybridization h_per",
        "no_hybridization mean", "percent hybridization mean", "percent agg_hybridization mean",
        #="no_hybridization std", "percent hybridization std", "percent agg_hybridization std",=#
        "no_hybridization cv", "hybridization cv", "agg_hybridization cv",
        "no_hybridization per_non_root_time", "hybridization per_non_root_time", "agg_hybridization per_non_root_time",
        "best_variant"
    )

    sort!(df, "no_hybridization mean"; rev = true)
    for col in names(df)
        (col == "instance_name" || col == "best_variant") && continue
        if col == "no_hybridization mean"
            df[!, col] = number2latex.(round.(df[!, col]; digits = 2); enclose = false)
        else
            df[!, col] = number2latex.(round.((Int,), df[!, col]); enclose = false)
        end
        #replace!(df[!, col], "0" => ">1")
    end
    for col in ("no_hybridization per_non_root_time", "hybridization per_non_root_time", "agg_hybridization per_non_root_time")
        replace!(df[!, col], "100" => ">99")
    end
    # Apply the best time bold just before returning.
    final_time_cols = ["no_hybridization mean", "percent hybridization mean", "percent agg_hybridization mean"]
    for (i, t) in enumerate(df[!, :best_variant])
        best_col = final_time_cols[t]
        df[i, best_col] = "\bestcolumnemph{$(df[i, best_col])}"
    end
    select!(df, Not(:best_variant))
    pretty_table(
        df; backend = :latex, nosubheader = true,
        alignment = vcat([:l], repeat([:r], ncol(df) - 1))
    )
end

# %%
let df = deepcopy(raw_oppdf)
    df[!, :instance_name] = basename.(df[!, :instance_path])
    flags2variant!(df)
    replace!(
        df[!, :model_variant],
        "no_hybridization" => "Base", "hybridization" => "Hyb.",
        "agg_hybridization" => "Agg. H."
    )
    df[!, :order] = replace(df[!, :model_variant],
        "Base" => 1, "Hyb." => 2, "Agg. H." => 3
    )
    #df[!, :discrete_seed] = string.(df[!, :solver_seed])
    sort!(df, [:order, :solver_seed])
    
    # #=
    ps = map(instance -> plot(
        filter(:instance_name => ==(instance), df), Guide.title(instance),
        Guide.ylabel("Time (s)"; orientation = :vertical),
        Guide.xlabel("Variant"; orientation = :horizontal),
        Scale.y_continuous(; format = :plain),
        #color = "discrete_seed",
        x = "model_variant", y = "build_and_solve_time", Geom.point
    ), unique(df[!, :instance_name]))
    vs = vstack(ps...)
    #gp = gridstack(reshape(ps, (2, 2)))
    vs |> PS("../plots/hyb_g2opp_seeded.ps")
    # =#
    #=
    plot(
        df, xgroup = "instance_name", x = "build_and_solve_time", y = "model_variant",
        Geom.subplot_grid(Geom.point)
    )
    =#
end
