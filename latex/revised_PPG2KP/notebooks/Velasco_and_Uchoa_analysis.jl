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
#     display_name: Julia 1.4.1
#     language: julia
#     name: julia-1.4
# ---

# %%
# Ensure your working directory is: https://github.com/henriquebecker91/phd/tree/master/latex/revised_PPG2KP
include("notebook_setup.jl")

# %%
# Read the data, show nothing.
raw_vudf_path = "./data/vu_experiment.csv"
raw_vudf = DataFrame(CSV.File(raw_vudf_path))
@show collect(zip(names(raw_vudf), typeof.(eachcol(raw_vudf))))
#showtable(raw_vudf)
nothing

# %%
# Just assert some basic things, and show the solved runs.
let raw_vudf = deepcopy(raw_vudf)
    @show size(raw_vudf)
    #describe(raw_vudf)
    @assert @with(raw_vudf, all(:had_timeout .== :run_ended_by_exception))
    #@assert @with(raw_vudf, .!((isnan.(:solution_profit) .& isnan.(:obj_value)) .|
    #    (:solution_profit .≈ :obj_value)))
    @show 160 - sum(raw_vudf[!, :had_timeout])
    # Should be an empty table. If it is not, then the solution assembled from model variables has
    # a different profit value from the objective function.  
    @assert isempty(@where(raw_vudf, .!((isnan.(:solution_profit) .& isnan.(:obj_value)) .|
        (:solution_profit .≈ :obj_value))))
    solved = @linq raw_vudf |> select(:instance_name, :pricing_method, :run_total_time, :solution_profit, :had_timeout) |>
        where(.!:had_timeout)
end

# %%
# Clean the data: replace sentinel values by missing; replace full paths by only the filename.
# No need to create config id (because there is only one parameter that changes between them);
# no need to find the right number of vars/plates (because for both is just the after_purge).
vudf = let vudf = deepcopy(raw_vudf)
    # Keep only the instance name (not the path).
    @with(vudf, :instance_name .= basename.(:instance_name))
    @with(vudf, :datafile .= basename.(:datafile))
    colnames = names(vudf)
    for name in colnames
        column = vudf[!, name]
        if eltype(column) <: AbstractFloat
            vudf[!, name] = ifelse.(isnan.(column), missing, column)
        elseif eltype(column) <: Integer
            vudf[!, name] = ifelse.(column .< 0, missing, column)
        end
    end
    
    vudf
end
showtable(vudf)
#nothing

# %%
@show names(vudf)
@show unique(vudf[!, :build_stop_reason])
@show @where(vudf, :had_timeout .!= :run_ended_by_exception)


# %%
vu2019_path = "./data/vu_2019.csv"
vu2019 = DataFrame(CSV.File(vu2019_path))
#nothing
showtable(vu2019)

# %%
# Get some statistics regarding the heuristic.
import Statistics: mean
stats_heuristic = let vudf = deepcopy(vudf)
    # Normalizes how some information appears in the restricted and unrestricted steps.
    vudf[!, :restricted_had_timeout] = occursin.("TIME_LIMIT", vudf[!, :restricted_stop_reason]) .|
        occursin.("TIME_LIMIT", vudf[!, :restricted_lp_stop_reason])
    vudf[!, :built_restricted_model] = .!occursin.("TIME_LIMIT", vudf[!, :restricted_lp_stop_reason])
    vudf[!, :built_model] = vudf[!, :build_stop_reason] .== "BUILT_MODEL"
    # Sets of columns used extensively below.
    shared_columns = [:instance_name, :pricing_method]
    restricted_columns = [
        :heuristic_lb, :restricted_obj_value, :restricted_had_timeout, :built_restricted_model
    ]
    unrestricted_columns = [
        :heuristic_lb, :obj_value, :had_timeout, :built_model
    ]
    # Clean the dataframe.
    used_columns = unique(vcat(shared_columns, restricted_columns, unrestricted_columns))
    select!(vudf, used_columns) # Get only the used_columns
    # Create a restricted dataframe, with the same column names as the unrestricted counterpart.
    restricted = select(vudf, vcat(shared_columns, restricted_columns))
    filter!("pricing_method" => isequal("furini"), restricted)
    select!(restricted, vcat(shared_columns, restricted_columns .=> unrestricted_columns)...)
    restricted[!, :pricing_method] .= "restricted" 
    # Create a single dataframe with restricted data as extra rows, not extra columns in
    # priced rows.
    select!(vudf, vcat(shared_columns, unrestricted_columns))
    vudf = vcat(vudf, restricted)
    filter!(:built_model => identity, vudf)
    #@show vudf
    select!(vudf, Not(:pricing_method), :pricing_method => :variant)
    vudf[!, :heuristic_gap] = @with vudf (100.0 .- (:heuristic_lb ./ :obj_value) .* 100.0)
    @show(@where(vudf, iszero.(:heuristic_gap)))
    #filter!([:instance_name, :variant] => !((i, v) -> i == "P3_250_250_25_3" && v == "none"), vudf)
    @show maximum(vudf[!, :heuristic_gap])
    @show minimum(vudf[!, :heuristic_gap])
    @show mean(vudf[!, :heuristic_gap])
    @linq select(vudf, :heuristic_gap, :variant) |> groupby(:variant) |>
        based_on(; mean = mean(:heuristic_gap), max = maximum(:heuristic_gap), min = minimum(:heuristic_gap))
    #sort!(vudf, [:instance_name, :variant])
end
#nothing
showtable(stats_heuristic)

# %%
# Summary dataframe. Dataframe with summarized time, #variables, etc...
function bitmask_of_rows_with_missing_values(df)
    return reduce(broadwrap(|), (map(broadwrap(ismissing), eachcol(df))))
end
sdf = let vudf = deepcopy(vudf)
    # Normalizes how some information appears in the restricted and unrestricted steps.
    vudf[!, :restricted_had_timeout] = occursin.("TIME_LIMIT", vudf[!, :restricted_stop_reason]) .|
        occursin.("TIME_LIMIT", vudf[!, :restricted_lp_stop_reason])
    vudf[!, :restricted_early_opt] = .!vudf[!, :solved_priced_restricted_model] .&
        .!occursin.("TIME_LIMIT", vudf[!, :restricted_stop_reason]) .&
        .!occursin.("TIME_LIMIT", vudf[!, :restricted_lp_stop_reason])
    vudf[!, :early_opt] = vudf[!, :build_stop_reason] .== "FOUND_OPTIMUM"
    
    # Sets of columns used extensively below.
    shared_columns = [:instance_name, :pricing_method]
    restricted_columns = [
        :qt_pevars_priced_restricted,
        :qt_cmvars_priced_restricted,
        :qt_plates_priced_restricted,
        :restricted_pricing_time,
        :restricted_had_timeout,
        :restricted_early_opt,
    ]
    unrestricted_columns = [
        :qt_pevars_after_purge,
        :qt_cmvars_after_purge,
        :qt_plates_after_purge,
        :build_and_solve_time,
        :had_timeout,
        :early_opt,
    ]
    # Clean the dataframe.
    used_columns = vcat(shared_columns, restricted_columns, unrestricted_columns)
    select!(vudf, used_columns) # Get only the used_columns
    # Create a restricted dataframe, with the same column names as the unrestricted counterpart.
    restricted = select(vudf, vcat(shared_columns, restricted_columns))
    filter!("pricing_method" => isequal("furini"), restricted)
    select!(restricted, vcat(shared_columns, restricted_columns .=> unrestricted_columns)...)
    restricted[!, :pricing_method] .= "restricted" 
    # Create a single dataframe with restricted data as extra rows, not extra columns in
    # priced rows.
    select!(vudf, vcat(shared_columns, unrestricted_columns))
    vudf = vcat(vudf, restricted)
    # Create the instance groups and remove/rename some columns.
    vudf[!, :group] = (m -> parse(Int, m.captures[1])).(
        match.(r"^P([1-4]).*", vudf[!, :instance_name])
    )
    select!(vudf, Not(:instance_name))
    select!(vudf, Not(:pricing_method), :pricing_method => :variant)
    # If no instance makes the pricing return early, then we can exclude this from the comparison. 
    @assert iszero(sum(vudf[!, :early_opt]))
    select!(vudf, Not(:early_opt))
    @assert all(ifelse.(
        bitmask_of_rows_with_missing_values(vudf), # If a row has some missing value
        vudf[!, :had_timeout], # then the same row has an indication of timeout
        true # if there is not missing value, then the run can have hit the time limit or not
    ))
    #@show maximum(skipmissing(vudf[!, :build_and_solve_time]))
    #showtable(@where(vudf, (:group .== 4) .& (:variant .== "furini")))
    vudf = @linq vudf |> groupby([:group, :variant]) |>
        based_on(;
            solved = sum(.!:had_timeout),
            variables = sum(skipmissing(:qt_pevars_after_purge)) + sum(skipmissing(:qt_cmvars_after_purge)),
            plates = sum(skipmissing(:qt_plates_after_purge)),
            built = sum(.!ismissing.(:qt_cmvars_after_purge)),
            total_time = sum(replace(
                :build_and_solve_time, missing => (3 * 60 * 60)
            ))
        ) |> sort([:group, :variant])
    vudf
end
#nothing
#showtable(@where(sdf, .!reduce(broadwrap(|), (map(broadwrap(ismissing), eachcol(sdf))))))

# %%
st = let sdf = deepcopy(sdf)
    variant_order = Dict{String, Int}(
        "none"       => 1,
        "restricted" => 2,
        "furini"     => 3,
    )
    sdf[!, :variant_order] = (e -> variant_order[e]).(sdf[!, :variant])
    sort!(sdf, [:group, :variant_order])
    select!(sdf, Not(:variant_order))
    pretty_names = Dict{String, String}(
        "none"       => "Not Priced",
        "restricted" => "Restricted Priced",
        "furini"     => "Priced",
    )
    replace!((e -> pretty_names[e]), sdf[!, :variant])
    # Transform some values in averages
    sdf[!, :variables] ./= sdf[!, :built]
    sdf[!, :plates] ./= sdf[!, :built]
    sdf[!, :solved_time] = sdf[!, :total_time] .- ((20 .- sdf[!, :solved]) .* (3 * 60 * 60))
    sdf[!, :solved_time] ./= sdf[!, :solved]
    sdf[!, :total_time] = trunc.(Int, sdf[!, :total_time])
    # define header names and column order
    select!(sdf,
        :group => "G.",
        :variant => "Variant",
        :built => "#m",
        :variables => #=(c -> c ./ sdf[!, :built]) =>=# "Avg. #v",
        :plates => #=(c -> c ./ sdf[!, :built]) =>=# "Avg. #p",
        :total_time => "T. T.",
        :solved => "#s",
        :solved_time => #=(c -> c ./ sdf[!, :solved]) =>=# "Avg. S. T.",
    )
    select!(sdf, names(sdf) .=> esc_latex.(names(sdf)))
    for column in names(sdf) # latexify all numeric columns
        if nonmissingtype(eltype(sdf[!, column])) <: Number
            sdf[!, column] = number2latex.(sdf[!, column]; enclose = false)
        end
    end
    sdf
end

pretty_table(
    st; backend = :latex, nosubheader = true, alignment = vcat([:l], repeat([:r], ncol(st) - 1))
)

# %%
# New Results DataFrame. Dataframe detailing new results obtained.
nrdf = let vudf = deepcopy(vudf), vu2019 = select(vu2019, :instance_name, :lb => :vu_lb, :ub => :vu_ub)
    used_columns = [
        :instance_name, :pricing_method, :heuristic_lb, :restricted_obj_value, :restricted_obj_bound,
        :obj_value, :obj_bound
    ]
    select!(vudf, used_columns) # Get only the used_columns
    for column in names(vudf) # Integralize all numerical columns
        column in ("instance_name", "pricing_method") && continue
        #@assert all(v -> abs(round(v, RoundNearest) - v) < 1e-6, filter(!ismissing, vudf[!, column]))
        vudf[!, column] = trunc.(Union{Int,Missing}, vudf[!, column] .+ 1e-6)
    end
    # Separate the priced and non-priced rows to after join them.
    priced = filter!(:pricing_method => isequal("furini"), deepcopy(vudf))
    select!(priced,
        :instance_name, :heuristic_lb => :ph_lb, :restricted_obj_value => :r_lb,
        :restricted_obj_bound => :r_ub, :obj_value => :p_lb, :obj_bound => :p_ub
    )
    non_priced = filter!(:pricing_method => isequal("none"), deepcopy(vudf))
    select!(non_priced,
        :instance_name, :heuristic_lb => :nph_lb, :obj_value => :np_lb, :obj_bound => :np_ub
    )
    vudf = innerjoin(priced, non_priced, vu2019; on = :instance_name)
    @assert @with(vudf, all(:ph_lb .== :nph_lb))
    select!(vudf, Not(:nph_lb), :ph_lb => :h_lb)
    @show names(vudf)
    vudf = @linq vudf |> where(
            (.!ismissing.(:p_lb) .& ((:p_lb .== :p_ub) .| (:p_lb .> :vu_lb) .| (:p_ub .< :vu_ub))) .|
            (.!ismissing.(:r_lb) .& ((:r_lb .== :r_ub) .| (:r_lb .> :vu_lb))) .|
            (.!ismissing.(:np_lb) .& ((:np_lb .== :np_ub) .| (:np_lb .> :vu_lb) .| (:np_ub .< :vu_ub)))
        )
    
    p_h_gap = @with vudf skipmissing(100.0 .- (:h_lb ./ :p_lb) .* 100.0)
    @show minimum(p_h_gap)
    @show maximum(p_h_gap)
    @show mean(p_h_gap)
    
    r_h_gap = @with vudf skipmissing(100.0 .- (:h_lb ./ :r_lb) .* 100.0)
    @show minimum(r_h_gap)
    @show maximum(r_h_gap)
    @show mean(r_h_gap)
    
    np_h_gap = @with vudf skipmissing(100.0 .- (:h_lb ./ :np_lb) .* 100.0)
    @show minimum(np_h_gap)
    @show maximum(np_h_gap)
    @show mean(np_h_gap)

    vudf
end

showtable(nrdf)

# %%
nrt = let nrdf = deepcopy(nrdf)
    select!(nrdf, Not(:h_lb)) # For now, let us disconsider the heuristic column.
    nrdf[!, :instance_name] = esc_latex.(nrdf[!, :instance_name])
    # Sometimes the restricted priced was stopped too early and gave terrible upper bounds,
    # this breaks the column layout and bring no useful information, better to display a dash instead.
    nrdf[!, :r_ub] .= replace(x -> (ismissing(x) ? missing : x > 100000 ? missing : x), nrdf[!, :r_ub])
    for column in names(nrdf) # latexify all numeric columns
        column == "instance_name" && continue
        nrdf[!, column] = number2latex.(nrdf[!, column]; enclose = false)
    end
    # First let us hightlight (with textbf), the best values among the LB, and the best values
    # among the UB.
    highlight_best_values!(nrdf;
        columns = [#=:h_lb,=# :r_lb, :p_lb, :np_lb, :vu_lb], chooser = find_index_allmax,
        cleaner = dirt2number
    )
    highlight_best_values!(nrdf;
        columns = [:p_ub, :np_ub, :vu_ub], chooser = find_index_allmin, cleaner = dirt2number
    )
    # Then, let us italicize the UBs and LBs from the same method, if they have the same value.
    select_all_if_equal(a) = all(e -> e == a[1], a) ? collect(Int, keys(a)) : Int[]
    latex_wrap(command, s) = "\\$(command){$(s)}"
    struct LatexWrapper; command :: String; end
    function (w::LatexWrapper)(s); return "\\$(w.command){$(s)}"; end
    latex_wrap(command) = LatexWrapper(command)
    highlight_best_values!(nrdf;
        columns = [:r_lb, :r_ub], chooser = select_all_if_equal,
        cleaner = dirt2number, changer = latex_wrap("underline")
    )
    highlight_best_values!(nrdf;
        columns = [:p_lb, :p_ub], chooser = select_all_if_equal,
        cleaner = dirt2number, changer = latex_wrap("underline")
    )
    highlight_best_values!(nrdf;
        columns = [:np_lb, :np_ub], chooser = select_all_if_equal,
        cleaner = dirt2number, changer = latex_wrap("underline")
    )
    highlight_best_values!(nrdf;
        columns = [:vu_lb, :vu_ub], chooser = select_all_if_equal,
        cleaner = dirt2number, changer = latex_wrap("underline")
    )
    sort!(nrdf, :instance_name)
    select!(nrdf,
        :instance_name => "Instance",
        #:h_lb => "H LB",
        :r_lb => "R LB",
        :p_lb => "P LB",
        :np_lb => "NP LB",
        :vu_lb => "VU LB",
        :r_ub => "R UB",
        :p_ub => "P UB",
        :np_ub => "NP UB",
        :vu_ub => "V UB",
    )
    select!(nrdf, names(nrdf) .=> esc_latex.(names(nrdf)))
    nrdf
end

pretty_table(
    nrt; backend = :latex, nosubheader = true, alignment = vcat([:l], repeat([:r], ncol(nrt) - 1))
)

# %%
