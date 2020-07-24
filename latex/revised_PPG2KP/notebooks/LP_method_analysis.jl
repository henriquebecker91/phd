# ---
# jupyter:
#   jupytext:
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
lp_method_csv_path = "./data/lp_method.csv"
lp_method_df = DataFrame(CSV.File(lp_method_csv_path))
showtable(lp_method_df)
#nothing

# %%
# Clean the data a little. Try to keep this safe to re-apply to a cleaned dataframe, if possible.
lp_method_df = let
    # Keep only the instance name (not the path).
    @with(lp_method_df, :instance_name .= basename.(:instance_name))
    @with(lp_method_df, :datafile .= basename.(:datafile))
    # Replace parameter numeric values with more descriptive strings.
    lp_method_code2name = Dict{Int, String}(1 => "simplex", 2 => "barrier")
    lp_switch_code2name = Dict{Int, String}(-2 => "only", 1 => "switch")
    if eltype(lp_method_df.lp_method) == Int
        lp_method_df.lp_method = getindex.((lp_method_code2name,), lp_method_df.lp_method)
    end
    if eltype(lp_method_df.lp_method_switch) == Int
        lp_method_df.lp_method_switch = getindex.((lp_switch_code2name,), lp_method_df.lp_method_switch)
    end
    # Change the not-a-number values in final_root_time to zero. They happen  
    lp_method_df
end
nothing

# %%
let lp_method_df = deepcopy(lp_method_df)
    finished = @linq lp_method_df |> select(:instance_name, :solution_profit, :finished)
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
# Shows the cleaned data.
println(names(lp_method_df))
showtable(lp_method_df)

# %%
# Just checking which instances ended by time limit or memory limit.
@linq lp_method_df |>
    where(.!:finished) |>
    select(:instance_name, :lp_method, :lp_method_switch, :pricing_method, :total_instance_time)

# %%
lp_params_df = let lp_method_df = deepcopy(lp_method_df)
    @show length(join.(zip(
        lp_method_df.lp_method, lp_method_df.lp_method_switch, lp_method_df.pricing_method
    ), "_"))
    @show length(lp_method_df.lp_method)
    lp_method_df[!, "params"] = join.(zip(
        lp_method_df.lp_method, lp_method_df.lp_method_switch, lp_method_df.pricing_method
    ), "_")
    pretty_params = Dict{String, String}(
        "barrier_only_none" => "NP B",
        "barrier_switch_furini" => "FP DS+B",
        "simplex_only_furini" => "FP DS",
        "barrier_only_furini" => "FP B",
        "simplex_only_none" => "NP DS"
    )
    lp_method_df[!, "params"] = getindex.((pretty_params,), lp_method_df.params)
    lp_method_df
end

# %%
lttt = let lp_params_df = lp_params_df
    sel_columns = [
        :instance_name, :params, :build_and_solve_time, :final_root_time
    ]
    # long (opposed to wide) table of total time
    lttt = select(lp_params_df, sel_columns)
    percents = @linq lttt |>
        where((:params .== "NP B") .| (:params .== "NP DS")) |>
        transform(value = (:final_root_time ./ :build_and_solve_time) .* 100.0) |>
        transform(variable = ifelse.(:params .== "NP B", "% B RNR", "% DS RNR")) |>
        select!(:instance_name, :variable, :value)
    
    select!(lttt, :instance_name, :params => :variable, :build_and_solve_time => :value)
    lttt = vcat(lttt, percents; cols = :setequal)

    lttt.value = (x -> @sprintf("%.2f", x)).(lttt.value)
    lttt.value = ifelse.(
        lttt.value .== "NaN", "--", string.(lttt.value)
    )

    # Table with instance, params, and time, in long format.
    lttt
end
nothing
showtable(lttt)

# %%
# Table with the build+solve time for each instance (rows) and parameter combination (columns).
# Note this counts the theoretically unnecessary root node solving of the priced models.
# wide (opposed to long) table of total time
wttt = unstack(lttt)

let
    latex_table = deepcopy(wttt)
    rename!(latex_table, :instance_name => :Name)
    desired_col_order = ["NP DS", "\\% DS RNR", "FP DS", "NP B", "\\% B RNR", "FP B", "FP DS+B"]
    rename!(latex_table, "% DS RNR" => "\\% DS RNR", "% B RNR" => "\\% B RNR")
    select!(latex_table, :Name, desired_col_order)
    highlight_best_values!(latex_table; columns = ["FP B", "FP DS", "FP DS+B", "NP B", "NP DS"])
    pretty_table(
        latex_table; backend = :latex, nosubheader = true, alignment = [:l; repeat([:r], 7)]
    )
end


# %%
# Now let us do the same but without the root node time.
lcttt = let
    sel_columns = [
        :instance_name, :lp_method, :lp_method_switch, :pricing_method, :build_and_solve_time,
        :final_root_time
    ]
    # long (opposed to wide) corrected table of total time
    lcttt = select(lp_method_df, sel_columns)
    lcttt.corrected_time = lcttt.build_and_solve_time .- lcttt.final_root_time
    lcttt.both_timings = tuple.(lcttt.build_and_solve_time, lcttt.corrected_time)
    lcttt.params = join.(zip(lcttt.lp_method, lcttt.lp_method_switch, lcttt.pricing_method), "_")
    pretty_params = Dict{String, String}(
        "barrier_only_none" => "NP B",
        "barrier_switch_furini" => "FP DS+B",
        "simplex_only_furini" => "FP DS",
        "barrier_only_furini" => "FP B",
        "simplex_only_none" => "NP DS"
    )
    lcttt.params = getindex.((pretty_params,), lcttt.params)
    lcttt = select(lcttt, :instance_name, :params, :both_timings)
    lcttt
end
showtable(lcttt)

# %%
# Table with both the build+solve time and this time minus the root node solving for each instance (rows)
# and parameter combination (columns).
# wide (opposed to long) corrected table of total time
wcttt = unstack(lcttt, :instance_name, :params, :both_timings)
clean_col(df, col) = setproperty!(df, Symbol(col), (both -> round.(both; digits = 2)).(getproperty(df, Symbol(col))))
clean_col.((wcttt,), ["FP B", "FP DS", "FP DS+B", "NP B", "NP DS"])

pretty_table(wcttt; backend = :latex)
#showtable(wcttt)


# %%
?setproperty!
