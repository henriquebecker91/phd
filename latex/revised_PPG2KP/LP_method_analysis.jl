# %%
# Ensure your working directory is: https://github.com/henriquebecker91/phd/tree/master/latex/revised_PPG2KP
include("notebook_setup.jl")

# %%
# Read the data, show nothing.
lp_method_csv_path = "./data/lp_method.csv"
lp_method_df = DataFrame(CSV.File(lp_method_csv_path))
#showtable(lp_method_df)
nothing

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
# Shows the cleaned data.
println(names(lp_method_df))
showtable(lp_method_df)

# %%
# Just checking which instances ended in timeout.
@linq lp_method_df |>
    where(:finished .== false) |>
    select(:instance_name, :lp_method, :lp_method_switch, :pricing_method, :total_instance_time)

# %%
sel_columns = [:instance_name, :lp_method, :lp_method_switch, :pricing_method, :build_and_solve_time]
# long (opposed to wide) table of total time
lttt = select(lp_method_df, sel_columns)
lttt.params = join.(zip(lttt.lp_method, lttt.lp_method_switch, lttt.pricing_method), "_")
pretty_params = Dict{String, String}(
    "barrier_only_none" => "NP B",
    "barrier_switch_furini" => "FP DS+B",
    "simplex_only_furini" => "FP DS",
    "barrier_only_furini" => "FP B",
    "simplex_only_none" => "NP DS"
)
lttt.params = getindex.((pretty_params,), lttt.params)
lttt = select(lttt, :instance_name, :params, :build_and_solve_time)
lttt.build_and_solve_time = (x -> @sprintf("%.2f", x)).(lttt.build_and_solve_time)
lttt.build_and_solve_time = ifelse.(
    lttt.build_and_solve_time .== "NaN", "--", lttt.build_and_solve_time
)
# Table with instance, params, and time, in long format.
showtable(lttt)

# %%
# Table with the build+solve time for each instance (rows) and parameter combination (columns).
# Note this counts the theoretically unnecessary root node solving of the priced models.
# wide (opposed to long) table of total time
wttt = unstack(lttt, :instance_name, :params, :build_and_solve_time)

showtable(wttt)


# %%
#pretty_table(wide_ttt; backend = :latex)

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
