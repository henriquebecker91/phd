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
# Ensure your working directory is: https://github.com/henriquebecker91/phd/tree/master/latex/thesis/notebooks
include("notebook_setup.jl")

# %%
ais_path = "./data/G2CSP_A_instances_stats.csv"
ais = DataFrame(CSV.File(ais_path))
rename!(ais, :filename => :instance_name)
#nothing
showtable(ais)

# %%
raw_adf_path = "./data/G2CSP_A.csv"
raw_adf = DataFrame(CSV.File(raw_adf_path))
#nothing
showtable(raw_adf)

# %%
# No run was done with the faithful algorithm.
@assert all(iszero, raw_adf.faithful)
# No pricing was used.
@assert all(==("none"), raw_adf.pricing)
# All runs used round2disc
@assert all(isone, raw_adf.round2disc)
# Every run with rotation has mirror_plates, and every run without rotation does not have mirror plates.
@assert all(raw_adf.rotation .== raw_adf.mirror_plates)
# Every instance appear at most two times.
@assert all(<=(2), count.(.==(unique(raw_adf.instance_path)), (raw_adf.instance_path,)))
# No interrupted instances. FALSE: many did fail because memory exhaution.
#@assert all(!isnan, raw_adf.build_and_solve_time)
#@assert all(==("BUILT_MODEL"), raw_adf.build_stop_reason)
# The value of the objective function and of the extracted solution match. 
@assert all(isapprox(replace(raw_adf.obj_value, NaN => -1.0), raw_adf.solution_value))

# %%
# Create adf (the version of raw_adf) with some pre-treatment
adf = let df = deepcopy(raw_adf)
    df[!, :instance_name] = basename.(df[!, :instance_path])
    select!(df, Not(:stop_reason), :stop_reason => broadwrap(x -> last(split(x, '.'))) => :stop_reason)
    # Bring instance_name to be the first column
    select(df, :instance_name, Not(:instance_name))

    @assert all((df.qt_pevars_after_preprocess .== df.length_pe_after_pricing) .| (df.length_pe_after_pricing .== -1))
    @assert all((df.qt_pevars_after_preprocess .== df.qt_pevars_after_purge) .| (df.qt_pevars_after_purge .== -1))
    select!(df, Not([:qt_pevars_after_purge, :length_pe_after_pricing]))
    rename!(df, :qt_pevars_after_preprocess => :qt_pe_vars)

    @assert all((df.qt_cmvars_after_preprocess .== df.length_cm_after_pricing) .| (df.length_cm_after_pricing .== -1))
    @assert all((df.qt_cmvars_after_preprocess .== df.qt_cmvars_after_purge) .| (df.qt_cmvars_after_purge .== -1))
    select!(df, Not([:qt_cmvars_after_purge, :length_cm_after_pricing]))
    rename!(df, :qt_cmvars_after_preprocess => :qt_cm_vars)

    @assert all((df.qt_plates_after_preprocess .== df.length_pc_after_pricing) .| (df.length_pc_after_pricing .== -1))
    @assert all((df.qt_plates_after_preprocess .== df.qt_plates_after_purge) .| (df.qt_plates_after_purge .== -1))
    select!(df, Not([:qt_plates_after_purge, :length_pc_after_pricing]))
    rename!(df, :qt_plates_after_preprocess => :qt_plates)

    df
end

showtable(sort(adf, :instance_name))

# %%
# Used to find the files for specific instances, no changes done
let df = deepcopy(adf)
    select(filter(:instance_name => in(("A8",)), df), :instance_name, :rotation, :this_data_file)
end

# %%
nonunique(x) = [k for (k, v) in countmap(x) if v > 1]

# These results are a little more messy, let us try to make sense of all these runs.
# No change is done here.
cadf = let df = deepcopy(adf)
    # cleaned df (some runs were repeated by mistake, any of the runs is valid, but we choose
    # the last one in this case)
    primary_keys = [:instance_name, :rotation, :mirror_plates]
    cdf = combine(groupby(df, primary_keys), names(df) .=> last .=> names(df))

    names_df = DataFrame(:instance_name => ais.instance_name)

    no_rotation_runs_df = copy(names_df)
    no_rotation_runs_df[!, :rotation] = zeros(Int, nrow(names_df))
    no_rotation_runs_df[!, :mirror_plates] .= zeros(Int, nrow(names_df))
    rotation_runs_df = copy(names_df)
    rotation_runs_df[!, :rotation] .= ones(Int, nrow(names_df))
    rotation_runs_df[!, :mirror_plates] .= ones(Int, nrow(names_df))

    expected_rows_df = vcat(no_rotation_runs_df, rotation_runs_df)

    # complete cdf
    ccdf = leftjoin(expected_rows_df, cdf; on = primary_keys)

    @show(filter(((k, v),) -> v > 1, countmap(filter(:rotation => isone, ccdf).instance_name)))
    
    ccdf
end
showtable(cadf)

# %%
# Gets some statistics comparing rotated and unrotated variants.
let df = deepcopy(cadf)
    nrdf = filter(:rotation => iszero, df) # no rotation
    rmdf = filter([:rotation, :mirror_plates] => (r, m) -> isone(r) && isone(m), df) # rotation and mirror

    select!.((nrdf, rmdf), ([:instance_name, :build_and_solve_time, :solution_value, :num_nonzeros, :stop_reason],))
    rename!(nrdf,
        :build_and_solve_time => :nr_time, :solution_value => :nr_profit, :num_nonzeros => :nr_nonzeros,
        :stop_reason => :nr_stop_reason
    )
    rename!(rmdf,
        :build_and_solve_time => :rm_time, :solution_value => :rm_profit, :num_nonzeros => :rm_nonzeros,
        :stop_reason => :rm_stop_reason
    )

    jdf = join(nrdf, rmdf; on = :instance_name)

    jdf[!, :nz_growth] = jdf[!, :rm_nonzeros] ./ jdf[!, :nr_nonzeros]
    jdf[!, :nz_abs_growth] = jdf[!, :rm_nonzeros] .- jdf[!, :nr_nonzeros]
    
    # Get the differences in number of nonzeros between rotation and no rotation variants.
    showtable(sort!(select(jdf, :instance_name, :nz_growth, :nz_abs_growth, :nr_nonzeros, :rm_nonzeros), :nz_abs_growth; rev = true))

    #=
    # Profit growth by allowing rotation.
    jdf[!, :p_growth] = ((jdf[!, :rm_profit] ./ jdf[!, :nr_profit]) .- 1.0) .* 100
    showtable(sort!(
            select!(
                filter([:nr_stop_reason, :rm_stop_reason] => (nr, rm) -> !ismissing(nr) && !ismissing(rm) && nr == "OPTIMAL" && rm == "OPTIMAL", jdf),
            :instance_name, :nr_profit, :rm_profit, :p_growth),
    :p_growth; rev = true))

    # Order instances by qt_piece_types
    jjdf = join(jdf, ais; on = :instance_name)
    showtable(sort!(
            select!(
                jjdf,
            vcat(names(ais), ["nr_time", "rm_time"])),
    :qt_piece_types; rev = true))

    # Order instances by how many times more piece area there is compared
    # to original plate area.
    jjdf[!, :pi_area_div_pl_area] = (jjdf[!, :total_piece_area]) ./
        (jjdf[!, :qt_original_plates] .* jjdf[!, :L] .* jjdf[!, :W])
    showtable(sort!(
        select!(
            jjdf,
        :instance_name, :pi_area_div_pl_area, :nr_time, :rm_time),
    :pi_area_div_pl_area; rev = true))
    =#
end

# %%
# Plot nonzeros to time BUT only exclude instances without nonzeros
nonzeros_cadf = let df = deepcopy(cadf)
    df = filter(:stop_reason => c -> !ismissing(c) && c in("OPTIMAL", "TIME_LIMIT"), df)
    df[!, :rotation] = ifelse.(isone.(df[!, :rotation]), "yes", "no")
    df
end

plot(
    nonzeros_cadf, x = "num_nonzeros", y = "build_and_solve_time", color="rotation",
    Geom.point,  Scale.x_log10, Scale.y_log10
)


# %%
# Plot (total piece area)/(stock plate area) to time
let df = deepcopy(cadf)
    jdf = leftjoin(df, ais; on = :instance_name)
    jdf[!, :pi_area_div_pl_area] = jdf[!, :total_piece_area] ./ (jdf[!, :L] .* jdf[!, :W])
    jdf = filter(:stop_reason => c -> !ismissing(c) && c in("OPTIMAL", "TIME_LIMIT"), jdf)
    jdf[!, :rotation] = ifelse.(isone.(jdf[!, :rotation]), "yes", "no")
    p2 = plot(jdf, x = "pi_area_div_pl_area", y = "build_and_solve_time", color="rotation", Geom.point, Scale.y_log10);
end

# %%
# Plot percentage of barrier time to time BUT exclude hard fail instances
per_barrier_time_cadf = let df = deepcopy(cadf)
    df = filter(:stop_reason => c -> !ismissing(c) && c in("OPTIMAL", "TIME_LIMIT"), df)
    @assert all(!isnan, df[!, :root_node_time])
    df[!, :percent_barrier] = df[!, :root_node_time] ./ df[!, :build_and_solve_time]
    df
end

plot(per_barrier_time_cadf, x = "percent_barrier", color="stop_reason", Geom.histogram(;bincount = 10))

# %%
# Plot percentage of barrier time to time (not histogram but density) BUT exclude hard fail instances
per_barrier_time_cadf = let df = deepcopy(cadf)
    all_rows = nrow(df)
    df = filter(:stop_reason => c -> !ismissing(c) && c in("OPTIMAL", "TIME_LIMIT"), df)
    @assert all(!isnan, df[!, :root_node_time])
    df[!, :percent_barrier] = df[!, :root_node_time] ./ df[!, :build_and_solve_time]
    df[!, :bin_percent_barrier] = trunc.(df[!, :percent_barrier]; digits = 1) .+ 0.05
    sdf = combine((:instance_name => length => :qt_instances), groupby(df, [:bin_percent_barrier, :stop_reason]))
    sdf[!, :percent_instances] = (sdf[!, :qt_instances] ./ all_rows) .* 100
    sdf
end

plot(per_barrier_time_cadf, x = "bin_percent_barrier", y = "percent_instances", color = "stop_reason", Geom.bar())

# %%
# Plot pergentage rotation improvement bars
cadf_per_rot_imp_df = let df = deepcopy(cadf)
    filter!(:stop_reason => isequal("OPTIMAL"), df)
    nrdf = filter(:rotation => iszero, df) # no rotation
    rmdf = filter([:rotation, :mirror_plates] => (r, m) -> isone(r) && isone(m), df) # rotation and mirror

    select!.((nrdf, rmdf), ([:instance_name, :solution_value],))
    rename!(nrdf, :solution_value => :nr_profit)
    rename!(rmdf, :solution_value => :rm_profit)
    jdf = join(nrdf, rmdf; on = :instance_name)
    jdf[!, :rot_gap] = ((jdf[!, :nr_profit] ./ jdf[!, :rm_profit]) .- 1.0) .* 100
    jdf = leftjoin(jdf, ais; on = "instance_name")
    jdf[!, :avg_piece_per_of_plate_area] = ((jdf[!, :total_piece_area] ./ jdf[!, :qt_pieces]) ./ (jdf[!, :L] .* jdf[!, :W])) .* 100
    jdf[!, "Instance Name (how much % of stock plate the avg piece cover)"] = jdf[!, :instance_name] .* " (" .* number2latex.(jdf[!, :avg_piece_per_of_plate_area]; enclose = false) .* ")"
    sort!(jdf, :avg_piece_per_of_plate_area)
end

plot(cadf_per_rot_imp_df, x = "Instance Name (how much % of stock plate the avg piece cover)", y = "rot_gap", Guide.ylabel("Improvement (%) by rotation"), Geom.bar())



# %%
# get the files from the instances finished by timeout
let df = deepcopy(cadf)
    filter!(:stop_reason => (sr -> !ismissing(sr) && sr == "TIME_LIMIT"), df)
    foreach(s -> print("./$s "), chop.(df.this_data_file; head = length("./experiments_outputs/")))
end

# %%
# Prints a medium table with all runs (latex format). The difference from the long table
# is that there is a single row per instance, instead of two (which were a single row
# for each configuration). To do this we report nonzero amounts, instead of variables and plates.
let df = deepcopy(cadf)
    select!(df,
        :instance_name, :rotation, :mirror_plates, :num_nonzeros, :solution_value,
        :stop_reason, :build_and_solve_time
    )

    # Change the numbers to latex formatted strings.
    df[!, :solution_value] = number2latex.(df.solution_value; enclose = false)
    df[!, :num_nonzeros] = number2latex.(df.num_nonzeros; enclose = false)
    df[!, :build_and_solve_time] = number2latex.(df.build_and_solve_time; enclose = false)
    # Change the numbers (already strings) to indicators of time limit, or memory limit.
    replace!(df[!, :stop_reason], missing => "REASON_NOT_FOUND")
    replace!(df[!, :solution_value], "-1" => "--")
    @show df[!, :stop_reason]
    df[!, :build_and_solve_time] .= ifelse.(
        df[!, :stop_reason] .== "TIME_LIMIT", ">3600", df[!, :build_and_solve_time]
    )
    df[!, :build_and_solve_time] .= ifelse.(
        df[!, :stop_reason] .== "MEMORY_LIMIT", "--"#="MEM\\textsuperscript{1}"=#, df[!, :build_and_solve_time]
    )
    df[!, :build_and_solve_time] .= ifelse.(
        df[!, :stop_reason] .== "REASON_NOT_FOUND", "--"#="MEM\\textsuperscript{2}"=#, df[!, :build_and_solve_time]
    )
    
    # Escape the name of the instances in latex (they have underlines), and change the nonzeros
    # to show "\ditto" if it is not the first variant of the instance (all variants of the same
    # instance have the same number of nonzeros).
    df[!, :instance_name] .= esc_latex.(df[!, :instance_name])

    nrdf = filter(:rotation => iszero, df) # no rotation
    rmdf = filter([:rotation, :mirror_plates] => (r, m) -> isone(r) && isone(m), df) # rotation and mirror

    select!.((nrdf, rmdf), ([:instance_name, :build_and_solve_time, :solution_value, :num_nonzeros],))
    rename!(nrdf,
        :build_and_solve_time => :nr_time, :solution_value => :nr_profit, :num_nonzeros => :nr_nonzeros
    )
    rename!(rmdf,
        :build_and_solve_time => :rm_time, :solution_value => :rm_profit, :num_nonzeros => :rm_nonzeros
    )

    jdf = join(nrdf, rmdf; on = :instance_name)
    
    # Order the instances by the ones that only have one variant first, then two variants, and then
    # three variants.
    jdf[!, :instance_number] = parse.((Int,), (s -> s[2:end]).(jdf[!, :instance_name]))
    sort!(jdf, :instance_number)

    # Remove instances for which no interesting result is presented.
    filter!([:nr_nonzeros, :rm_nonzeros] => ((nr, rm) -> nr != "-1" || rm != "--"), jdf)

    final_headers = [
        :instance_name => "Inst.",
        :nr_profit     => "G. F",
        :nr_nonzeros   => "G. #nz",
        :nr_time       => "G. T (s)",
        :rm_profit     => "G. R. F",
        :rm_nonzeros   => "G. R. M. #nz",
        :rm_time       => "G. R. M. (s)"
    ]

    select!(jdf, first.(final_headers)) # This is done to re-order the columns.
    rename!(jdf, final_headers...)

    rename!(jdf, names(jdf) .=> esc_latex.(names(jdf)))

    pretty_table(
        jdf; backend = :latex, nosubheader = true, alignment = vcat([:l], repeat([:r], ncol(jdf) - 1))
    )
end

# %%
raw_cldf_path = "./data/G2CSP_CLASS.csv"
raw_cldf = DataFrame(CSV.File(raw_cldf_path))
#nothing
showtable(sort(raw_cldf, :instance_path))

# %%
clis_path = "./data/G2CSP_CLASS_instances_stats.csv"
clis = DataFrame(CSV.File(clis_path))
rename!(clis, :filename => :instance_name)
#nothing
showtable(clis)

# %%
# No run was done with the faithful algorithm.
@assert all(iszero, raw_cldf.faithful)
# No pricing was used.
@assert all(==("none"), raw_cldf.pricing)
# All runs used round2disc
@assert all(isone, raw_cldf.round2disc)
# Every run with rotation has mirror_plates, and every run without rotation does not have mirror plates.
@assert all(raw_cldf.rotation .== raw_cldf.mirror_plates)
# Every instance appear just two times.
@assert all(==(2), count.(.==(unique(raw_cldf.instance_path)), (raw_cldf.instance_path,)))
@assert nrow(raw_cldf) == 1000
@assert sum(raw_cldf.rotation) == 500
@assert all(!isnan, raw_cldf.build_and_solve_time)
@assert all(==("BUILT_MODEL"), raw_cldf.build_stop_reason)
#@assert all(isapprox.(raw_cldf.obj_value, raw_cldf.solution_value))
@assert all(
    (((obj, sol),) -> (isnan(obj) && sol == -1) || isapprox(obj, sol; atol = 1e-6)),
    collect(zip(raw_cldf.obj_value, raw_cldf.solution_value))
)

# %%
# Compute how much of total time was spent on enumeration and purge.
let df = deepcopy(raw_cldf)
    total_time = sum(df.build_and_solve_time)
    enumeration_percent = sum(df.enumeration_time) / total_time
    purge_percent = sum(df.purge_unreachable_time) / total_time
    enumeration_largest_percent = maximum(df.purge_unreachable_time ./ df.build_and_solve_time)
    purge_largest_percent = maximum(df.purge_unreachable_time ./ df.build_and_solve_time)
    @show enumeration_largest_percent
    @show purge_largest_percent
    @show enumeration_percent
    @show purge_percent
    nothing
end

# %%
# Create cdf (the version of raw_cdf) with some pre-treatment
cldf = let df = deepcopy(raw_cldf)
    df[!, :instance_name] = basename.(df[!, :instance_path])
    select!(df, Not(:stop_reason), :stop_reason => broadwrap(x -> last(split(x, '.'))) => :stop_reason)
    
    @assert all(df.qt_pevars_after_preprocess .== df.length_pe_after_pricing)
    @assert all(df.qt_pevars_after_preprocess .== df.qt_pevars_after_purge)
    select!(df, Not([:qt_pevars_after_preprocess, :length_pe_after_pricing]))
    rename!(df, :qt_pevars_after_purge => :qt_pe_vars)
    
    @assert all(df.qt_cmvars_after_preprocess .== df.length_cm_after_pricing)
    @assert all(df.qt_cmvars_after_preprocess .== df.qt_cmvars_after_purge)
    select!(df, Not([:qt_cmvars_after_preprocess, :length_cm_after_pricing]))
    rename!(df, :qt_cmvars_after_purge => :qt_cm_vars)
    
    @assert all(df.qt_plates_after_preprocess .== df.length_pc_after_pricing)
    @assert all(df.qt_plates_after_preprocess .== df.qt_plates_after_purge)
    select!(df, Not([:qt_plates_after_preprocess, :length_pc_after_pricing]))
    rename!(df, :qt_plates_after_purge => :qt_plates)
    
    df
end

# %%
# Used to find the files for specific instances, no changes done
let df = deepcopy(cldf)
    select(filter(:instance_name => in(("cl_09_020_08",)), df), :instance_name, :rotation, :this_data_file)
end

# %%
# Plot nonzeros to time
nonzeros_cldf = let df = deepcopy(cldf)
    df = filter(:stop_reason => c -> !ismissing(c) && c in("OPTIMAL", "TIME_LIMIT"), df)
    df[!, :rotation] = ifelse.(isone.(df[!, :rotation]), "yes", "no")
    df
end

plot(
    nonzeros_cldf, x = "num_nonzeros", y = "build_and_solve_time", color="rotation",
    Geom.point,  Scale.x_log10, Scale.y_log10
)


# %%
let df1 = deepcopy(nonzeros_cadf), df2 = deepcopy(nonzeros_cldf)
    jdf = vcat(df1, df2)
    jdf[!, "Dataset"] = ifelse.(occursin.(("cl_",), jdf[!, :instance_name]), "CLASS", "A")
    plot(
        jdf, x = "num_nonzeros", y = "build_and_solve_time",
        Scale.x_log10, Scale.y_log10, Guide.xlabel("Time (s)"), Guide.ylabel("Number of nonzeros"),
        color = "Dataset",
    )
end

# %%
# Plot (total piece area)/(stock plate area) to time
let df = deepcopy(cldf)
    jdf = leftjoin(df, clis; on = :instance_name)
    jdf[!, :pi_area_div_pl_area] = jdf[!, :total_piece_area] ./ (jdf[!, :L] .* jdf[!, :W])
    jdf = filter(:stop_reason => c -> !ismissing(c) && c in("OPTIMAL", "TIME_LIMIT"), jdf)
    jdf[!, :rotation] = ifelse.(isone.(jdf[!, :rotation]), "yes", "no")
    p2 = plot(jdf, x = "pi_area_div_pl_area", y = "build_and_solve_time", color="rotation", Geom.point, Scale.y_log10);
end

# %%
# Plot percentage of barrier time to time BUT exclude hard fail instances
per_barrier_time_cldf = let df = deepcopy(cldf)
    df = filter(:stop_reason => c -> !ismissing(c) && c in("OPTIMAL", "TIME_LIMIT"), df)
    @assert all(!isnan, df[!, :root_node_time])
    df[!, :percent_barrier] = df[!, :root_node_time] ./ df[!, :build_and_solve_time]
    df[!, :percent_barrier] = ifelse.(
        df[!, :percent_barrier] .> 1.0,
        ifelse.(df[!, :root_node_time] .- df[!, :build_and_solve_time] .<= 0.001, 1.0, NaN),
        df[!, :percent_barrier]
    )
    #@show df[findfirst(>(1.0), df[!, :percent_barrier]), :]
    df
end
plot(per_barrier_time_cldf, x = "percent_barrier", color="stop_reason", Geom.histogram(;bincount = 10))


# %%
# Plot percentage of barrier time to time (not histogram but density) BUT exclude hard fail instances
per_barrier_time_cldf = let df = deepcopy(cldf)
    all_rows = nrow(df)
    df = filter(:stop_reason => c -> !ismissing(c) && c in("OPTIMAL", "TIME_LIMIT"), df)
    @assert all(!isnan, df[!, :root_node_time])
    df[!, :percent_barrier] = df[!, :root_node_time] ./ df[!, :build_and_solve_time]
    df[!, :percent_barrier] = ifelse.(
        df[!, :percent_barrier] .> 1.0,
        ifelse.(df[!, :root_node_time] .- df[!, :build_and_solve_time] .<= 0.001, 0.95, NaN),
        df[!, :percent_barrier]
    )
    df[!, :bin_percent_barrier] = trunc.(df[!, :percent_barrier]; digits = 1) .+ 0.05
    sdf = combine((:instance_name => length => :qt_instances), groupby(df, [:bin_percent_barrier, :stop_reason]))
    sdf[!, :percent_instances] = (sdf[!, :qt_instances] ./ all_rows) .* 100
    sdf
end

plot(
    per_barrier_time_cldf,
    x = "bin_percent_barrier", y = "percent_instances", color = "stop_reason",
    Geom.bar()
)

# %%
# TODO: number of instances is too distinct, need to free x axis
let df1 = deepcopy(per_barrier_time_cadf), df2 = deepcopy(per_barrier_time_cldf)
    df1[:dataset] = "A"
    df2[:dataset] = "CLASS"
    jdf = vcat(df1, df2)
    replace!(jdf[!, :stop_reason], "OPTIMAL" => "Solved")
    replace!(jdf[!, :stop_reason], "TIME_LIMIT" => "Timeout")
    # Change from ratio to percentage.
    jdf[!, :bin_percent_barrier] .= jdf[!, :bin_percent_barrier] .* 100
    # The sort is just a stupid workaround to guarantee that the timeout bars
    # will be stacked over the solved bars.
    sort!(jdf, :stop_reason; rev = true)
    rename!(jdf, "stop_reason" => "Stop Reason")
    plot(
        jdf, xgroup="dataset", x = "bin_percent_barrier", y = "percent_instances",
        Guide.xlabel("Barrier time percentage (10% bins) by dataset"),
        Guide.ylabel("Percentage of dataset instances"),
        color = "Stop Reason", Geom.subplot_grid(Geom.bar()),
        Gadfly.Scale.color_discrete_manual("gray50", "gray10"),
        Guide.colorkey(; pos = [0.76, -0.3])
    ) |> PDF("../plots/g2csp_a_class_histogram_barrier_percent_time.pdf")
end

# %%
# Plot pergentage rotation improvement bars
cldf_per_rot_imp_df = let df = deepcopy(cldf)
    filter!(:stop_reason => isequal("OPTIMAL"), df)
    nrdf = filter(:rotation => iszero, df) # no rotation
    rmdf = filter([:rotation, :mirror_plates] => (r, m) -> isone(r) && isone(m), df) # rotation and mirror

    select!.((nrdf, rmdf), ([:instance_name, :solution_value],))
    rename!(nrdf, :solution_value => :nr_profit)
    rename!(rmdf, :solution_value => :rm_profit)
    jdf = join(nrdf, rmdf; on = :instance_name)
    jdf[!, :rot_gap] = ((jdf[!, :nr_profit] ./ jdf[!, :rm_profit]) .- 1.0) .* 100
    jdf = leftjoin(jdf, clis; on = "instance_name")
    jdf[!, :avg_piece_per_of_plate_area] = ((jdf[!, :total_piece_area] ./ jdf[!, :qt_pieces]) ./ (jdf[!, :L] .* jdf[!, :W])) .* 100
    jdf[!, "Instance Name (how much % of stock plate the avg piece cover)"] = jdf[!, :instance_name] .* " (" .* number2latex.(jdf[!, :avg_piece_per_of_plate_area]; enclose = false) .* ")"
    sort!(jdf, :avg_piece_per_of_plate_area)
end

plot(cldf_per_rot_imp_df, x = "Instance Name (how much % of stock plate the avg piece cover)", y = "rot_gap", Guide.ylabel("Improvement (%) by rotation"), Geom.bar())

# %%
# Plot pergentage rotation improvement bars
# WRITE ABOUT: probably ordering by waste (%) in the solution without rotation
# would give a better correlation, but this assumes we have this information
# beforehand and does not helps with a priori guessing how much a instance will
# benefit from allowing rotation
let df1 = cadf_per_rot_imp_df, df2 = deepcopy(cldf_per_rot_imp_df)
    jdf = vcat(df1, df2)
    jdf[!, "Dataset"] = ifelse.(occursin.(("cl_",), jdf[!, :instance_name]), "CLASS", "A")
    sort!(jdf, :avg_piece_per_of_plate_area)
    plot(jdf,
        x = "Instance Name (how much % of stock plate the avg piece cover)", y = "rot_gap",
        Guide.ylabel("Improvement (%) by rotation"), color = "Dataset", Geom.bar()
    )
end


# %%
# Gets some statistics comparing rotated and unrotated variants.
let df = deepcopy(cldf)
    nrdf = filter(:rotation => iszero, df) # no rotation
    rmdf = filter([:rotation, :mirror_plates] => (r, m) -> isone(r) && isone(m), df) # rotation and mirror

    select!.((nrdf, rmdf), ([:instance_name, :build_and_solve_time, :solution_value, :num_nonzeros, :stop_reason],))
    rename!(nrdf,
        :build_and_solve_time => :nr_time, :solution_value => :nr_profit, :num_nonzeros => :nr_nonzeros,
        :stop_reason => :nr_stop_reason
    )
    rename!(rmdf,
        :build_and_solve_time => :rm_time, :solution_value => :rm_profit, :num_nonzeros => :rm_nonzeros,
        :stop_reason => :rm_stop_reason
    )

    jdf = join(nrdf, rmdf; on = :instance_name)

    jdf[!, :nz_growth] = jdf[!, :rm_nonzeros] ./ jdf[!, :nr_nonzeros]
    jdf[!, :nz_abs_growth] = jdf[!, :rm_nonzeros] .- jdf[!, :nr_nonzeros]
    
    # Get the differences in number of nonzeros between rotation and no rotation variants.
    showtable(sort!(select(jdf, :instance_name, :nz_growth, :nz_abs_growth, :nr_nonzeros, :rm_nonzeros), :nz_growth; rev = true))

    # Profit growth by allowing rotation.
    jdf[!, :p_growth] = ((jdf[!, :rm_profit] ./ jdf[!, :nr_profit]) .- 1.0) .* 100
    showtable(sort!(
        select!(
            filter([:nr_stop_reason, :rm_stop_reason] => (nr, rm) -> !ismissing(nr) && !ismissing(rm) && nr == "OPTIMAL" && rm == "OPTIMAL", jdf),
        :instance_name, :nr_profit, :rm_profit, :p_growth),
    :p_growth; rev = true))
    
    jjdf = join(jdf, clis; on = :instance_name)
    jjdf[!, :pi_area_div_pl_area] = jjdf[!, :total_piece_area] ./ (jjdf[!, :L] .* jjdf[!, :W])

    showtable(sort!(
        select!(
            jjdf,
        :instance_name, :pi_area_div_pl_area, :nr_time, :rm_time),
    :pi_area_div_pl_area; rev = true))
    
end

# %%
# Sum all solutions from all instances optimally solved in both datasets.
let df = deepcopy(cldf)
    select!(df, :instance_name, :rotation, :stop_reason, :solution_value)
    filter!(:stop_reason => ==("OPTIMAL"), df)
    select!(df, Not(:stop_reason))

    nrdf = select!(filter(:rotation => iszero, df), Not(:rotation)) # no rotation
    @show nrow(nrdf)
    rmdf = select!(filter(:rotation => isone, df), Not(:rotation)) # rotation and mirror

    rename!(nrdf, :solution_value => :nr_solution)
    rename!(rmdf, :solution_value => :rm_solution)
    jdf = innerjoin(nrdf, rmdf; on = :instance_name)
    
    nrow(jdf)
    @show sum(jdf.nr_solution)
    @show sum(jdf.rm_solution)
    @show (1.0 - sum(jdf.rm_solution)/sum(jdf.nr_solution)) * 100
end

# %%
# Check if we have memory exhaustion without considering rotation
# Sum all solutions from all instances optimally solved in both datasets.
let df = deepcopy(cldf)
    select!(df, :instance_name, :rotation, :stop_reason, :solution_value)

    nrdf = select!(filter(:rotation => iszero, df), Not(:rotation)) # no rotation
    rmdf = select!(filter(:rotation => isone, df), Not(:rotation)) # rotation and mirror
    
    @show nrow(nrdf)
    @show unique(nrdf.stop_reason)
    @show nrow(rmdf)
    @show unique(rmdf.stop_reason)

    rename!(nrdf, :solution_value => :nr_solution)
    rename!(rmdf, :solution_value => :rm_solution)
    jdf = innerjoin(nrdf, rmdf; on = :instance_name)
    
    nrow(jdf)
    @show sum(jdf.nr_solution)
    @show sum(jdf.rm_solution)
    @show (1.0 - sum(jdf.rm_solution)/sum(jdf.nr_solution)) * 100
end

# %%
# Check the times and solutions of a subset of the CLASS dataset.
let df = deepcopy(cldf)
    filter!(:instance_name => (name -> occursin("cl_10_040", name)), df)
    foreach(s -> print("./$s "), chop.(df.this_data_file; head = length("./experiments_outputs/"), tail = 0))
    #@show names(df)
    select!(df, :instance_name, :rotation, :this_data_file #=:num_nonzeros, :run_total_time, :obj_bound, :solution_value, :stop_reason=#)
    sort!(df, [:instance_name, :rotation])
    showtable(df)
end

# %%
let df = deepcopy(cldf)
    select!(df,
        :instance_name, :rotation, :qt_pe_vars, :qt_cm_vars, :qt_plates,
        :stop_reason, :build_and_solve_time
    )

    df[!, :build_and_solve_time] = number2latex.(df.build_and_solve_time; enclose = false)
    sort!(df, [:instance_name, :rotation])
    replace!(df[!, :stop_reason], "OPTIMAL" => "F", "INFEASIBLE" => "N", "TIME_LIMIT" => "X")
    
    final_headers = [
        :instance_name         => "Inst.",
        :rotation              => "R.",
        :qt_pe_vars            => "#E. vars",
        :qt_cm_vars            => "#C. vars",
        :qt_plates             => "#Plates",
        :stop_reason           => "G. F.",
        :build_and_solve_time  => "Time (s)"
    ]

    select!(df, first.(final_headers)) # This is done to re-order the columns.
    rename!(df, final_headers...)
    rename!(df, names(df) .=> esc_latex.(names(df)))

    pretty_table(
        df; backend = :latex, nosubheader = true, alignment = vcat([:l], repeat([:r], ncol(df) - 1))
    )
end

# %%
# Prints a medium table with all runs (latex format). The difference from the long table
# is that there is a single row per instance, instead of two (which were a single row
# for each configuration), and that we summarise the results for 10 instance groups
# (what leads to not presenting "solution value" but instead, we present number of
# solved instances)
let df = deepcopy(cldf)
    select!(df,
        :instance_name, :rotation, :mirror_plates, :num_nonzeros, :solution_value,
        :stop_reason, :build_and_solve_time
    )
    df[!, :instance_group] = join.((a -> a[begin:(end-1)]).(split.(df[!, :instance_name], ('_',))), ('_',))
    # NOTE, if I use mean (instead of sum), then I need to use the qt_solved or filter to have the
    # correct denominator
    sdf = combine(
        groupby(df, [:instance_group, :rotation]),
        :stop_reason => (x -> sum("OPTIMAL" .== x)) => :qt_solved,
        [:stop_reason, :num_nonzeros] => 
            ((sr, nnz) -> (sum(ifelse.(sr .== "OPTIMAL", nnz, 0)) / sum(==("OPTIMAL"), sr))) => :sum_built_nnz,
        [:stop_reason, :build_and_solve_time] => 
            ((sr, bst) -> (sum(ifelse.(sr .== "OPTIMAL", bst, 0)) / sum(==("OPTIMAL"), sr))) => :sum_solve_time,
    )

    # We need to convert NaN to -1 and "-1" to "--" otherwise we cannot apply round in the meantime.
    replace!(sdf[!, :sum_built_nnz], NaN => -1)
    replace!(sdf[!, :sum_solve_time], NaN => -1)
    
    sdf[!, :sum_built_nnz] = number2latex.(round.((Int,), sdf.sum_built_nnz); enclose = false)
    @show sdf[!, :sum_built_nnz]
    sdf[!, :sum_solve_time] = number2latex.(sdf.sum_solve_time; enclose = false)
    sdf[!, :instance_group] .= esc_latex.(sdf[!, :instance_group])

    replace!(sdf[!, :sum_built_nnz], "-1" => "--")
    replace!(sdf[!, :sum_solve_time], "-1.00" => "--")

    nrdf = filter(:rotation => iszero, sdf) # no rotation
    rmdf = filter(:rotation => isone, sdf) # rotation and mirror

    select!.((nrdf, rmdf), ([:instance_group, :qt_solved, :sum_built_nnz, :sum_solve_time],))
    rename!(nrdf,
        :sum_solve_time => :nr_time, :qt_solved => :nr_qt_solved, :sum_built_nnz => :nr_nonzeros
    )
    rename!(rmdf,
        :sum_solve_time => :rm_time, :qt_solved => :rm_qt_solved, :sum_built_nnz => :rm_nonzeros
    )

    jdf = join(nrdf, rmdf; on = :instance_group)
    sort!(jdf, :instance_group)
    
    final_headers = [
        :instance_group => "Group",
        :nr_qt_solved   => "G. F",
        :nr_nonzeros    => "G. #nz",
        :nr_time        => "G. T (s)",
        :rm_qt_solved   => "G. R. F",
        :rm_nonzeros    => "G. R. M. #nz",
        :rm_time        => "G. R. M. (s)"
    ]

    select!(jdf, first.(final_headers)) # This is done to re-order the columns.
    rename!(jdf, final_headers...)

    rename!(jdf, names(jdf) .=> esc_latex.(names(jdf)))

    pretty_table(
        jdf; backend = :latex, nosubheader = true, alignment = vcat([:l], repeat([:r], ncol(jdf) - 1))
    )
end
