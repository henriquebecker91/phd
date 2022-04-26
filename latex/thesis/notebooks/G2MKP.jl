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
#     display_name: Julia 1.5 1.5.4
#     language: julia
#     name: julia-1.5-1.5
# ---

# %%
# Ensure your working directory is: https://github.com/henriquebecker91/phd/tree/master/latex/thesis/notebooks
include("notebook_setup.jl")

# %%
ais_path = "./data/G2MKP_A_instances_stats.csv"
ais = DataFrame(CSV.File(ais_path))
rename!(ais, :filename => :instance_name)
#nothing
showtable(ais)

# %%
raw_adf_path = "./data/G2MKP_A.csv"
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
# Every instance appear just two times. FALSE: we repeated some runs inadvertently.
#@assert all(==(2), count.(.==(unique(raw_adf.instance_path)), (raw_adf.instance_path,)))
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
    select(filter(:instance_name => in(("A13_M2", "A13_M4")), adf), :instance_name, :rotation, :this_data_file)
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

    names_df = DataFrame(CSV.File("./data/G2MKP_instances.csv"))

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
function row_color_MKP(instance_name)
    qt_knapsacks = parse(Int, last(instance_name))
    @assert qt_knapsacks in (2, 4, 8)
    return if qt_knapsacks == 2
        instance_name
    elseif qt_knapsacks == 4
        "\rowcolor{gray-inner-row} $instance_name"
    elseif qt_knapsacks == 8
        "\rowcolor{gray-table-row} $instance_name"
    end
end

function non_zeros_MKP(instance_name, nonzero)
    qt_knapsacks = parse(Int, last(instance_name))
    @assert qt_knapsacks in (2, 4, 8)
    return if qt_knapsacks == 2
        string(nonzero)
    else
        "\\ditto"
    end
end

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
            ["instance_name", "qt_pieces", "nr_time", "rm_time"]),
    :qt_pieces))

    # Order instances by how many times more piece area there is compared
    # to original plate area.
    #=
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
# Plot extra piece area to time
let df = deepcopy(cadf)
    jdf = leftjoin(df, ais; on = :instance_name)
    jdf[!, :pi_area_div_pl_area] = (jdf[!, :total_piece_area]) ./
        (jdf[!, :qt_original_plates] .* jdf[!, :L] .* jdf[!, :W])
    jdf = filter(:stop_reason => c -> !ismissing(c) && c in("OPTIMAL", "TIME_LIMIT"), jdf)
    jdf[!, :rotation] = ifelse.(isone.(jdf[!, :rotation]), "yes", "no")
    p2 = plot(jdf, x = "pi_area_div_pl_area", y = "build_and_solve_time", color="rotation", Geom.point, Scale.y_log10);
end

# %%
# Plot percentage of barrier time to time BUT exclude hard fail instances
per_barrier_time_cadf = let df = deepcopy(cadf)
    df = filter(:stop_reason => c -> !ismissing(c) && c in("OPTIMAL", "TIME_LIMIT"), df)
    replace!(df[!, :stop_reason], "OPTIMAL" => "Solved")
    replace!(df[!, :stop_reason], "TIME_LIMIT" => "Timeout")
    @assert all(!isnan, df[!, :root_node_time])
    df[!, :percent_barrier] = (df[!, :root_node_time] ./ df[!, :build_and_solve_time]) .* 100
    df
end

plot(
    per_barrier_time_cadf, x = "percent_barrier", color="stop_reason", Geom.histogram(; bincount = 10, limits = (; max = 100)),
    Gadfly.Guide.xlabel("Percentage of the total time spent solving barrier (10% bins)"),
    Gadfly.Guide.ylabel("Number of runs/instances"),
    Guide.colorkey(; title = "Stop reason"),
    Gadfly.Scale.color_discrete_manual("gray10", "gray50")
) |> PDF("../plots/g2mkp_a_histogram_barrier_percent_time.pdf")

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
    jdf[!, :rot_gap] = (jdf[!, :rm_profit] ./ jdf[!, :nr_profit] .- 1.0) * 100
    jdf = leftjoin(jdf, ais; on = "instance_name")
    jdf[!, :avg_piece_per_of_plate_area] = ((jdf[!, :total_piece_area] ./ jdf[!, :qt_pieces]) ./ (jdf[!, :L] .* jdf[!, :W])) .* 100
    jdf[!, "Instance Name (how much % of stock plate the avg piece cover)"] = jdf[!, :instance_name] .* " (" .* number2latex.(jdf[!, :avg_piece_per_of_plate_area]; enclose = false) .* ")"
    sort!(jdf, :avg_piece_per_of_plate_area)
end

plot(cadf_per_rot_imp_df, x = "Instance Name (how much % of stock plate the avg piece cover)", y = "rot_gap", Guide.ylabel("Improvement (%) by rotation"), Geom.bar())


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
        df[!, :stop_reason] .== "MEMORY_LIMIT", "MEM\\textsuperscript{1}", df[!, :build_and_solve_time]
    )
    df[!, :build_and_solve_time] .= ifelse.(
        df[!, :stop_reason] .== "REASON_NOT_FOUND", "MEM\\textsuperscript{2}", df[!, :build_and_solve_time]
    )
    
    # Escape the name of the instances in latex (they have underlines), and change the nonzeros
    # to show "\ditto" if it is not the first variant of the instance (all variants of the same
    # instance have the same number of nonzeros).
    df[!, :instance_name] .= esc_latex.(df[!, :instance_name])
    transform!(df, [:instance_name, :num_nonzeros] => ByRow(non_zeros_MKP) => :num_nonzeros)

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
    jdf[!, :parent_instance] = first.(split.(jdf[!, :instance_name], '_'))
    jdf[!, :qt_variants] = count.(.==(jdf[!, :parent_instance]), (jdf[!, :parent_instance],))
    sort!(jdf, [:qt_variants, :instance_name])

    # Remove instances for which no interesting result is presented.
    #filter!([:nr_profit, :rm_profit] => (nr, rm) -> !(nr in ("--", "0") && rm in ("--", "0")), jdf)

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

    replace!(row_color_MKP, jdf[!, "Inst."])

    pretty_table(
        jdf; backend = :latex, nosubheader = true, alignment = vcat([:l], repeat([:r], ncol(jdf) - 1))
    )
end

# %%
raw_cwdf_path = "./data/G2MKP_CW.csv"
raw_cwdf = DataFrame(CSV.File(raw_cwdf_path))
#nothing
showtable(raw_cwdf)

# %%
cwis_path = "./data/G2MKP_CW_instances_stats.csv"
cwis = DataFrame(CSV.File(cwis_path))
rename!(cwis, :filename => :instance_name)
#nothing
showtable(cwis)

# %%
# No run was done with the faithful algorithm.
@assert all(iszero, raw_cwdf.faithful)
# No pricing was used.
@assert all(==("none"), raw_cwdf.pricing)
# All runs used round2disc
@assert all(isone, raw_cwdf.round2disc)
# Every run with rotation has mirror_plates, and every run without rotation does not have mirror plates.
@assert all(raw_cwdf.rotation .== raw_cwdf.mirror_plates)
# Every instance appear just two times.
@assert all(==(2), count.(.==(unique(raw_cwdf.instance_path)), (raw_cwdf.instance_path,)))
@assert nrow(raw_cwdf) == 56
@assert sum(raw_cwdf.rotation) == 28
@assert all(!isnan, raw_cwdf.build_and_solve_time)
@assert all(==("BUILT_MODEL"), raw_cwdf.build_stop_reason)
@assert all(isapprox(raw_cwdf.obj_value, raw_cwdf.solution_value))

# %%
# Compute how much of total time was spent on enumeration and purge.
let df = deepcopy(raw_cwdf)
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
cwdf = let df = deepcopy(raw_cwdf)
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
# Plot nonzeros to time
nonzeros_cwdf = let df = deepcopy(cwdf)
    df = filter(:stop_reason => c -> !ismissing(c) && c in("OPTIMAL", "TIME_LIMIT"), df)
    df[!, :rotation] = ifelse.(isone.(df[!, :rotation]), "yes", "no")
    df
end

plot(
    nonzeros_cwdf, x = "num_nonzeros", y = "build_and_solve_time", color="rotation",
    Geom.point,  Scale.x_log10, Scale.y_log10
)


# %%
let df1 = deepcopy(nonzeros_cadf), df2 = deepcopy(nonzeros_cwdf)
    jdf = vcat(df1, df2)
    jdf[!, "Dataset"] = ifelse.(occursin.(("CW",), jdf[!, :instance_name]), "CW", "A")
    plot(
        jdf, x = "num_nonzeros", y = "build_and_solve_time",
        Scale.x_log10, Scale.y_log10, Guide.xlabel("Number of nonzeros"), Guide.ylabel("Time (s)"),
        color = "Dataset",
    )
end

# %%
# Plot extra piece area to time
let df = deepcopy(cwdf)
    jdf = leftjoin(df, cwis; on = :instance_name)
    jdf[!, :pi_area_div_pl_area] = (jdf[!, :total_piece_area]) ./
        (jdf[!, :qt_original_plates] .* jdf[!, :L] .* jdf[!, :W])
    jdf = filter(:stop_reason => c -> !ismissing(c) && c in("OPTIMAL", "TIME_LIMIT"), jdf)
    jdf[!, :rotation] = ifelse.(isone.(jdf[!, :rotation]), "yes", "no")
    p2 = plot(jdf, x = "pi_area_div_pl_area", y = "build_and_solve_time", color="rotation", Geom.point, Scale.y_log10);
end

# %%
# Plot qt_knapsacks to time
let df = deepcopy(cwdf)
    jdf = leftjoin(df, cwis; on = :instance_name)
    jdf = filter(:stop_reason => c -> !ismissing(c) && c in("OPTIMAL", "TIME_LIMIT"), jdf)
    jdf[!, :rotation] = ifelse.(isone.(jdf[!, :rotation]), "yes", "no")
    p2 = plot(jdf, x = "qt_original_plates", y = "build_and_solve_time", color="rotation", Geom.point, Scale.y_log10);
end

# %%

# %%
# Plot percentage of barrier time to time BUT exclude hard fail instances
per_barrier_time_cwdf = let df = deepcopy(cwdf)
    df = filter(:stop_reason => c -> !ismissing(c) && c in("OPTIMAL", "TIME_LIMIT"), df)
    replace!(df[!, :stop_reason], "OPTIMAL" => "Solved")
    replace!(df[!, :stop_reason], "TIME_LIMIT" => "Timeout")
    @assert all(!isnan, df[!, :root_node_time])
    df[!, :percent_barrier] = (df[!, :root_node_time] ./ df[!, :build_and_solve_time]) .* 100
    df
end

plot(
    per_barrier_time_cwdf, x = "percent_barrier", color="stop_reason", Geom.histogram(; bincount = 10, limits = (; max = 100)),
    Gadfly.Guide.xlabel("Percentage of the total time spent solving barrier (10% bins)"),
    Gadfly.Guide.ylabel("Number of runs/instances"),
    Guide.colorkey(; title = "Stop reason"),
    Gadfly.Scale.color_discrete_manual("gray10", "gray50")
) |> PDF("../plots/g2mkp_cw_histogram_barrier_percent_time.pdf")


# %%
let
    df = DataFrame(
        :percent_barrier => sqrt.(rand(100)) .* 100,
        :stop_reason => ifelse.(rand(100) .> 0.75, "Solved", "Timeout")
    )
    p = plot(
        df, x = "percent_barrier", color="stop_reason", Geom.histogram(;bincount = 10),
        Gadfly.Guide.xlabel("Time spent solving barrier (10% bins)"),
        Gadfly.Guide.ylabel("Amount of instances"),
        Guide.colorkey(; title = "Stop reason")
    )
    p |> PNG("example.png")
end

# %%
let df1 = deepcopy(per_barrier_time_cadf), df2 = deepcopy(per_barrier_time_cwdf)
    jdf = vcat(df1, df2)
    jdf[!, "Dataset"] = ifelse.(occursin.(("CW",), jdf[!, :instance_name]), "CW", "A")
    rename!(jdf, "stop_reason" => "Stop Reason")
    plot(
        jdf, xgroup="Dataset", x = "percent_barrier", Guide.xlabel("Barrier time percentage (10% bins) by dataset"),
        Guide.ylabel("Number of instances"),
        color = "Stop Reason", Geom.subplot_grid(Geom.histogram(;bincount = 10, limits = (xmin = 0.0, xmax = 1.0)))
    )
end

# %%
# Plot pergentage rotation improvement bars
# WRITE ABOUT: probably ordering by waste (%) in the solution without rotation
# would give a better correlation, but this assumes we have this information
# beforehand and does not helps with a priori guessing how much a instance will
# benefit from allowing rotation
cwdf_per_rot_imp_df = let df = deepcopy(cwdf)
    filter!(:stop_reason => isequal("OPTIMAL"), df)
    nrdf = filter(:rotation => iszero, df) # no rotation
    rmdf = filter([:rotation, :mirror_plates] => (r, m) -> isone(r) && isone(m), df) # rotation and mirror

    select!.((nrdf, rmdf), ([:instance_name, :solution_value],))
    rename!(nrdf, :solution_value => :nr_profit)
    rename!(rmdf, :solution_value => :rm_profit)
    jdf = join(nrdf, rmdf; on = :instance_name)
    jdf[!, :rot_gap] = (jdf[!, :rm_profit] ./ jdf[!, :nr_profit] .- 1.0) * 100
    jdf = leftjoin(jdf, cwis; on = "instance_name")
    jdf[!, :avg_piece_per_of_plate_area] = ((jdf[!, :total_piece_area] ./ jdf[!, :qt_pieces]) ./ (jdf[!, :L] .* jdf[!, :W])) .* 100
    jdf[!, "Instance Name (how much % of stock plate the avg piece cover)"] = jdf[!, :instance_name] .* " (" .* number2latex.(jdf[!, :avg_piece_per_of_plate_area]; enclose = false) .* ")"
    sort!(jdf, :avg_piece_per_of_plate_area)
end

plot(cwdf_per_rot_imp_df, x = "Instance Name (how much % of stock plate the avg piece cover)", y = "rot_gap", Guide.ylabel("Improvement (%) by rotation"), Geom.bar())


# %%
# Plot pergentage rotation improvement bars
# WRITE ABOUT: probably ordering by waste (%) in the solution without rotation
# would give a better correlation, but this assumes we have this information
# beforehand and does not helps with a priori guessing how much a instance will
# benefit from allowing rotation]
# IMPORTANT NOTE: complete bogus ordering, the number of pieces per plate means nothing
# because different instances have different amounts of unused piece area (ufulfilled
# demand), ideally we would need to have the average piece area by the stock plate area.
let df1 = cadf_per_rot_imp_df, df2 = deepcopy(cwdf_per_rot_imp_df)
    jdf = vcat(df1, df2)
    jdf[!, "Dataset"] = ifelse.(occursin.(("CW",), jdf[!, :instance_name]), "CW", "A")
    sort!(jdf, :avg_piece_per_of_plate_area)
    plot(jdf,
        x = "Instance Name (how much % of stock plate the avg piece cover)", y = "rot_gap",
        Guide.ylabel("Improvement (%) by rotation"), color = "Dataset", Geom.bar()
    )
end


# %%
# Gets some statistics comparing rotated and unrotated variants.
let df = deepcopy(cwdf)
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
    
    jjdf = join(jdf, cwis; on = :instance_name)
    jjdf[!, :pi_area_div_pl_area] = (jjdf[!, :total_piece_area]) ./
        (jjdf[!, :qt_original_plates] .* jjdf[!, :L] .* jjdf[!, :W])

    showtable(sort!(
        select!(
            jjdf,
        :instance_name, :pi_area_div_pl_area, :nr_time, :rm_time),
    :pi_area_div_pl_area; rev = true))
    
end

# %%
let df = deepcopy(cwdf)
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
function row_color_MKP(instance_name)
    qt_knapsacks = parse(Int, last(instance_name))
    @assert qt_knapsacks in (2, 4, 8)
    return if qt_knapsacks == 2
        instance_name
    elseif qt_knapsacks == 4
        "\rowcolor{gray-inner-row} $instance_name"
    elseif qt_knapsacks == 8
        "\rowcolor{gray-table-row} $instance_name"
    end
end

function non_zeros_MKP(instance_name, nonzero)
    qt_knapsacks = parse(Int, last(instance_name))
    @assert qt_knapsacks in (2, 4, 8)
    return if qt_knapsacks == 2
        string(nonzero)
    else
        "\\ditto"
    end
end


# Prints a medium table with all runs (latex format). The difference from the long table
# is that there is a single row per instance, instead of two (which were a single row
# for each configuration). To do this we report nonzero amounts, instead of variables and plates.
let df = deepcopy(cwdf)
    select!(df,
        :instance_name, :rotation, :mirror_plates, :num_nonzeros, :solution_value,
        :stop_reason, :build_and_solve_time
    )

    df[!, :num_nonzeros] = number2latex.(df.num_nonzeros; enclose = false)
    df[!, :build_and_solve_time] = number2latex.(df.build_and_solve_time; enclose = false)
    df[!, :build_and_solve_time] .= ifelse.(
        df[!, :stop_reason] .== "TIME_LIMIT", ">3600", df[!, :build_and_solve_time]
    )
    df[!, :instance_name] .= esc_latex.(df[!, :instance_name])
    transform!(df, [:instance_name, :num_nonzeros] => ByRow(non_zeros_MKP) => :num_nonzeros)

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
    sort!(jdf, :instance_name)
    
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

    replace!(row_color_MKP, jdf[!, "Inst."])

    pretty_table(
        jdf; backend = :latex, nosubheader = true, alignment = vcat([:l], repeat([:r], ncol(jdf) - 1))
    )
end

# %%
let df = deepcopy(cwdf)
    select!(df, :instance_name, :rotation, :this_data_file)
    interesting_cases = [
        ("CW01_M2", 0),
        ("CW01_M4", 0),
        ("CW03_M2", 1),
        ("CW03_M4", 1),
        ("CW05_M2", 1),
        ("CW05_M4", 1),
        ("CW08_M2", 0),
        ("CW08_M4", 0),
        ("CW09_M2", 0),
        ("CW09_M2", 1),
        ("CW09_M4", 0),
        ("CW11_M2", 0),
        ("CW11_M4", 0),
        ("CW11_M8", 0),
    ]
    filter!([:instance_name, :rotation] => (name, r) -> (name, r) in interesting_cases, df)
    sort!(df, [:instance_name, :rotation])
    println(collect(zip(first.(interesting_cases), basename.(df[!, :this_data_file]))))
    foreach(x -> print(join(split(x, '/')[3:end], '/'), " "), df[!, :this_data_file])
    showtable(df)
end

# %%
let df = deepcopy(adf)
    select!(df, :instance_name, :rotation, :this_data_file)
    interesting_cases = [
        ("A13_M2", 0),
        ("A13_M4", 0),
    ]
    filter!([:instance_name, :rotation] => (name, r) -> (name, r) in interesting_cases, df)
    sort!(df, [:instance_name, :rotation])
    println(collect(zip(first.(interesting_cases), basename.(df[!, :this_data_file]))))
    foreach(x -> print(join(split(x, '/')[3:end], '/'), " "), df[!, :this_data_file])
    showtable(df)
end
