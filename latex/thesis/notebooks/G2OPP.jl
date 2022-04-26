# ---
# jupyter:
#   jupytext:
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
raw_cdf_path = "./data/G2OPP_hopper_c.csv"
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
# Every run with rotation has mirror_plates, and every run without rotation does not have mirror plates.
@assert all(raw_cdf.rotation .== raw_cdf.mirror_plates)
# Every instance appear just two times.
@assert all(==(2), count.(.==(unique(raw_cdf.instance_path)), (raw_cdf.instance_path,)))
@assert nrow(raw_cdf) == 42
@assert sum(raw_cdf.rotation) == 21
@assert all(!isnan, raw_cdf.build_and_solve_time)
@assert all(==("BUILT_MODEL"), raw_cdf.build_stop_reason)

# %%
# Compute how much of total time was spent on enumeration and purge.
let df = deepcopy(raw_cdf)
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
cdf = let df = deepcopy(raw_cdf)
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
# Plot percentage of barrier time to time
per_barrier_time_c = let df = deepcopy(cdf)
    #filter!(:rotation => iszero, df)
    #@assert all(in(["OPTIMAL", "INFEASIBLE"]), df[!, :stop_reason])
    @assert all(!isnan, df[!, :root_node_time])
    df[!, :percent_barrier] = (df[!, :root_node_time] ./ df[!, :build_and_solve_time]) .* 100
    df
end
plot(per_barrier_time_c, x = "percent_barrier", color="stop_reason", Geom.histogram(;bincount = 10))

# %%
let df = deepcopy(cdf)
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
# for each configuration). To do this we report nonzero amounts, instead of variables and plates.
let df = deepcopy(cdf)
    select!(df,
        :instance_name, :rotation, :mirror_plates, :num_nonzeros, :stop_reason, :build_and_solve_time
    )

    df[!, :num_nonzeros] = number2latex.(df[!, :num_nonzeros]; enclose = false)
    df[!, :build_and_solve_time] = ifelse.(
        df[!, :build_and_solve_time] .> 3600.0,
        ">3600",
        number2latex.(df.build_and_solve_time; enclose = false)
    )
    replace!(df[!, :stop_reason], "OPTIMAL" => "F", "INFEASIBLE" => "N", "TIME_LIMIT" => "X")

    nrdf = filter(:rotation => iszero, df) # no rotation
    rmdf = filter([:rotation, :mirror_plates] => (r, m) -> isone(r) && isone(m), df) # rotation and mirror

    select!.((nrdf, rmdf), ([:instance_name, :build_and_solve_time, :stop_reason, :num_nonzeros],))
    rename!(nrdf,
        :build_and_solve_time => :nr_time, :stop_reason => :nr_status, :num_nonzeros => :nr_nonzeros
    )
    rename!(rmdf,
        :build_and_solve_time => :rm_time, :stop_reason => :rm_status, :num_nonzeros => :rm_nonzeros
    )

    jdf = join(nrdf, rmdf; on = :instance_name)
    sort!(jdf, :instance_name)
    
    final_headers = [
        :instance_name              => "Inst.",
        :nr_status                  => "G. F",
        :nr_nonzeros                => "G. #nz",
        :nr_time                    => "G. T (s)",
        :rm_status                  => "G. R. F",
        :rm_nonzeros                => "G. R. M. #nz",
        :rm_time                    => "G. R. M. (s)"
    ]

    select!(jdf, first.(final_headers)) # This is done to re-order the columns.
    rename!(jdf, final_headers...)

    rename!(jdf, names(jdf) .=> esc_latex.(names(jdf)))

    pretty_table(
        jdf; backend = :latex, nosubheader = true, alignment = vcat([:l], repeat([:r], ncol(jdf) - 1))
    )
end

# %%
raw_df42_path = "./data/G2OPP_clautiaux42.csv"
raw_df42 = DataFrame(CSV.File(raw_df42_path))
#nothing
showtable(raw_df42)

# %%
# These results are a little more messy, let us try to make sense of all these runs.
# No change is done here.
let df = deepcopy(raw_df42)
    sdf = groupby(df, [:rotation, :mirror_plates, :faithful, :round2disc, :pricing]) |> 
        x -> combine(x, :instance_path => length => :nrows)
    show(sdf)

    # For some reason, we did run rotation+mirror two times, let us examine the differences between
    # these runs.
    rm_df = filter([:rotation, :mirror_plates] => (r, m) -> isone(r) && isone(m), df) |>
        x -> select!(x, :instance_path, :build_and_solve_time, :stop_reason, :this_data_file) |>
        x -> sort!(x, :instance_path)
    #show(rm_df)
    # it seems that we just run the same configuration two times (the is no significant changes to the
    # results), keeping the same recent runs
    df = filter!(:this_data_file => x -> !occursin("2021-06-08T18:31", x), df)
    
    # The below is just to confirm the right runs were removed.
    sdf = groupby(df, [:rotation, :mirror_plates, :faithful, :round2disc, :pricing]) |> 
        x -> combine(x, :instance_path => length => :nrows)
    show(sdf)
end

# %%
# Let us confirm some things:
# No faithful, no pricing, always round2disc.
@assert all(iszero, raw_df42.faithful)
@assert all(==("none"), raw_df42.pricing)
@assert all(isone, raw_df42.round2disc)
# Every run with mirror-plates is a rotation one.
@assert all(raw_df42.rotation .>= raw_df42.mirror_plates)
# The purge removed nothing, because no pricing was done.
@assert all(iszero, raw_df42.qt_unreachable_plate_types)
# The times are all present, and all models were built.
@assert all(!isnan, raw_df42.build_and_solve_time)
@assert all(==("BUILT_MODEL"), raw_df42.build_stop_reason)

# %%
# Create df42 (the version of raw_df42) with some pre-treatment
df42 = let df = deepcopy(raw_df42)
    filter!(:this_data_file => x -> !occursin("2021-06-08T18:31", x), df)
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
    
    @show names(df)
    showtable(df)
    df
end

# %%
# Check number of high percentage barrier times.
let df1 = deepcopy(cdf), df2 = deepcopy(df42)
    #filter!(:rotation => iszero, df)
    #@assert all(in(["OPTIMAL", "INFEASIBLE"]), df[!, :stop_reason])
    df = vcat(df1, df2)
    @assert all(!isnan, df[!, :root_node_time])
    df[!, :percent_barrier] = (df[!, :root_node_time] ./ df[!, :build_and_solve_time]) .* 100
    @show sum(df[!, :percent_barrier] .> 12.5) / nrow(df)
    df
end

# %%
let df = deepcopy(df42)
    df_norotation = filter(:rotation => iszero, df)
    df_nomirror = filter([:rotation, :mirror_plates] => ((r, m) -> isone(r) && iszero(m)), df)
    df_mirror = filter([:rotation, :mirror_plates] => ((r, m) -> isone(r) && isone(m)), df)
    @show mean(df_norotation.build_and_solve_time)
    @show mean(df_nomirror.build_and_solve_time)
    @show mean(df_mirror.build_and_solve_time)
    @show mean(df_norotation.num_nonzeros)
    @show mean(df_nomirror.num_nonzeros)
    @show mean(df_mirror.num_nonzeros)
    @show sum(==("OPTIMAL"), df_norotation.stop_reason)
    @show sum(==("OPTIMAL"), df_nomirror.stop_reason)
    @show sum(==("OPTIMAL"), df_mirror.stop_reason)
    qt_nomirror_nonzeros = sum(df_nomirror.num_nonzeros)
    qt_mirror_nonzeros = sum(df_mirror.num_nonzeros)
    percent_nonzeros_kept = qt_mirror_nonzeros / qt_nomirror_nonzeros
    @show percent_nonzeros_kept
    all_percentages = filter([:rotation, :mirror_plates] => ((r, m) -> isone(r) && isone(m)), df).num_nonzeros ./
        filter([:rotation, :mirror_plates] => ((r, m) -> isone(r) && iszero(m)), df).num_nonzeros
    @show minimum(all_percentages)
    @show maximum(all_percentages)
end

# %%
# Plot percentage of barrier time to time
per_barrier_time_df42 = let df = deepcopy(df42)
    #filter!(:rotation => iszero, df)
    @assert all(in(["OPTIMAL", "INFEASIBLE"]), df[!, :stop_reason])
    @assert all(!isnan, df[!, :root_node_time])
    df[!, :percent_barrier] = (df[!, :root_node_time] ./ df[!, :build_and_solve_time]) .* 100
    df
end
plot(per_barrier_time_df42, x = "percent_barrier", color="stop_reason", Geom.histogram(;bincount = 10))

# %%
# Prints a long table with all runs (latex format).
let df = deepcopy(df42)
    select!(df,
        :instance_name, :rotation, :mirror_plates, :qt_pe_vars, :qt_cm_vars,
        :qt_plates, :stop_reason, :build_and_solve_time
    )

    #foreach(println, unique!(sort(df.instance_name)))
    df[!, :build_and_solve_time] = number2latex.(df.build_and_solve_time; enclose = false)
    sort!(df, [:instance_name, :rotation, :mirror_plates])
    replace!(df[!, :stop_reason], "OPTIMAL" => "F", "INFEASIBLE" => "N")
    
    # Add non-guillotine feasibility column from another dataframe.
    ngf_df = DataFrame(CSV.File("./data/G2OPP_clautiaux42_feasibility.csv"))
    df = leftjoin(df, ngf_df; on = :instance_name)
    
    final_headers = [
        :instance_name         => "Inst.",
        :non_guillotine_feasibility => "NG. F.",
        :rotation              => "R.",
        :mirror_plates         => "M.",
        :qt_pe_vars            => "#E. vars",
        :qt_cm_vars            => "#C. vars",
        :qt_plates             => "#Plates",
        :stop_reason           => "G. F.",
        :build_and_solve_time  => "Time (s)"
    ]

    select!(df, first.(final_headers)) # This is done to re-order the columns.
    rename!(df, final_headers...)

    rename!(df, names(df) .=> esc_latex.(names(df)))

    #showtable(df)

    pretty_table(
        df; backend = :latex, nosubheader = true, alignment = vcat([:l], repeat([:r], ncol(df) - 1))
    )
end

# %%
# Prints a medium table with all runs (latex format). The difference from the long table
# is that there is a single row per instance, instead of three (which were a single row
# for each configuration). To do this we report nonzero amounts, instead of variables and plates.
let df = deepcopy(df42)
    select!(df,
        :instance_name, :rotation, :mirror_plates, :num_nonzeros, :stop_reason, :build_and_solve_time
    )

    df[!, :num_nonzeros] = number2latex.(df[!, :num_nonzeros]; enclose = false)
    df[!, :build_and_solve_time] = number2latex.(df.build_and_solve_time; enclose = false)
    replace!(df[!, :stop_reason], "OPTIMAL" => "F", "INFEASIBLE" => "N")

    nrdf = filter(:rotation => iszero, df) # no rotation
    rdf = filter([:rotation, :mirror_plates] => (r, m) -> isone(r) && iszero(m), df) # rotation
    rmdf = filter([:rotation, :mirror_plates] => (r, m) -> isone(r) && isone(m), df) # rotation and mirror

    select!.((nrdf, rdf, rmdf), ([:instance_name, :build_and_solve_time, :stop_reason, :num_nonzeros],))
    #@show names.((nrdf, rdf, rmdf))
    rename!(nrdf,
        :build_and_solve_time => :nr_time, :stop_reason => :nr_status, :num_nonzeros => :nr_nonzeros
    )
    rename!(rdf,
        :build_and_solve_time => :r_time , :stop_reason => :r_status , :num_nonzeros => :r_nonzeros
    )
    rename!(rmdf,
        :build_and_solve_time => :rm_time, :stop_reason => :rm_status, :num_nonzeros => :rm_nonzeros
    )
    
    sort!.((nrdf, rdf, rmdf), (:instance_name,))
    @assert all(rdf.r_status .== rmdf.rm_status)
    
    select!(rmdf, Not(:rm_status))

    jdf = join(nrdf, rdf, rmdf; on = :instance_name)

    # Add non-guillotine feasibility column from another dataframe.
    #ngf_df = DataFrame(CSV.File("./data/G2OPP_clautiaux42_feasibility.csv"))
    #jdf = leftjoin(jdf, ngf_df; on = :instance_name)
    
    final_headers = [
        :instance_name              => "Inst.",
        #:non_guillotine_feasibility => "N. G. F",
        :nr_status                  => "G. F",
        :nr_nonzeros                => "G. #nz",
        :nr_time                    => "G. T (s)",
        :r_status                   => "G. R. F",
        :r_nonzeros                 => "G. R. #nz",
        :r_time                     => "G. R. T (s)",
        :rm_nonzeros                => "G. R. M. #nz",
        :rm_time                    => "G. R. M. (s)"
    ]

    select!(jdf, first.(final_headers)) # This is done to re-order the columns.
    rename!(jdf, final_headers...)

    rename!(jdf, names(jdf) .=> esc_latex.(names(jdf)))

    #showtable(df)

    pretty_table(
        jdf; backend = :latex, nosubheader = true, alignment = vcat([:l], repeat([:r], ncol(jdf) - 1))
    )
end

# %%
let
    df = DataFrame(CSV.File("./data/G2OPP_clautiaux42_feasibility.csv"))
    num_fs = count(==("F"), df.non_guillotine_feasibility)
    num_ns = count(==("N"), df.non_guillotine_feasibility)
    @show num_fs
    @show num_ns
    @assert num_fs + num_ns == 42
end

# %%
let df = deepcopy(df42)
    nrdf = filter(:rotation => iszero, df) # no rotation
    rdf = filter([:rotation, :mirror_plates] => (r, m) -> isone(r) && iszero(m), df) # rotation
    rmdf = filter([:rotation, :mirror_plates] => (r, m) -> isone(r) && isone(m), df) # rotation and mirror

    select!.((nrdf, rdf, rmdf), ([:instance_name, :build_and_solve_time, :stop_reason, :num_nonzeros],))
    rename!(nrdf,
        :build_and_solve_time => :nr_time, :stop_reason => :nr_status, :num_nonzeros => :nr_nonzeros
    )
    rename!(rdf,
        :build_and_solve_time => :r_time , :stop_reason => :r_status , :num_nonzeros => :r_nonzeros
    )
    rename!(rmdf,
        :build_and_solve_time => :rm_time, :stop_reason => :rm_status, :num_nonzeros => :rm_nonzeros
    )

    jdf = join(nrdf, rdf, rmdf; on = :instance_name)
    
    @assert nrow(jdf) == 42
    
    # #=
    nr_faster_than_rm = filter([:nr_time, :rm_time] => (<), jdf)
    @show nrow(nr_faster_than_rm)
    @show nr_faster_than_rm.instance_name
    showtable(nr_faster_than_rm)
    #=
    showtable(
        select!(
            sort!(
                filter(:instance_name => in(nr_faster_than_rm.instance_name), df42),
            [:instance_name, :rotation, :mirror_plates]),
        [:instance_name, :rotation, :mirror_plates, :qt_pe_vars, :qt_cm_vars, :qt_plates])
    )
    =#
    
    nr_less_nz_than_rm = filter([:nr_nonzeros, :rm_nonzeros] => (<), jdf)
    @show nrow(nr_less_nz_than_rm)
    @show nr_less_nz_than_rm.instance_name
    showtable(nr_less_nz_than_rm)
    
    # =#

    # showtable ony works if it is the last thing returned by the cell, so we need to selectively
    # comment these blocks to see the table in an interactive format
    #=
    nr_faster_than_r = filter([:nr_time, :r_time] => (<), jdf)
    @show nr_faster_than_r
    showtable(filter(:instance_name => in(nr_faster_than_r.instance_name), df42))
    =#
end

# %%
?outerjoin

# %%
