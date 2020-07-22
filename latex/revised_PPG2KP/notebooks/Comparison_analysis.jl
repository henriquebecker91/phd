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
raw_cdf_path = "./data/comparison.csv"
raw_cdf = DataFrame(CSV.File(raw_cdf_path))
#nothing
showtable(raw_cdf)

# %%
# Clean the data a little. Try to keep this safe to re-apply to a cleaned dataframe, if possible.
cdf = let cdf = deepcopy(raw_cdf)
    # Keep only the instance name (not the path).
    @with(cdf, :instance_name .= basename.(:instance_name))
    @with(cdf, :datafile .= basename.(:datafile))
    colnames = names(cdf)
    for name in colnames
        column = cdf[!, name]
        if eltype(column) <: AbstractFloat
            cdf[!, name] = ifelse.(isnan.(column), missing, column)
        end
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
# Check how many not finished in each group.
@linq cdf |> groupby(:model_variant) |>
    based_on(; qt_not_finished = sum(.!:finished))

# %%
# Check which instances did not finish.
@linq cdf |> where(.!:finished) |>
    select(:instance_name, :model_variant, :total_pricing_time, :datafile) |>
    sort([:instance_name, :model_variant])

# %%
# Create the comparison time table
ctt = let cdf = deepcopy(cdf)
    used_columns = [
        :model_variant, :restricted_pricing_time, :iterated_pricing_time, :final_pricing_time,
        :final_solving_time, :total_instance_time, :finished
    ]
    select!(cdf, used_columns)
    filter!(:finished => identity, cdf)
    select!(cdf, Not(:finished))
    for column in names(cdf)
        if eltype(cdf[!, column]) == Union{Missing, Float64}
            col_len = length(cdf[!, column])
            cdf[!, column] = ifelse.(ismissing.(cdf[!, column]), zeros(col_len), cdf[!, column])
        end
    end
    # length(unique!(sort(raw_cdf[!, "instance_name"])))
    cdf = @linq cdf |> groupby(:model_variant) |> based_on(;
        qt_solved = length(:model_variant),
        restricted_pricing_time = sum(:restricted_pricing_time),
        iterated_pricing_time = sum(:iterated_pricing_time),
        final_pricing_time = sum(:final_pricing_time),
        final_solving_time = sum(:final_solving_time),
        total_instance_time = sum(:total_instance_time)
    )
    for percent_col in [:restricted_pricing_time, :iterated_pricing_time, :final_pricing_time, :final_solving_time]
        cdf[!, percent_col] ./= cdf[!, :total_instance_time]
        cdf[!, percent_col] .*= 100.0
    end
    cdf = @transform(cdf, avg_instance_time .= :total_instance_time ./ :qt_solved)
    cdf[!, :total_instance_time] = trunc.(Int, cdf[!, :total_instance_time])
    qt_instances = 59
    time_limit_secs = 3 * 60 * 60 # three hours
    cdf[!, :total_with_timeouts] = cdf[!, :total_instance_time] .+
        (qt_instances .- cdf[!, :qt_solved]) .* time_limit_secs
    
    select!(cdf,
        :model_variant => "Variant",
        :total_with_timeouts => "T. (s)",
        :qt_solved => "#",
        :avg_instance_time => "S. A. (s)",
        :total_instance_time => "S. T. (s)",
        :restricted_pricing_time => "RP (%)",
        :iterated_pricing_time => "IP (%)",
        :final_pricing_time => "FP (%)",
        :final_solving_time => "FS (%)"
    )

    @show names(cdf)
    #@show unique(cdf[!, "Variant"])
    pretty_variant_names = Dict{String, NamedTuple{(:pretty_name,:order), Tuple{String, Int}}}(
        "faithful" => (pretty_name = "Original", order = 1),
        "simple_revised" => (pretty_name = "Enhanced", order = 2),
        "rounded_faithful" => (pretty_name = "O. +Rounding", order = 3),
        "rounded_revised" => (pretty_name = "E. +Rounding", order = 4),
        "warmed_rounded_faithful" => (pretty_name = "O. +R. +Warming", order = 5),
        "warmed_rounded_revised" => (pretty_name = "E. +R. +Warming", order = 6),
        "priced_faithful" => (pretty_name = "P. Faithful +R. +W.", order = 7),
        "priced_revised" => (pretty_name = "Priced E. +R. +W.", order = 8),
        "no_purge_priced_revised" => (pretty_name = "P. E. +R. +W. -Purge", order = 9)
    )
    sort!(cdf, "Variant"; by = (name -> pretty_variant_names[name].order))
    cdf[!, "Variant"] = getindex.(getindex.((pretty_variant_names,), cdf[!, "Variant"]), :pretty_name)
    
    function num2latex(num)
        if ismissing(num)
            "--"
        elseif isa(num, Integer)
            "\\(" * sprintf1("%'d", num) * "\\)" 
        elseif isa(num, AbstractFloat)
            "\\(" * sprintf1("%'.2f", num) * "\\)"
        else
            error("Unexpected type.")
        end
    end
    numeric_columns = [
        "#", "T. (s)", "S. A. (s)", "S. T. (s)", "RP (%)", "IP (%)", "FP (%)", "FS (%)"
    ]
    for col in numeric_columns
        cdf[!, col] = num2latex.(cdf[!, col])
    end
    select!(cdf, names(cdf) .=> esc_latex.(names(cdf)))
    cdf
end

pretty_table(
    ctt; backend = :latex, nosubheader = true, alignment = [:l, :r, :r, :r, :r, :r, :r, :r, :r]
)

showtable(ctt)

# %%
?filter!
