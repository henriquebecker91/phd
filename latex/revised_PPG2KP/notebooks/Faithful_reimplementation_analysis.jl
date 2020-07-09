# -*- coding: utf-8 -*-
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
fr_csv_path = "./data/faithful_reimplementation.csv"
fr_df = DataFrame(CSV.File(fr_csv_path))
nothing
#showtable(fr_df)

# %%
# Clean the data a little. Try to keep this safe to re-apply to a cleaned dataframe, if possible.
fr_df = let
    # Keep only the instance name (not the path).
    @with(fr_df, :instance_name .= basename.(:instance_name))
    @with(fr_df, :datafile .= basename.(:datafile))
    colnames = names(fr_df)
    for name in colnames
        column = fr_df[!, name]
        if eltype(column) <: AbstractFloat
            fr_df[!, name] = ifelse.(isnan.(column), missing, column)
        elseif eltype(column) <: Integer && occursin(r"^qt_", string(name))
            fr_df[!, name] = ifelse.(column .< 0, missing, column)
        end
    end
    if any(name -> name ∈ colnames, ["pricing_method", "disabled_redundant_cut", "disabled_cut_position"])
        args2id = Dict{NTuple{3, Bool}, Symbol}(
            (false, true, true) => :complete,
            (false, false, true) => :only_redundant_cut,
            (false, true, false) => :only_cut_position,
            (false, false, false) => :both_reductions,
            (true, false, false) => :priced
        )
        fr_df.model_variant = getindex.((args2id,), tuple.(
            fr_df.pricing_method .== "furini", fr_df.disabled_redundant_cut, fr_df.disabled_cut_position
        ))
        select!(fr_df, Not([:pricing_method, :disabled_redundant_cut, :disabled_cut_position]))
    end
    
    fr_df
end
nothing

# %%
# Shows the cleaned data.
println(names(fr_df))
showtable(fr_df)

# %%
@linq fr_df |>
    where((:instance_name .== "gcut2") .& (:model_variant .== ^(:priced))) |>
    select(:instance_name, :model_variant, :finished, :total_instan)

# %%
# Shows info about the unfinished runs.
@linq fr_df |>
    where((:finished .== false) .| isnan.(:total_instance_time)) |>
    select(:instance_name, :model_variant, :datafile, :total_instance_time, :finished)

# %%
# Check if the number of lines per model variant is correct.
@linq fr_df |>
    groupby(:model_variant) |>
    based_on(; qt = length(:model_variant))

# %%
# The data extracted from DOI 10.6092/unibo/amsdottorato/7399 by hand. The amount of variables in many
# model variants for the 59 instances used by Thomopulos in their 2016 thesis.
# The columns 'vars' and 'plates' have absolute values for rows in which 'model_variant' column is
# 'complete', and percentages for the rest, as in the thesis tables.
th59_csv_path = "./data/thomopulos_2016.csv"
th59 = DataFrame(CSV.File(th59_csv_path))
showtable(th59)
#nothing

# %%
const THOMOPULOS_THESIS_INSTANCES = vcat(
  # unweighted
  "gcut" .* string.(1:12)
  , String.(split("wang20 2s 3s A1s A2s STS2s STS4s"))
  , ["OF1", "OF2", "W", "CHL1s", "CHL2s"]
  , "A" .* string.(3:5)
  , "CHL" .* string.(5:7)
  , ["CU1", "CU2"]
  , "Hchl" .* split("3s 4s 6s 7s 8s")
  # weighted
  , "cgcut" .* string.(1:3)
  , "okp" .* string.(1:5)
  , String.(split("HH 2 3 A1 A2 STS2 STS4 CHL1 CHL2"))
  , "CW" .* string.(1:3)
  , "Hchl" .* ["2", "9"]
) :: Vector{String}
nothing

# %%
# Changes 'instance_name' and 'model_variant' to symbols (more semantically appropriate).
# Defines the 'origin' column. 
th59 = let
    th59[!, :origin] .= :thomopulosthesis
    th59.instance_name = Symbol.(th59.instance_name)
    th59.model_variant = Symbol.(th59.model_variant)
    th59
end

# %%
# Absolutize the vars and plates columns, i.e., now thay all refer to absolute
# quantities of variables and plates, instead of having all rows with
# 'model_variant != complete' with percentages (as it is in the CSV file).
th59 = let th59 = deepcopy(th59)
    if eltype(th59[!, :vars]) != Int
        th59_complete = @where(th59, :model_variant .== ^(:complete))
        stats_type = NamedTuple{(:vars, :plates),Tuple{Int64,Int64}}
        th59_complete_dict = Dict{Symbol, stats_type}()
        @byrow! th59_complete begin
            th59_complete_dict[:instance_name] = (vars = :vars, plates = :plates)
        end
        for column in [:vars, :plates]
            complete_stats = getindex.(getindex.((th59_complete_dict,), th59.instance_name), column)
            th59[!, column] = ifelse.(
                th59.model_variant .== :complete,
                complete_stats,
                round.(Union{Int64, Missing}, complete_stats .* th59[!, column] ./ 100.0)
            )
        end
    end
    th59
end

# %%
# Put the faithful reimplementation table in the same format as the data extracted from
# Thomopulos 2016.
short_fr_df = let
    sel_cols = [
        :instance_name, :model_variant, :qt_piece_types, :qt_cmvars_pre_purge,
        :qt_plates_pre_purge
    ]
    short_fr_df = select(fr_df, sel_cols)
    short_fr_df[!, :instance_name] = Symbol.(short_fr_df[!, :instance_name])
    short_fr_df[!, :origin] .= :reimplementation
    var_categories = [:qt_piece_types, :qt_cmvars_pre_purge] # basically the 'x' and 'y' of the model
    select!(short_fr_df, var_categories => ((a,b) -> a .+ b) => :vars, Not(var_categories))
    rename!(short_fr_df, :qt_plates_pre_purge => :plates)
    
    short_fr_df
end

# %%
# Create and add rows with model_variant equal to 'restricted' or 'after_restricted_pricing'
# to 'short_fr_df', those come from the 'priced' runs (i.e., the 'priced' runs execute
# the restricted pricing as a phase of their process, so we do not need runs just to
# get such values, they are just in another column of 'priced' rows).
short_fr_df = ("restricted" in names(short_fr_df)) ? short_fr_df : (let short_fr_df = deepcopy(short_fr_df)
    needed_cols = [
        :instance_name, :model_variant, :qt_piece_types, :qt_cmvars_restricted,
        :qt_plates_restricted, :qt_cmvars_before_iterated
    ]
    # Get just the needed columns and rows.
    cut = @linq fr_df |> select(needed_cols) |> where(:model_variant .== ^(:priced))
    # Some changes to normalize with short_fr_df.
    cut[!, :instance_name] = Symbol.(cut[!, :instance_name])
    cut[!, :origin] .= :reimplementation
    # We will use 'cut' for 'restricted' and ''
    cut_copy = deepcopy(cut)
    
    final_columns = [:instance_name, :model_variant, :vars, :plates, :origin]
    
    # Transform 'cut' into the 'restricted' set of rows.
    restricted_vars = [:qt_piece_types, :qt_cmvars_restricted] # basically the 'x' and 'y' of the model
    rename!(cut, :qt_plates_restricted => :plates)
    transform!(cut, restricted_vars => ((a,b) -> a .+ b) => :vars)
    select!(cut, final_columns)
    cut[!, :model_variant] .= :restricted

    # Transform 'cut_copy' into the 'after_restricted_pricing' set of rows.
    restricted_vars = [:qt_piece_types, :qt_cmvars_before_iterated]
    rename!(cut_copy, :qt_plates_restricted => :plates)
    transform!(cut_copy, restricted_vars => ((a,b) -> a .+ b) => :vars)
    select!(cut_copy, final_columns)
    cut_copy[!, :model_variant] .= :after_restricted_pricing
    
    vcat(short_fr_df, cut, cut_copy; cols = :setequal)
end)

# %%
# Just show how much data is missing for each variant.
combine(
    groupby(th59, :model_variant),
    #groupby(short_fr_df, :model_variant),
    :model_variant => length,
    :vars => (sum ∘ broadwrap(ismissing)),
    :plates => (sum ∘ broadwrap(ismissing))
)

# %%
# Removes from both dataframes ('th59' and 'short_fr_df') all rows that have 'missing'
# in the 'vars' column of 'short_fr_df' 
let
    expected_cols = sort!(["instance_name", "model_variant", "origin", "vars", "plates"])
    @assert sort!(names(th59)) == expected_cols
    @assert sort!(names(short_fr_df)) == expected_cols
    table_keys = ["instance_name", "model_variant", "origin"]
    sort!(th59, table_keys)
    sort!(short_fr_df, table_keys)
    to_delete = ismissing.(short_fr_df[!, "vars"])
    delete!(th59, to_delete)
    delete!(short_fr_df, to_delete)
    @assert th59[!, "instance_name"] == short_fr_df[!, "instance_name"]
    @assert th59[!, "model_variant"] == short_fr_df[!, "model_variant"]
    nothing
end

# %%
# Create a dataframe with just the purged data.
purged_df = let 
    sel_cols = [
        :instance_name, :model_variant, :qt_piece_types, :qt_cmvars_pos_purge,
        :qt_plates_pre_purge
    ]
    purged_df = select(fr_df, sel_cols)
    filter!(:model_variant => (mv -> mv == :priced), purged_df)
    purged_df[!, :model_variant] .= :purged
    purged_df[!, :instance_name] = Symbol.(purged_df[!, :instance_name])
    purged_df[!, :origin] .= :reimplementation
    var_categories = [:qt_piece_types, :qt_cmvars_pos_purge] # basically the 'x' and 'y' of the model
    select!(purged_df, var_categories => ((a,b) -> a .+ b) => :vars, Not(var_categories))
    rename!(purged_df, :qt_plates_pre_purge => :plates)
    
    purged_df
end

# %%
# Concatenate the purged rows in the dataframe.
#short_fr_df = vcat(short_fr_df, purged_df; cols = :setequal)

# %%
# Band-aid while the data of the faithful reimplementation has missing values.
# NOT TO BE USED IN THE PAPER.
short_fr_df = let short_fr_df = deepcopy(short_fr_df)
    # Create the :restricted and :after_restricted_pricing rows while the new round of experiments
    # is not done.
    mock_data = @where(th59, in.(:model_variant, ([^(:restricted), ^(:after_restricted_pricing)],)))
    mock_data[!, :origin] .= :reimplementation
    short_fr_df = vcat(short_fr_df, mock_data; cols = :setequal)
    # Solve the problem of the ones missing by timeout.
    sort_keys = [:model_variant, :instance_name]
    sort!(short_fr_df, sort_keys)
    sorted_th59 = sort(th59, sort_keys)
    @assert short_fr_df[!, :instance_name] == sorted_th59[!, :instance_name]
    @assert short_fr_df[!, :model_variant] == sorted_th59[!, :model_variant]
    incomplete_vars = short_fr_df[!, :vars]
    short_fr_df[!, :vars] = ifelse.(ismissing.(incomplete_vars), th59[!, :vars], incomplete_vars)

    short_fr_df
end

# %%
full_df = vcat(th59, short_fr_df; cols = :setequal)
nothing
#showtable(fulldf)

# %%
@linq full_df |>
    where((:instance_name .== ^(:gcut2)) .& (:model_variant .== ^(:restricted)))

# %%
# Transforms the long format into a wide one and rename many columns and
# string cells to the final name they will have in the table.
function make_header(origin, variable)
    @assert origin ∈ [:reimplementation, :thomopulosthesis]
    @assert variable ∈ ["vars", "plates"]
    prefix = origin == :reimplementation ? "R. %" : "O. #"
    suffix = variable == "vars" ? "v" : "p" 
    return prefix * suffix
end
make_headers(origins, variables) = make_header.(origins, variables)

fr_table = let
    fr_table = @linq full_df |>
        select(:model_variant, :origin, :vars, :plates) |>
        groupby([:model_variant, :origin]) |>
        based_on(vars = sum(:vars), plates = sum(:plates)) |>
        stack([:vars, :plates]) |> # vars and plates become values in column 'variable', values in 'value'
        unstack([:origin, :variable], :model_variant, :value) |>
        select!([:origin, :variable] => make_headers => :header, Not([:origin, :variable])) |>
        stack(Not(:header)) |>
        unstack(:variable, :header, :value)
    
    fr_table[!, "R. %p"] = round.(100.0 .* fr_table[!, "R. %p"] ./ fr_table[!, "O. #p"]; digits = 2)
    fr_table[!, "R. %v"] = round.(100.0 .* fr_table[!, "R. %v"] ./ fr_table[!, "O. #v"]; digits = 2)
    #fr_table[:, ["variable", "O. #v", "R. %v", "O. #p", "R. %p"]]
    select!(fr_table, ["variable", "O. #v", "R. %v", "O. #p", "R. %p"])
    select!(fr_table, :variable => "Variant", Not(:variable))
    pretty_variant_names = Dict{String, Any}(
        "complete" => (longname = "Complete PP-G2KP", order = 1),
        "only_cut_position" => (longname = "Complete +Cut-Position", order = 2),
        "only_redundant_cut" => (longname = "Complete +Redundant-Cut", order = 3),
        "both_reductions" => (longname = "PP-G2KP (CP + RC)", order = 4),
        "priced" => (longname = "Priced PP-G2KP", order = 5),
        "restricted" => (longname = "Restricted PP-G2KP", order = 6),
        "after_restricted_pricing" => (longname = "After Restricted Pricing", order = 7)
    )
    sort!(fr_table, "Variant"; by = (name -> pretty_variant_names[name].order))
    fr_table[!, "Variant"] = getindex.(getindex.((pretty_variant_names,), fr_table[!, "Variant"]), :longname)

    fr_table
end

# %%
# Finally, the dataframe is transformed into a latex table.
let
    latex_df = copy(fr_table)
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
    for col in ["O. #v", "R. %v", "O. #p", "R. %p"]
        latex_df[!, col] = num2latex.(latex_df[!, col])
    end
    # Needs to escape some symbols for latex.
    rename!(latex_df,
        ["Variant", "O. \\#v", "R. \\%v",  "O. \\#p", "R. \\%p"]
    )
    #=rename!(
        "Variant", "O. #v" => "O. \\#v", "R. %v" => "R. \\%v", "O. #p" => "O. \\#p", 
        "R. %p" => "R. \\%p"
    )=#
    #=
    for int_col in ["O. #v", "O. #p"]
        latex_df[!, int_col] = 
    end
    for float_col in ["R. %v", "R. %p"]
    end=#
    pretty_table(
        latex_df; backend = :latex, nosubheader = true, alignment = [:l, :r, :r, :r, :r]
    )
end

# %%
?pretty_table

