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
        args2id = Dict{NTuple{3, Bool}, String}(
            (false, true, true) => "complete",
            (false, false, true) => "only_redundant_cut",
            (false, true, false) => "only_cut_position",
            (false, false, false) => "both_reductions",
            (true, false, false) => "priced"
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
# Shows info about the unfinished runs. Theoretically those three conditions should always
# appear together in a row, of some row do not have all of them, then something is wrong.
@linq fr_df |>
    where(.!:finished .| :had_timeout .| (:build_stop_reason .== "NOT_REACHED")) |>
    select(:instance_name, :model_variant, :datafile, :finished, :had_timeout, :build_stop_reason)

# %%
# Shows info about OPTIMAL_FOUND runs
@linq fr_df |>
    where(:build_stop_reason .== "FOUND_OPTIMUM") |>
    select(:instance_name, :model_variant, :datafile, :finished, :build_stop_reason) |>
    (tmp_df -> for s in sort(tmp_df[!, :instance_name]); println(s); end)()


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
# Defines the 'origin' column. 
th59 = let th59 = th59
    th59[!, :origin] .= "thomopulosthesis"
    th59
end
unique(th59[!, :model_variant])

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
# Absolutize the vars and plates columns, i.e., now thay all refer to absolute
# quantities of variables and plates, instead of having all rows with
# 'model_variant != complete' with percentages (as it is in the CSV file).
th59 = let th59 = deepcopy(th59)
    # Fix the missing plates. The number of plates is not described in the thesis but
    # a safe interpretation is that the number stayed the same from the previous step.
    if any(ismissing, th59[!, :plates])
        reduced_plates = Dict{String, Float64}()
        restricted_plates = Dict{String, Float64}()

        @byrow! th59 begin
            if :model_variant == "both_reductions"
                reduced_plates[:instance_name] = :plates
            elseif :model_variant == "restricted"
                restricted_plates[:instance_name] = :plates
            end
        end
        th59 = @byrow! th59 begin
            if :model_variant == "priced"
                :plates = reduced_plates[:instance_name]
            elseif :model_variant == "after_restricted_pricing"
                :plates = restricted_plates[:instance_name]
            end
        end
    end
    show(@where(th59, ismissing.(:plates)))
    th59[!, :plates] = Float64.(th59[!, :plates])
    if eltype(th59[!, :vars]) != Int
        th59_complete = @where(th59, :model_variant .== "complete")
        stats_type = NamedTuple{(:vars, :plates),Tuple{Int64,Int64}}
        th59_complete_dict = Dict{String, stats_type}()
        @byrow! th59_complete begin
            th59_complete_dict[:instance_name] = (vars = :vars, plates = :plates)
        end
        for column in [:vars, :plates]
            complete_stats = getindex.(getindex.((th59_complete_dict,), th59.instance_name), column)
            th59[!, column] = ifelse.(
                th59.model_variant .== "complete",
                complete_stats,
                round.(Int64, complete_stats .* th59[!, column] ./ 100.0)
            )
        end
    end
    th59
end
showtable(th59)

# %%
@show unique(th59[!, :model_variant])
@show names(fr_df)[.!isnothing.(match.(r"qt_cmvars.*", names(fr_df)))]
nothing

# %%
# Put the faithful reimplementation table in the same format as the data extracted from
# Thomopulos 2016 (initially, only include the rows with model_variant in "complete", "only_cut_position",
# "only_redundant_cut" and "both_reductions").
short_fr_df = let
    sel_cols = [
        :instance_name, :model_variant, :qt_pevars_after_preprocess, :qt_cmvars_after_preprocess, 
        :qt_plates_after_preprocess
    ]
    short_fr_df = select(fr_df, sel_cols)
    short_fr_df[!, :origin] .= "reimplementation"
    variants = ("complete", "only_cut_position", "only_redundant_cut", "both_reductions")
    short_fr_df = @where(short_fr_df, in.(:model_variant, (variants,)))
    # basically the 'x' and 'y' of the model
    var_categories = [:qt_pevars_after_preprocess, :qt_cmvars_after_preprocess]
    select!(short_fr_df, var_categories => ((a,b) -> a .+ b) => :vars, Not(var_categories))
    rename!(short_fr_df, :qt_plates_after_preprocess => :plates)
    
    short_fr_df
end

# %%
# Add thee rows with model_variant equal to 'priced', 'restricted', 'after_restricted_pricing',
# and purged to 'short_fr_df', those come from the 'priced' runs (i.e., the 'priced' runs execute
# the restricted pricing as a phase of their process, so we do not need runs just to
# get such values, they are just in another column of 'priced' rows).
short_fr_df = if "restricted" in short_fr_df[!, :model_variant]
    short_fr_df
else let short_fr_df = deepcopy(short_fr_df)
    needed_cols = [
        :instance_name, :model_variant,
        # restricted
        :qt_pevars_restricted, :qt_cmvars_restricted, :qt_plates_restricted,
        # priced restricted
        :qt_pevars_priced_restricted, :qt_cmvars_priced_restricted, :qt_plates_priced_restricted,
        # after purge
        :qt_pevars_after_purge, :qt_cmvars_after_purge, :qt_plates_after_purge,
        # after pricing
        :qt_pevars_after_pricing, :qt_cmvars_after_pricing, :qt_plates_after_pricing
    ]
    # Get just the needed columns and rows.
    cut = @linq fr_df |> select(needed_cols) |> where(:model_variant .== "priced")
    # Create the origin column.
    cut[!, :origin] .= "reimplementation"

    function to_th59_format!(df, col_suffix, final_cols, variant_name)
        plate_colname = Symbol("qt_plates_" * string(col_suffix))
        var_colnames = Symbol[
            Symbol("qt_pevars_" * string(col_suffix)),
            Symbol("qt_cmvars_" * string(col_suffix))
        ]
        rename!(df, plate_colname => :plates)
        transform!(df, var_colnames => ((a,b) -> a .+ b) => :vars)
        df[!, :model_variant] .= variant_name
        select!(df, final_cols)
        df
    end

    final_columns = [:instance_name, :model_variant, :vars, :plates, :origin]
    restricted = to_th59_format!(deepcopy(cut), "restricted", final_columns, "restricted")
    priced = to_th59_format!(deepcopy(cut), "after_pricing", final_columns, "priced")
    priced_r = to_th59_format!(deepcopy(cut), "priced_restricted", final_columns, "after_restricted_pricing")
    purged = to_th59_format!(deepcopy(cut), "after_purge", final_columns, "purged")
    
    vcat(short_fr_df, restricted, priced_r, priced, purged; cols = :setequal)
end #= else =# end #= let =#

# %%
# Create a mock "purged" rowset for th59 that is a copy of the
# "priced" to allow comparison with the extra data generated in
# our experiments.
th59 = let
    mock_purged = @where(th59, :model_variant .== "priced")
    mock_purged[!, :model_variant] .= "purged"
    vcat(th59, mock_purged; cols = :setequal)
end

# %%
# Just show how much data is missing for each variant.
missing_from_fr = combine(
    #groupby(th59, :model_variant),
    groupby(short_fr_df, :model_variant),
    :model_variant => length => :num_rows,
    :vars => (sum ∘ broadwrap(ismissing)) => :missing_vars,
    :plates => (sum ∘ broadwrap(ismissing) => :missing_plates)
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
full_df = vcat(th59, short_fr_df; cols = :setequal)
nothing
#showtable(fulldf)

# %%
showtable(
    #@linq @where(full_df, (:model_variant .== "priced")) |>
    @linq @where(full_df, (:instance_name .== "gcut1")) |>
        sort!([:instance_name, :origin])
)

# %%
@linq full_df |>
    where((:instance_name .== "gcut2") .& (:model_variant .== "restricted"))

# %%
# Transforms the long format into a wide one and rename many columns and
# string cells to the final name they will have in the table.
function make_header(origin, variable)
    @assert origin ∈ ["reimplementation", "thomopulosthesis"]
    @assert variable ∈ ["vars", "plates"]
    prefix = origin == "reimplementation" ? "R. %" : "O. #"
    suffix = variable == "vars" ? "v" : "p" 
    return prefix * suffix
end
make_headers(origins, variables) = make_header.(origins, variables)

#transform!(full_df, :origin => broadwrap(Symbol) => :origin)
#transform!(full_df, :model_variant => broadwrap(Symbol) => :model_variant)

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
                "restricted" => (longname = "Restricted PP-G2KP", order = 5),
        "after_restricted_pricing" => (longname = "Priced Restricted PP-G2KP", order = 6),
        "priced" => (longname = "(no purge) Priced PP-G2KP", order = 7),
        "purged" => (longname = "Priced PP-G2KP", order = 8),
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

