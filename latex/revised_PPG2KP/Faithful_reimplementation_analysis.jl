# -*- coding: utf-8 -*-
# %%
# Ensure your working directory is: https://github.com/henriquebecker91/phd/tree/master/latex/revised_PPG2KP
include("notebook_setup.jl")

# %%
# Read the data, show nothing.
fr_csv_path = "./data/faithful_reimplementation.csv"
fr_df = DataFrame(CSV.File(fr_csv_path))
#showtable(lp_method_df)
nothing

# %%
# Clean the data a little. Try to keep this safe to re-apply to a cleaned dataframe, if possible.
fr_df = let
    # Keep only the instance name (not the path).
    @with(fr_df, :instance_name .= basename.(:instance_name))
    @with(fr_df, :datafile .= basename.(:datafile))
    colnames = names(fr_df)
    if any(name -> name ∈ colnames, ["pricing_method", "disabled_redundant_cut", "disabled_cut_position"])
        args2id = Dict{NTuple{3, Bool}, Symbol}(
            (false, true, true) => :complete,
            (false, false, true) => :only_redundant_cut,
            (false, true, false) => :only_cut_position,
            (false, false, false) => :both_reductions,
            (true, false, false) => :priced
        )
        fr_df.params = getindex.((args2id,), tuple.(
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
#showtable(fr_df)
nm_fr_df = dropmissing(fr_df)

# %%
@linq fr_df |>
    where((:finished .== false) .| isnan.(:total_instance_time)) |>
    select(:instance_name, :params, :datafile, :total_instance_time, :finished)

# %%
@linq fr_df |>
    groupby(:params) |>
    based_on(; qt = length(:params))

# %%
# The data extracted from DOI 10.6092/unibo/amsdottorato/7399 by hand. The amount of variables in many
# model variants for the 59 instances used by Thomopulos in their 2016 thesis.
# The columns 'vars' and 'plates' have absolute values for rows in which 'model_variant' column is
# 'complete', and percentages for the rest, as in the thesis tables.
th59_csv_path = "./data/thomopulos_2016.csv"
th59 = DataFrame(CSV.File(th59_csv_path))
#showtable(th59)
nothing

# %%
# Just some checking if did not forgot anything.
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
let
    already_extracted = unique!(sort(th59.instance_name))
    display(setdiff(THOMOPULOS_THESIS_INSTANCES, already_extracted))
    nothing
end

# %%
# Changes 'instance_name' and 'model_variant' to symbols (more semantically appropriate).
th59 = let
    th59.instance_name = Symbol.(th59.instance_name)
    th59.model_variant = Symbol.(th59.model_variant)
    th59
end

# %%
# Absolutize the vars and plates columns, i.e., now thay all refer to absolute
# quantities of variables and plates, instead of having all rows with
# 'model_variant != complete' with percentages (as it is in the CSV file).
let
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
            round.(Union{Int64, Missing}, complete_stats .* th59[!, column])
        )
    end
    th59
end

# %%
short_fr_df = deepcopy(th59)
short_fr_df[!, :origin] .= :reimplementation
short_fr_df.vars = round.(Int64, short_fr_df.vars .* 0.95)
th59[!, :origin] .= :thomopulosthesis
full_df = vcat(short_fr_df, th59)
function make_header(origin, variable)
    @assert origin ∈ [:reimplementation, :thomopulosthesis]
    @assert variable ∈ ["vars", "plates"]
    prefix = origin == :reimplementation ? "R. %" : "O. #"
    suffix = variable == "vars" ? "v" : "p" 
    return prefix * suffix
end
make_headers(origins, variables) = make_header.(origins, variables)

fr_table = @linq full_df |>
    select(:model_variant, :origin, :vars, :plates) |>
    groupby([:model_variant, :origin]) |>
    based_on(vars = sum(:vars), plates = sum(:plates)) |>
    stack([:vars, :plates]) |>
    unstack([:origin, :variable], :model_variant, :value) |>
    select!([:origin, :variable] => make_headers => :header, Not([:origin, :variable])) |>
    stack(Not(:header)) |>
    unstack(:variable, :header, :value)

fr_table[!, "R. %p"] = fr_table[!, "R. %p"] ./ fr_table[!, "O. #p"]
fr_table[!, "R. %v"] = fr_table[!, "R. %v"] ./ fr_table[!, "O. #v"]
#fr_table[:, ["variable", "O. #v", "R. %v", "O. #p", "R. %p"]]
select!(fr_table, ["variable", "O. #v", "R. %v", "O. #p", "R. %p"])
select!(fr_table, :variable => "Variant", Not(:variable))
pretty_variant_names = Dict{String, Any}(
    "complete" => (longname = "Complete PP-G2KP", order = 1),
    "only_cut_position" => (longname = "Complete +Cut-Position", order = 2),
    "only_redundant_cut" => (longname = "Complete +Redundant-Cut", order = 3),
    "both_reductions" => (longname = "PP-G2KP (CP + RC)", order = 4),
    "priced" => (longname = "Priced PP-G2KP", order = 5)
)
sort!(fr_table, "Variant"; by = (name -> pretty_variant_names[name].order))
fr_table[!, "Variant"] = getindex.(getindex.((pretty_variant_names,), fr_table[!, "Variant"]), :longname)
fr_table

# %%
?unstack
