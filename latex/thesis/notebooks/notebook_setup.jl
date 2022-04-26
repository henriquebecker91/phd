# Necessary to avoid packages failing because they tried to use the system library (of the wrong version) instead of the artifact that comes with the package.
ENV["LD_LIBRARY_PATH"] = ""
import Pkg
Pkg.activate(".")
Pkg.instantiate()
using WebIO
WebIO.install_jupyter_nbextension()
using DataFrames
using DataFramesMeta
using Gadfly
using IJulia
using Weave
using CSV
using Revise
using TableView
using Printf
using PrettyTables
using Formatting
using StatsBase
import Cairo
import Fontconfig

# ==================== FILLING BASE HOLES ====================

function Base.parse(::Type{Rational{T}}, s::AbstractString) where {T}
    ns, ds = split(s, '/'; keepempty = false)
    n = parse(T, ns)
    d = parse(T, ds)
    return n//d
end

# ==================== GENERAL UTILITY FUNCTIONS ====================

broadwrap(f) = function (args...) broadcast(f, args...) end
esc_latex(s) = replace(s, r"(#|%|&|_)" => s"\\\1")
function number2latex(num; enclose = true)
    if ismissing(num)
        "--"
    elseif isa(num, Integer)
        (s -> enclose ? "\\($s\\)" : s)(sprintf1("%'d", num))
    elseif isa(num, AbstractFloat)
        (s -> enclose ? "\\($s\\)" : s)(sprintf1("%'.2f", num))
    elseif isa(num, NTuple{2, Number})
        "$(number2latex(num[1]; enclosed = enclose)) ($(number2latex(num[2]; enclosed = enclose)))"
    else
        error("Unexpected type.")
    end
end

# ==================== METHODS FOR highlight_best_values! ====================
wrap_in_textbf(str) = "\\textbf{$str}"
wrap_in_textit(str) = "\\textit{$str}"
# See https://www.pcre.org/current/doc/html/pcre2syntax.html#SEC2
const MATH_MODE_REGEX = Regex(
    "^(?:\\Q\$\\E|\\Q\\(\\E)(.+)(?:\\Q\$\\E|\\Q\\)\\E)\$"
)
function dirt2number(
  dirt, decimal_separator = '.', failure = NaN
) :: Float64
		gold = filter(c -> isdigit(c) | (c == decimal_separator), dirt)
    number = tryparse(Float64, gold)
    isnothing(number) && return failure
    return number
end
function latex2number(
  str, failure = NaN; separator = nothing
) :: Float64
    str = isnothing(separator) ? str : replace(str, separator => "")
    v = tryparse(Float64, str)
    if isnothing(v)
        m = match(MATH_MODE_REGEX, str)
        if !isnothing(m)
          v = tryparse(Float64, m.captures[1])
        end
    end
    return isnothing(v) ? failure : v
end
find_index_min(a) = isempty(a) ? Int[] : Int[findmin(a)[2]]
find_index_max(a) = isempty(a) ? Int[] : Int[findmax(a)[2]]
find_index_allmin(a) = isempty(a) ? Int[] : findall(isequal(findmin(a)[1]), a)
find_index_allmax(a) = isempty(a) ? Int[] : findall(isequal(findmax(a)[1]), a)

function highlight_best_values!(
    df;
    # Set of columns to be considered.
    columns = 1:ncol(df),
    # Transforms the cell body into another object.
    cleaner = latex2number,
    # Return true if the object returned by cleaner should be ignored.
    ignorer = isnan,
    # Takes a row of clean non-ignored values and returns the
    # indexes of the ones that should be highlighted.
    chooser = find_index_min,
    # Takes the cell body and return the highlighted cell body.
    changer = wrap_in_textbf
)
    for row_id in 1:nrow(df)
        row = df[row_id, :]
        selected_values = []
        selected_columns = empty(columns)
        for col_id in columns
            value = cleaner(row[col_id])
            if !ignorer(value)
                push!(selected_values, value)
                push!(selected_columns, col_id)
            end
        end
        chosen = chooser(selected_values)
        for col_id in selected_columns[chosen]
            df[row_id, col_id] = changer(df[row_id, col_id])
        end
    end
    return
end

