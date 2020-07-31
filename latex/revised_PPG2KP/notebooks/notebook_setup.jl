import Pkg
Pkg.activate(".")
Pkg.instantiate()
using WebIO
WebIO.install_jupyter_nbextension()
using DataFrames
using DataFramesMeta
# using Gadfly # Removed for now because downgrades DataFrames version
using IJulia
using Weave
using CSV
using Revise
using TableView
using Printf
using PrettyTables
using Formatting

# ==================== GENERAL UTILITY FUNCTIONS ====================

broadwrap(f) = function (args...) broadcast(f, args...) end
esc_latex(s) = replace(s, r"(#|%|&)" => s"\\\1")
function number2latex(num)
		if ismissing(num)
				"--"
		elseif isa(num, Integer)
				"\\(" * sprintf1("%'d", num) * "\\)"
		elseif isa(num, AbstractFloat)
				"\\(" * sprintf1("%'.2f", num) * "\\)"
		elseif isa(num, NTuple{2, Number})
        "$(number2latex(num[1])) ($(number2latex(num[2])))"
		else
				error("Unexpected type.")
		end
end

# ==================== METHODS FOR highlight_best_values! ====================
wrap_in_textbf(str) = "\\textbf{$str}"
# See https://www.pcre.org/current/doc/html/pcre2syntax.html#SEC2
const MATH_MODE_REGEX = Regex(
    "^(?:\\Q\$\\E|\\Q\\(\\E)(.+)(?:\\Q\$\\E|\\Q\\)\\E)\$"
)
function latex2number(str, failure = NaN) :: Float64
    v = tryparse(Float64, str)
    if isnothing(v)
        m = match(MATH_MODE_REGEX, str)
        isnothing(m) || (v = tryparse(Float64, m.captures[1]))
    end
    return isnothing(v) ? failure : v
end
find_index_min(a) = isempty(a) ? Int[] : Int[findmin(a)[2]]

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

