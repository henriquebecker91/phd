#!/bin/bash
# -*- mode: julia -*-
#=
exec julia --project=@. --color=yes --startup-file=no -e "include(popfirst!(ARGS))" "${BASH_SOURCE[0]}" "$@"
=#

import DelimitedFiles

# ==================== GATHERING METHODS ====================

function gather_data_from_files(
	@nospecialize(filenames), @nospecialize(extractors)
)
	file_contents = String[]
	sizehint!(file_contents, length(filenames))
	for filename in filenames
		push!(file_contents, read(filename, String))
	end
	processed_data = []
	sizehint!(processed_data, length(extractors))
	for e in extractors
		push!(processed_data, map(e, file_contents))
	end
	return processed_data
end

function gather_data_from_folder(
	@nospecialize(folder_name), @nospecialize(extractors)
)
	filenames = readdir(folder_name; join = true, sort = false)
	return gather_data_from_files(filenames, extractors)
end

function gather_csv_from_folder(
	@nospecialize(folder_name), @nospecialize(extractors);
	delim = ';', column_names = missing
) :: String
	data = gather_data_from_folder(folder_name, extractors)
	iob = IOBuffer()
	if column_names !== missing
		@assert length(column_names) == length(extractors)
		println(iob, join(column_names, delim))
	end
	DelimitedFiles.writedlm(iob, zip(data...), delim)
	return read(seekstart(iob), String)
end

# ==================== MY_PARSE METHODS ====================

# An extensible parsing method.
function my_parse(::Type{String}, str :: AbstractString) :: String
	m = match(r"^\"(.*)\"$", str)
	if isnothing(m)
		return String(str)
	else
		return String(m.captures[1])
	end
end

function my_parse(::Type{N}, str :: AbstractString) :: N where {N <: Number}
	return parse(N, str)
end

# ==================== EXTRACTOR METHODS ====================

struct NoDefault{T} end

struct RegexExtractor{T}
	regex :: Regex
	default :: Union{T, NoDefault{T}}
end

function (re::RegexExtractor{T})(data :: AbstractString) :: T where {T}
	m = match(re.regex, data)
	isnothing(m) || return my_parse(T, only(m.captures))
	if re.default === NoDefault{T}()
		error("The regex $(re.regex) has found no matches. No default available.")
	end
	return re.default
end

function key_equals_extractor(
	key :: AbstractString, default :: Union{T, NoDefault{T}}
) where {T}
	return RegexExtractor{T}(
		Regex("^$key = (.*)\$", "m"), default
	)
end

function p_args_key_extractor(
	key :: AbstractString, default :: Union{T, NoDefault{T}}
) where {T}
	return RegexExtractor{T}(
		Regex("^p_args = .*\"$key\" => ([^,)]*)", "m"), default
	)
end

struct BooleanExtractor
	regex :: Regex
	match_answer :: Bool
end

function (be::BooleanExtractor)(data :: AbstractString)
	m = match(be.regex, data)
	isnothing(m) && return !be.match_answer
	return be.match_answer
end

function matches(r :: Regex)
	return BooleanExtractor(r, true)
end

function does_not_match(r :: Regex)
	return BooleanExtractor(r, false)
end

@enum RNStat Objective=1 Iterations=2 Time=3
Base.getindex(indexable :: AbstractArray, s :: RNStat) = indexable[Int(s)]
Base.getindex(indexable :: Tuple, s :: RNStat) = indexable[Int(s)]

struct RootNodeStatsExtractor
	marker :: String
	stat :: RNStat
end

# Return the objective value, number of iterations, and the time spent
# (in seconds) in the first Gurobi log node relaxation line it finds after
# the `marker` string (the `maker` should be a full line).
function (rnse::RootNodeStatsExtractor)(
	log :: AbstractString
) :: Union{Int, Float64}
	marker = rnse.marker
	# Gets the values in a line starting with "Root relaxation: objective"
	# after a line with only "MARK_$marker" but before the next line starting
	# with "MARK_". Uses non-greedy repetition, non-capture groups and
	# negative look-arounds, so it is basically black magic.
	mark_sentinel = "^MARK_$marker\$(?:(?!^MARK_).)+?"
	# The presolve can remove all rows and columns in some cases, and then
	# the extraction below will fail. So we check this first, and return zero
	# values in this case (that is the most sensible value).
	presolve_line = "^Presolve: All rows and columns removed\$"
	presolve_regex = Regex(mark_sentinel * presolve_line, "ms")
	if !isnothing(match(presolve_regex, log))
		return (0.0, 0, 0.0)[rnse.stat]
	end
	# If barrier is used, the barrier log line appears first and has the most
	# accurate time. If it isn't the root relaxation line is the only one
	# present and has the correct time. Unfortunately the order of the captures
	# in this segment does not match the RNStat enum.
	barrier_log = "^Barrier solved model in (\\d+) iterations and (\\S+) " *
		"seconds\nOptimal objective (\\S+)"
	# The order of captures in the line below match the enum RNStat.
	root_rel = "^Root relaxation: (?:objective (\\S+)|cutoff), (\\d+) " *
		"iterations, (\\S+) seconds"
	either_line = "(?:$barrier_log|$root_rel)"
	r = Regex(mark_sentinel * either_line, "ms")
	m = match(r, log)
	if nothing === m
		return (NaN, -1, NaN)[rnse.stat]
	else
		@assert 6 == length(m.captures)
		# If there is no barrie line.
		if isnothing(match(r"^Barrier solved model"ms, m.match))
			# A capture has index i if it is the i-esim capture to appear in the
			# string, even if previous captures cannot coexist with the current
			# one (i.e., one is before a '|' inside a group and the other after).
			# This way, all non-barrier captures appear after the barrier ones,
			# even if there were no barrier section (the first three captures
			# are nothing in this case).
			@assert m.captures[[1, 2, 3]] == [nothing, nothing, nothing]
			idx = 3 + Int(rnse.stat) # The enumeration has the right order.
			# The root node relaxation can present the string 'cutoff' instead
			# of the objective value, not sure, but it seems to happen when the
			# warm-start has the same value as the relaxation rounded down. In
			# this case, the script processing the data needs to search the
			# objective value in the restricted_LP key. We will not do this
			# here.
			if rnse.stat == Objective && isnothing(m.captures[idx])
				@assert !isnothing(match(r"^Root relaxation: cutoff,"ms, m.match))
				return NaN
			end
			return parse((Float64, Int, Float64)[rnse.stat], m.captures[idx])
		else # If there is a barrier line.
			@assert m.captures[[4, 5, 6]] == [nothing, nothing, nothing]
			# The enumeration has not the right order in this line.
			idx :: Int = rnse.stat == Iterations ? 1 : rnse.stat == Time ? 2 : 3
			return parse((Int, Float64, Float64)[rnse.stat], m.captures[idx])
		end
	end
end

# ==================== FUNCTION CALLS ====================

# Data extraction for the experiment comparing the revised model with
# the our reimplementation of the original model.
#=
csv = gather_csv_from_folder(
	"./finished_experiments/comparison_2020-07-08T19:53:19/",
	[
		key_equals_extractor("instfname", NoDefault{String}()),
		p_args_key_extractor("PPG2KP-pricing", NoDefault{String}()),
		p_args_key_extractor("PPG2KP-round2disc", NoDefault{String}()),
		p_args_key_extractor("PPG2KP-MIP-start", NoDefault{String}()),
		p_args_key_extractor(
			"PPG2KP-do-not-purge-unreachable", NoDefault{String}()
		),
		p_args_key_extractor(
			"PPG2KP-faithful2furini2016", NoDefault{String}()
		),
		key_equals_extractor("build_and_solve_time", NaN),
		key_equals_extractor("total_instance_time", NaN),
		key_equals_extractor("final_pricing_time", NaN),
		key_equals_extractor("iterative_pricing_time", NaN),
		key_equals_extractor("restricted_final_pricing_time", NaN),
		key_equals_extractor("pricing_time", NaN),
		key_equals_extractor("finished_model_solve", NaN),
		does_not_match(r"TimeoutError"),
		key_equals_extractor("this_data_file", NoDefault{String}())
	];
	column_names = [
		"instance_name",
		"pricing_method",
		"round2disc",
		"mip_start",
		"purge_disabled",
		"faithful",
		"build_and_solve_time",
		"total_instance_time",
		"final_pricing_time",
		"iterated_pricing_time",
		"restricted_pricing_time",
		"total_pricing_time",
		"final_solving_time",
		"finished",
		"datafile"
	]
)
print(csv)
=#

# Data extraction for the experiment related to barrier vs dual simplex
# and their effects in the furini pricing.
#=
csv = gather_csv_from_folder(
	"./finished_experiments/LP_method_2020-07-07T16:53:13/",
	[
		key_equals_extractor("instfname", NoDefault{String}()),
		p_args_key_extractor("Gurobi-LP-method", NoDefault{Int}()),
		p_args_key_extractor(
			"PPG2KP-Gurobi-LP-method-inside-furini-pricing", NoDefault{Int}()
		),
		p_args_key_extractor("PPG2KP-pricing", NoDefault{String}()),
		RootNodeStatsExtractor("FURINI_PRICING_RESTRICTED_MIP_SOLVE", Time),
		RootNodeStatsExtractor("FINAL_GENERIC_SOLVE", Time),
		key_equals_extractor("build_and_solve_time", NaN),
		key_equals_extractor("total_instance_time", NaN),
		key_equals_extractor("final_pricing_time", NaN),
		key_equals_extractor("iterative_pricing_time", NaN),
		key_equals_extractor("restricted_final_pricing_time", NaN),
		key_equals_extractor("pricing_time", NaN),
		does_not_match(r"TimeoutError"),
		key_equals_extractor("this_data_file", NoDefault{String}())
	];
	column_names = [
		"instance_name",
		"lp_method",
		"lp_method_switch",
		"pricing_method",
		"restricted_root_time",
		"final_root_time",
		"build_and_solve_time",
		"total_instance_time",
		"final_pricing_time",
		"iterated_pricing_time",
		"restricted_pricing_time",
		"total_pricing_time",
		"finished",
		"datafile"
	]
)
print(csv)
=#

csv = gather_csv_from_folder(
	"./finished_experiments/faithful_2020-07-06T19:05:19/",
	[
		# Primary key: these four identify a set of parameters and, as there
		# was no repetition, a row of the CSV.
		key_equals_extractor("instfname", NoDefault{String}()),
		p_args_key_extractor("PPG2KP-pricing", NoDefault{String}()),
		p_args_key_extractor("PPG2KP-no-redundant-cut", NoDefault{Bool}()),
		p_args_key_extractor("PPG2KP-no-cut-position", NoDefault{Bool}()),

		# The model stats at multiple stats.
		key_equals_extractor("qt_pevars_after_preprocess", -1),
		key_equals_extractor("qt_cmvars_after_preprocess", -1),
		key_equals_extractor("qt_plates_after_preprocess", -1),

		key_equals_extractor("qt_pevars_restricted", -1),
		key_equals_extractor("qt_cmvars_restricted", -1),
		key_equals_extractor("qt_plates_restricted", -1),

		key_equals_extractor("qt_pevars_priced_restricted", -1),
		key_equals_extractor("qt_cmvars_priced_restricted", -1),
		key_equals_extractor("qt_plates_priced_restricted", -1),

		key_equals_extractor("qt_pevars_before_iterated", -1),
		key_equals_extractor("qt_cmvars_before_iterated", -1),
		key_equals_extractor("qt_plates_before_iterated", -1),

		key_equals_extractor("qt_pevars_after_iterated", -1),
		key_equals_extractor("qt_cmvars_after_iterated", -1),
		key_equals_extractor("qt_plates_after_iterated", -1),

		key_equals_extractor("qt_pevars_after_final", -1),
		key_equals_extractor("qt_cmvars_after_final", -1),
		key_equals_extractor("qt_plates_after_final", -1),

		key_equals_extractor("length_pe_after_pricing", -1),
		key_equals_extractor("length_cm_after_pricing", -1),
		key_equals_extractor("length_pc_after_pricing", -1),

		key_equals_extractor("qt_pevars_after_purge", -1),
		key_equals_extractor("qt_cmvars_after_purge", -1),
		key_equals_extractor("qt_plates_after_purge", -1),

		key_equals_extractor("build_stop_reason", "NOT_REACHED"),
		does_not_match(r"TimeoutError"),
		key_equals_extractor("this_data_file", NoDefault{String}())
	];
	column_names = [
		# primary key
		"instance_name",
		"pricing_method",
		"disabled_redundant_cut",
		"disabled_cut_position",
		# The model stats at multiple stats.
		"qt_pevars_after_preprocess",
		"qt_cmvars_after_preprocess",
		"qt_plates_after_preprocess",

		"qt_pevars_restricted",
		"qt_cmvars_restricted",
		"qt_plates_restricted",

		"qt_pevars_priced_restricted",
		"qt_cmvars_priced_restricted",
		"qt_plates_priced_restricted",

		"qt_pevars_before_iterated",
		"qt_cmvars_before_iterated",
		"qt_plates_before_iterated",

		"qt_pevars_after_iterated",
		"qt_cmvars_after_iterated",
		"qt_plates_after_iterated",

		"qt_pevars_after_final",
		"qt_cmvars_after_final",
		"qt_plates_after_final",

		# The names of the three stats below are normalized.
		"qt_pevars_after_pricing", # "length_pe_after_pricing",
		"qt_cmvars_after_pricing", # "length_cm_after_pricing",
		"qt_plates_after_pricing", # "length_pc_after_pricing",

		"qt_pevars_after_purge",
		"qt_cmvars_after_purge",
		"qt_plates_after_purge",
		# extra info
		"build_stop_reason",
		"finished",
		"datafile"
	]
)
print(csv)

