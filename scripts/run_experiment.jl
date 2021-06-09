#!/bin/bash
# -*- mode: julia -*-
#=
exec julia --project=@. --color=yes --startup-file=no -e "include(popfirst!(ARGS))" "${BASH_SOURCE[0]}" "$@"
=#

# Not directly used in this script, but needed to enable the glue code in
# GuillotineModels and make it able to use the solver with Requires.jl.
import Gurobi
import TimerOutputs
import CPLEX

import Dates
import Dates: @dateformat_str

# If you want to just check if the calls are being assembled right
# comment the import below and uncomment the mock GuillotineModels
# module after it.
import GuillotineModels
# Mock GuillotineModels.run used to check if the script is correct.
#=
module GuillotineModels
	import TimerOutputs
	const TIMER = TimerOutputs.TimerOutput()

	struct TimeoutError <: Exception end
	function Base.showerror(io :: IO, e :: TimeoutError)
		print(io, "mock timeout error")
	end

	module CommandLine
		import ..TimeoutError
		function run(args; kwargs...)
			@show args
			@show kwargs
			if any(s -> !isnothing(match(r".*/Hchl6s", s)), args)
				throw(TimeoutError())
			end
			if any(s -> !isnothing(match(r".*/A5", s)), args)
				error("mock ErrorException")
			end
		end
	end
end
=#

function save_output_in_path(f, path :: AbstractString)
	open(path, "a+") do file
		save_output_in_file(f, file)
	end
end
function save_output_in_file(f, file :: IO)
  redirect_stdout(file) do
		redirect_stderr(file) do
			f()
			# See: discourse.julialang.org/t/cannot-capture-log-with-jump/55365
			# My thanks to Miles Lubin
			Base.Libc.flush_cstdio()
		end
	end
end

macro display_error(io, e)
	return quote
		println("run_ended_by_exception = true")
		showerror($(esc(io)), $(esc(e)))
		println($(esc(io)))
		show($(esc(io)), "text/plain", stacktrace(catch_backtrace()))
		println($(esc(io)))
	end
end

function safe_run(args, supported_solvers, implemented_models, verbose = true)
	try
		GuillotineModels.CommandLine.run(
			args; supported_solvers = supported_solvers,
			implemented_models = implemented_models
		)
		verbose && println("run_ended_by_exception = false")
	catch e
		# If `verbose == false` this is probably a mock run, and we do want to
		# have exception logging in a mock run because no exceptions should
		# happen inside a mock run (not even timeouts).
		if isa(e, GuillotineModels.TimeoutError)
			@display_error stderr e
		elseif isa(e, LoadError) && isa(e.error, GuillotineModels.TimeoutError)
			@display_error stderr e.error
		else
			@display_error stderr e
			open("UNEXPECTED_EXCEPTION_README.log", "a") do io
				println(io, Dates.format(
					Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS"
				))
				println(io, args)
				@display_error io e
			end
		end
	end
	# The GC should happen in a mock run, but not the printing.
	extra_gc_time = @elapsed GC.gc()
	if verbose
		@show extra_gc_time
		TimerOutputs.print_timer(GuillotineModels.TIMER; allocations = false)
		println()
	end
	TimerOutputs.reset_timer!(GuillotineModels.TIMER)
	return
end

function saved_run(
	args               :: Vector{String},
	supported_solvers  :: Vector{Symbol},
	implemented_models :: Vector{Symbol},
	output_folder      :: String,
	verbose            :: Bool = true
)
	format = dateformat"yyyy-mm-ddTHH:MM:SS"
	verbose && println("Next run arguments: $args")
	verbose && println("Started at: ", Dates.format(Dates.now(), format))
	filepath, io = mktemp(output_folder; cleanup = false)
	save_output_in_file(io) do
		_, t, bytes, gctime, _ = @timed safe_run(
			args, supported_solvers, implemented_models, true
		)
		println("run_total_time = $t")
		println("run_bytes_allocated = $bytes")
		println("run_time_spent_on_gc = $gctime")
		println("this_data_file = $filepath")
	end
	close(io)
	verbose && println("Finished at: ", Dates.format(Dates.now(), format))
	println("Output saved to: $(filepath)")
	return filepath
end

# The concept of a batch is a set of runs that share the same model, solver,
# and most options but may vary the instances and the random seeds.
# The main utility of a batch is that we solve a mock instance before the
# other instances and then we have the guarantee that everything that depends
# on model, solver, and most options will be already compiled by JIT before
# we make serious time measurements (theoretically, the change of seed and
# instance should not change the inner methods called).
function run_batch(
	problem           :: String
	, instance_format :: String
	, formulation     :: String
	, solver          :: String
	, mock_inst_path  :: String
	, instance_paths  :: AbstractVector{String}
	; options         :: Vector{String} = String[]
	, solver_seeds    :: Vector{Int} = Int[]
	, output_folder   :: String = "./"
)
	supported_solvers = [Symbol(solver)]
	implemented_models = [Symbol(formulation)]
	# The mock should not print anything, and even if it prints something,
	# it will not go into the file (it is not entirely supressed to make
	# the failure of these assumptions visible so we can correct them).
	mock_output_flags = [
		"--no-csv-output", "--$solver-no-output", "--$formulation-quiet"
	]
	tinkered_options = append!(["--warm-jit", "no"], options)
	mock_flags = vcat(tinkered_options, mock_output_flags)
	mock_args = append!(
		[problem, instance_format, formulation, solver, mock_inst_path],
		mock_flags
	)
	safe_run(mock_args, supported_solvers, implemented_models, false)
	for inst_path in instance_paths
		non_seed_args = append!(
			[problem, instance_format, formulation, solver, inst_path],
			tinkered_options
		)
		if isempty(solver_seeds)
			filename = saved_run(
				non_seed_args, supported_solvers, implemented_models, output_folder
			)
		else
			for seed in solver_seeds
				args_with_seed = append!(["--$solver-seed", "$seed"], non_seed_args)
				filename = saved_run(
					args_with_seed, supported_solvers, implemented_models, output_folder
				)
			end
		end
	end
	return
end

#=
function run_experiments(
	solver = "CPLEX"
	; instance_folder :: String = "../instances/"
	, output_folder   :: String = "./experiments_outputs/" * Dates.format(
		Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS"
	)
)
	isdir(output_folder) || mkpath(output_folder)
	all_instance_names = vcat(
		"gcut" .* string.(collect(2:12))
		, ["STS4", "STS4s", "A5"]
		, "Hchl" .* split("2 3s 4s 6s 7s")
		, "CW" .* split("1 2 3")
		, "CHL" .* split("1 1s 6 7")
		, "CU" .* split("1 2")
		, "okp" .* split("1 2 3 4 5")
		, "P1_100_200_25_" .* split("1 2 3 4 5")
	)
	all_instance_paths = instance_folder .* all_instance_names
	easy_instance_names = vcat(
		"okp" .* split("1 4 5")
		, ["CU1", "STS4", "STS4s", "gcut9"]
	)
	easy_instance_paths = instance_folder .* easy_instance_names

	tl = 60.0

	#=
	let # just to avoid the variable names from leaking
		common_flags = ["--do-not-solve", "--PPG2KP-faithful2furini2016"]
		option_sets = [
			String[]
			, ["--PPG2KP-no-redundant-cut", "--PPG2KP-no-cut-position"]
			, ["--PPG2KP-no-redundant-cut"]
			, ["--PPG2KP-no-cut-position"]
			, ["--PPG2KP-round2disc"]
		]
		for options in option_sets
			append!(options, common_flags) # note: this changes `option_sets`
			run_batch(
				"PPG2KP", solver, mock_inst_path, all_instance_paths, options;
				output_path = output_path
			)
		end
	end
	=#

	let # just to avoid the variable names from leaking
		common_options = [
			"--generic-time-limit", "$tl", "--PPG2KP-building-time-limit", "$tl",
			"--PPG2KP-verbose"
		]
		#=option_sets = [
			String[], ["--PPG2KP-faithful2furini2016"], #=["--PPG2KP-no-pricing"],
			["--PPG2KP-no-pricing", "--PPG2KP-faithful2furini2016"]=#
			["--PPG2KP-round2disc"], ["--PPG2KP-round2disc",
			"--PPG2KP-faithful2furini2016"]
			#["--PPG2KP-no-cut-position", "--PPG2KP-no-redundant-cut"]
		]=#
		option_sets = [
			["--PPG2KP-pricing", "none"], # revised, no pricing
			String[], # revised, becker pricing
			["--PPG2KP-faithful2furini2016", "--PPG2KP-pricing", "none"], # furini, no pricing
			["--PPG2KP-faithful2furini2016"] # furini, furini pricing
		]
		#solver_seeds = [1] #collect(1:10)
		for options in option_sets
			append!(options, common_options) # NOTE: changes `option_sets` elements
			run_batch(
				"PPG2KP", solver, easy_instance_paths[1], easy_instance_paths;
				options = options, #solver_seeds = solver_seeds,
				output_folder = output_folder
			)
		end
	end
end
=#

# The 7 instances that Furini 2016 method without pricing is able to solve
# within their time limit.
const EASY_SEVEN = vcat(
	"okp" .* split("1 4 5"), ["CU1", "STS4", "STS4s", "gcut9"]
)

function run_LP_method_experiment(
	solver = "CPLEX"
	; instance_folder :: String = "../instances/"
	, output_folder   :: String = "./experiments_outputs/" * Dates.format(
		Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS"
	)
)
	isdir(output_folder) || mkpath(output_folder)
	easy_instance_paths = instance_folder .* EASY_SEVEN

	time_limit = 3600.0

	common_options = [
		"--generic-time-limit", "$time_limit", "--PPG2KP-building-time-limit",
		"$time_limit", "--PPG2KP-verbose", "--PPG2KP-faithful2furini2016"
	]
	option_sets = [
		# No pricing, use dual simplex.
		["--PPG2KP-pricing", "none", "--Gurobi-LP-method", "1"],
		# No pricing, use barrier.
		["--PPG2KP-pricing", "none", "--Gurobi-LP-method", "2"],
		# Furini pricing, use dual simplex.
		["--PPG2KP-pricing", "furini", "--Gurobi-LP-method", "1"],
		# Furini pricing, use barrier.
		["--PPG2KP-pricing", "furini", "--Gurobi-LP-method", "2"],
		# Furini pricing, use barrier outside of the pricing, in the pricing use
		# dual simplex.
		[
			"--PPG2KP-pricing", "furini", "--Gurobi-LP-method", "2",
			"--PPG2KP-Gurobi-LP-method-inside-furini-pricing", "1"
		]
	]
	solver_seeds = [1]
	for options in option_sets
		append!(options, common_options) # NOTE: changes `option_sets` elements
		run_batch(
			# TODO: the mock instance probably should not be hardcoded below
			"PPG2KP", solver, instance_folder * "CU1", easy_instance_paths;
			options = options, solver_seeds = solver_seeds,
			output_folder = output_folder
		)
	end
	return
end

# The 59 instances presented in table 2.1, page 30,
# from DOI: 10.6092/unibo/amsdottorato/7399 (order was kept).
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

function run_faithful_reimplementation_experiment(
	solver = "CPLEX"
	; instance_folder :: String = "../instances/"
	, output_folder   :: String = "./experiments_outputs/" * Dates.format(
		Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS"
	)
)
	isdir(output_folder) || mkpath(output_folder)
	#HARD_SIX = ["gcut11", "gcut12", "okp3", "okp2", "Hchl6s", "Hchl7s"]
	instance_paths = instance_folder .* THOMOPULOS_THESIS_INSTANCES
	time_limit = 10800.0 # three hours

	common_options = [
		"--generic-time-limit", "$time_limit", "--PPG2KP-building-time-limit",
		"$time_limit", "--PPG2KP-verbose", "--PPG2KP-faithful2furini2016",
		"--do-not-solve"
	]
	option_sets = [
		# The Complete PP-G2KP (no reductions).
		["--PPG2KP-pricing", "none", "--PPG2KP-no-redundant-cut",
			"--PPG2KP-no-cut-position"],
		# The Complete PP-G2KP (only cut position).
		["--PPG2KP-pricing", "none", "--PPG2KP-no-redundant-cut"],
		# The Complete PP-G2KP (only redundant cut).
		["--PPG2KP-pricing", "none", "--PPG2KP-no-cut-position"],
		# The Complete PP-G2KP (no reductions, no pricing).
		["--PPG2KP-pricing", "none"],
		# The PP-G2KP (i.e., with pricing).
		["--PPG2KP-pricing", "furini",
			"--PPG2KP-Gurobi-LP-method-inside-furini-pricing", "1"]
	]
	solver_seeds = [1]
	for options in option_sets
		append!(options, common_options) # NOTE: changes `option_sets` elements
		# gcut6 was chosen as the mock instance because it has median
		# times in all option sets, so it probably enter the longer
		# computation paths on all of them.
		run_batch(
			"PPG2KP", solver, instance_folder * "gcut6", instance_paths;
			options = options, solver_seeds = solver_seeds,
			output_folder = output_folder
		)
	end
	return
end

function run_comparison_experiment(
	; instance_folder :: String = "../instances/"
	, output_folder   :: String = "./experiments_outputs/" * Dates.format(
		Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS"
	)
)
	isdir(output_folder) || mkpath(output_folder)
	instance_paths = instance_folder .* THOMOPULOS_THESIS_INSTANCES
	#instance_paths = [instance_folder .* "A5"]
	time_limit = 10800.0 # three hours

	common_options = [
		"--generic-time-limit", "$time_limit", "--PPG2KP-building-time-limit",
		"$time_limit", "--PPG2KP-verbose"
	]
	solver_options = Dict{String, Vector{String}}(
		"CPLEX" => ["--CPLEX-root-relax-method", "barrier"],
		"Gurobi" => [
			"--Gurobi-LP-method", "2" #= 2 == barrier =#,
			"--PPG2KP-Gurobi-LP-method-inside-furini-pricing", "1" #= 1 == dual =#
		]
	)
	option_sets = [
		# Just the revised model.
		#["--PPG2KP-pricing", "none"],
		# The revised model with our reduction.
		#["--PPG2KP-pricing", "none", "--PPG2KP-round2disc"],
		# The revised model with our reduction and warm-start.
		#["--PPG2KP-pricing", "none", "--PPG2KP-round2disc",
		#	"--PPG2KP-MIP-start", "guaranteed"],
		# The revised model with our reduction, pricing, and warm-start.
		#["--PPG2KP-pricing", "furini", "--PPG2KP-round2disc",
		#	"--PPG2KP-MIP-start", "guaranteed"],
		# The revised model with our reduction, pricing, and warm-start,
		# but removal of unreachable disabled.
		#["--PPG2KP-pricing", "furini", "--PPG2KP-round2disc",
		#	"--PPG2KP-MIP-start", "guaranteed",
		#	"--PPG2KP-do-not-purge-unreachable"],
		# The original model (no pricing, just their reductions).
		["--PPG2KP-faithful2furini2016", "--PPG2KP-pricing", "none"],
		# The original model with our reduction too.
		["--PPG2KP-faithful2furini2016", "--PPG2KP-pricing", "none",
			"--PPG2KP-round2disc"],
		# The orifinal model with our reduction and warm-start.
		["--PPG2KP-faithful2furini2016", "--PPG2KP-pricing", "none",
			"--PPG2KP-round2disc", "--PPG2KP-MIP-start", "guaranteed"],
		# The original model with furini pricing (Priced PPG2KP).
		#["--PPG2KP-faithful2furini2016", "--PPG2KP-pricing", "furini"]
		# The original model with our reduction plus furini pricing.
		#["--PPG2KP-faithful2furini2016", "--PPG2KP-pricing", "furini",
		#	"--PPG2KP-round2disc"]
	]
	solver_seeds = [1]#, 2, 3]
	for solver in [#="CPLEX",=# "Gurobi"]
		for options in option_sets
			append!(options, common_options) # NOTE: changes `option_sets` elements
			# gcut6 was chosen as the mock instance because it has median
			# times in all option sets, so it probably enter the longer
			# computation paths on all of them.
			run_batch(
				"PPG2KP", solver,
				instance_folder * "gcut6",
				instance_paths;
				options = vcat(options, solver_options[solver]),
				solver_seeds = solver_seeds,
				output_folder = output_folder
			)
		end
	end
	return
end

const VELASCO_AND_UCHOA_80 = String[
	"P$(i)_$(L)_$(W)_$(n)_$(s)" for i in 1:4 for (L, W) in
	([(100, 200), (100, 400)], [(200, 100), (400, 100)],
		[(150, 150), (250, 250)], [(150, 150), (250, 250)])[i] for n in (25, 50)
	for s in 1:5
]

function run_vel_uchoa_experiment(
	; instance_folder :: String = "../instances/"
	, output_folder   :: String = "./experiments_outputs/" * Dates.format(
		Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS"
	)
)
	isdir(output_folder) || mkpath(output_folder)
	instance_paths = instance_folder .* VELASCO_AND_UCHOA_80
	time_limit = 10800.0 # three hours

	common_options = [
		"--generic-time-limit", "$time_limit", "--PPG2KP-building-time-limit",
		"$time_limit", "--PPG2KP-verbose"
	]
	solver_options = Dict{String, Vector{String}}(
		#"CPLEX" => ["--CPLEX-root-relax-method", "barrier"],
		"Gurobi" => [
			"--Gurobi-LP-method", "2" #= 2 == barrier =#,
			"--PPG2KP-Gurobi-LP-method-inside-furini-pricing", "1" #= 1 == dual =#
		]
	)
	option_sets = [
		# Just the revised model.
		#["--PPG2KP-pricing", "none"],
		# The revised model with our reduction.
		#["--PPG2KP-pricing", "none", "--PPG2KP-round2disc"],
		# The revised model with our reduction, pricing, and warm-start.
		["--PPG2KP-pricing", "furini", "--PPG2KP-round2disc",
			"--PPG2KP-MIP-start", "guaranteed"],
		# The revised model with our reduction and warm-start.
		["--PPG2KP-pricing", "none", "--PPG2KP-round2disc",
			"--PPG2KP-MIP-start", "guaranteed"],
		# The revised model with our reduction, pricing, and warm-start,
		# but removal of unreachable disabled.
		#["--PPG2KP-pricing", "furini", "--PPG2KP-round2disc",
		#	"--PPG2KP-MIP-start", "guaranteed",
		#	"--PPG2KP-do-not-purge-unreachable"],
		# The original model (no pricing, just their reductions).
		#["--PPG2KP-faithful2furini2016", "--PPG2KP-pricing", "none"],
		# The original model with our reduction too.
		#["--PPG2KP-faithful2furini2016", "--PPG2KP-pricing", "none",
		#	"--PPG2KP-round2disc"],
		# The orifinal model with our reduction and warm-start.
		#["--PPG2KP-faithful2furini2016", "--PPG2KP-pricing", "none",
		#	"--PPG2KP-round2disc", "--PPG2KP-MIP-start", "guaranteed"],
		# The original model with furini pricing (Priced PPG2KP).
		#["--PPG2KP-faithful2furini2016", "--PPG2KP-pricing", "furini"]
		# The original model with our reduction plus furini pricing.
		#["--PPG2KP-faithful2furini2016", "--PPG2KP-pricing", "furini",
		#	"--PPG2KP-round2disc"]
	]
	solver_seeds = [1]#, 2, 3]
	for solver in [#="CPLEX",=# "Gurobi"]
		for options in option_sets
			append!(options, common_options) # NOTE: changes `option_sets` elements
			# gcut6 was chosen as the mock instance because we have no idea yet which
			# instance from the velasco and uchoa is a good mock
			run_batch(
				"PPG2KP", solver,
				instance_folder * "gcut6",
				instance_paths;
				options = vcat(options, solver_options[solver]),
				solver_seeds = solver_seeds,
				output_folder = output_folder
			)
		end
	end
	return
end

# DATASET C OF https://doi.org/10.1016/j.ejor.2010.01.039
const DATASET_C = vcat(
	"A" .* string.(1:5),
	"CHL" .* string.([1, 2, 5, 6, 7]),
	"CU" .* ["1", "2"],
	"CW" .* string.(1:3),
	"Hchl" .* ["2", "9"],
	"Hchl" .* string.([3, 4, 6, 7, 8]) .* "s",
	"HH",
	"OF" .* ["1", "2"],
	"STS" .* ["2", "4"],
	"W",
	"2",
	"3"
)

# DATASET D OF https://doi.org/10.1016/j.ejor.2010.01.039
const DATASET_D = "ATP" .* string.(30:49)

const GCUTS = "gcut" .* string.(1:12)

function run_lagos_experiment(
	all_instances :: AbstractVector{String}
	; instance_folder :: String = "../instances/"
	, output_folder   :: String = "./experiments_outputs/" * Dates.format(
		Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS"
	)
)
	isdir(output_folder) || mkpath(output_folder)
	instance_paths = instance_folder .* all_instances
	time_limit = 10800.0 # three hours

	common_options = [
		"--generic-time-limit", "$time_limit", "--PPG2KP-building-time-limit",
		"$time_limit", "--PPG2KP-verbose"
	]
	option_sets = [
		# The revised model, multiple-threads.
		[ # FOR YOJIMBO
			"--PPG2KP-pricing", "none", "--PPG2KP-round2disc",
			"--Gurobi-threads", "12", "--Gurobi-LP-method", "3", # 3 == concurrent
    ],
		# The revised model, single-thread.
		[ # FOR YOJIMBO
			"--PPG2KP-pricing", "none", "--PPG2KP-round2disc",
			"--Gurobi-LP-method", "2", # 2 == barrier
		],
		# The original model, multiple-threads.
		[ # FOR LEVIATHAN
			"--PPG2KP-pricing", "none",
			"--PPG2KP-faithful2furini2016",
			"--Gurobi-threads", "12", "--Gurobi-LP-method", "3", # 3 == concurrent
    ],
		# The original model, single-thread.
		[ # FOR LEVIATHAN
			"--PPG2KP-pricing", "none",
			"--PPG2KP-faithful2furini2016",
			"--Gurobi-LP-method", "2", # 2 == barrier
		],
		# The original model with plate-size normalization, multiple-threads.
		[ # FOR ODIN
			"--PPG2KP-pricing", "none", "--PPG2KP-round2disc",
			"--PPG2KP-faithful2furini2016",
			"--Gurobi-threads", "12", "--Gurobi-LP-method", "3", # 3 == concurrent
    ],
		# The original model with plate-size normalization, single-thread.
		[ # FOR ODIN
			"--PPG2KP-pricing", "none", "--PPG2KP-round2disc",
			"--PPG2KP-faithful2furini2016",
			"--Gurobi-LP-method", "2", # 2 == barrier
		],
	]
	solver_seeds = [1]
	for options in option_sets
		append!(options, common_options) # NOTE: changes `option_sets` elements
		run_batch(
			"PPG2KP", "Gurobi",
			instance_folder * "STS4", # This is the mock instance.
			instance_paths;
			options = options,
			solver_seeds = solver_seeds,
			output_folder = output_folder
		)
	end
	return
end

const CWs = "CW" .* string.(1:11)
const CUs = "CU" .* string.(1:11)

function run_rotation_experiment(
	solver :: String,
	all_instances :: AbstractVector{String},
	mock_instance :: String,
	; instance_folder :: String = "../instances/"
	, output_folder   :: String = "./experiments_outputs/" * Dates.format(
		Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS"
	)
)
	isdir(output_folder) || mkpath(output_folder)
	instance_paths = instance_folder .* all_instances
	time_limit = 7200.0 # two hours

	common_options = [
		"--generic-time-limit", "$time_limit", "--PPG2KP-building-time-limit",
		"$time_limit", "--PPG2KP-verbose", "--PPG2KP-round2disc",
		"--PPG2KP-pricing", "none"
	]
	option_sets = [
		["--PPG2KP-allow-rotation"]
	]
	@assert solver in ("Gurobi", "CPLEX")
	if solver == "Gurobi"
		append!.(option_sets, (["--Gurobi-LP-method", "2"],)) # barrier
	else
		append!.(option_sets, (["--CPLEX-LP-method", "4"],)) # barrier
	end
	solver_seeds = [1]
	for options in option_sets
		append!(options, common_options) # NOTE: changes `option_sets` elements
		run_batch(
			"G2KP", "Classic_G2KP", "PPG2KP", solver,
			instance_folder * mock_instance, # This is the mock instance.
			instance_paths;
			options = options,
			solver_seeds = solver_seeds,
			output_folder = output_folder
		)
	end
	return
end

function save_models(
	all_instances :: AbstractVector{String},
	mock_instance :: String,
	; instance_folder :: String = "../instances/"
	, output_folder   :: String = "./experiments_outputs/" * Dates.format(
		Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS"
	),
	extension = "mps",
	mps_folder = joinpath("./$extension/", last(splitpath(instance_folder)))
)
	# extension should be ascii utf8 so this should be safe
	@assert extension[1] != '.'
	isdir(output_folder) || mkpath(output_folder)
	isdir(mps_folder) || mkpath(mps_folder)
	instance_paths = joinpath.((instance_folder,), all_instances)

	common_options = [
		"--PPG2KP-verbose", "--PPG2KP-pricing", "none",
		"--do-not-solve", "--save-model", joinpath(
			mps_folder, "enhanced2-\$<instance_file>" *
			"-\$<PPG2KP-allow-rotation>.$extension"
		)
	]
	option_sets = [
		["--PPG2KP-round2disc", "--PPG2KP-allow-rotation",
			"--PPG2KP-mirror-plate"],
		#["--PPG2KP-faithful2furini2016", "--PPG2KP-allow-rotation",
		#	"--PPG2KP-mirror-plate"],
		["--PPG2KP-round2disc"], # no rotation
		#["--PPG2KP-faithful2furini2016"], # no rotation
	]
	for options in option_sets
		append!(options, common_options) # NOTE: changes `option_sets` elements
		run_batch(
			"G2KP", "Simple_CPG_SLOPP", "PPG2KP", "NoSolver",
			instance_folder * mock_instance, # This is the mock instance.
			instance_paths; options = options, output_folder = output_folder
		)
	end
	return
end

const HARD4 = String[
	"Hchl4s", "Hchl7s", "okp2", "okp3"
]

function run_hybridization_experiment(
	solver :: String,
	all_instances :: AbstractVector{String},
	mock_instance :: String,
	; instance_folder :: String = "../instances/"
	, output_folder   :: String = "./experiments_outputs/" * Dates.format(
		Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS"
	)
)
	isdir(output_folder) || mkpath(output_folder)
	instance_paths = instance_folder .* all_instances
	time_limit = 3600.0 * 3 # three hours

	common_options = [
		"--generic-time-limit", "$time_limit", "--PPG2KP-building-time-limit",
		"$time_limit", "--PPG2KP-verbose", "--PPG2KP-round2disc",
		"--PPG2KP-pricing", "none"
	]
	option_sets = [
		String[], # Only common options
		["--PPG2KP-hybridize-with-restricted"],
		["--PPG2KP-hybridize-with-restricted",
			"--PPG2KP-aggressive-hybridization"],
	]
	@assert solver in ("Gurobi", "CPLEX")
	if solver == "Gurobi"
		append!.(option_sets, (["--Gurobi-LP-method", "2"],)) # barrier
	else
		append!.(option_sets, (["--CPLEX-LP-method", "4"],)) # barrier
	end
	solver_seeds = [1]
	for options in option_sets
		append!(options, common_options) # NOTE: changes `option_sets` elements
		run_batch(
			"G2KP", "Classic_G2KP", "PPG2KP", solver,
			instance_folder * mock_instance, # This is the mock instance.
			instance_paths;
			options = options,
			solver_seeds = solver_seeds,
			output_folder = output_folder
		)
	end
	return
end

const Clautiaux42 = String[
	"E00N10", "E00N15", "E00N23", "E00X23", "E02F17", "E02F20", "E02F22",
	"E02N20", "E03N10", "E03N15", "E03N16", "E03N17", "E03X18", "E04F15",
	"E04F17", "E04F19", "E04F20", "E04N15", "E04N17", "E04N18", "E05F15",
	"E05F18", "E05F20", "E05N15", "E05N17", "E05X15", "E07F15", "E07N10",
	"E07N15", "E07X15", "E08F15", "E08N15", "E10N10", "E10N15", "E10X15",
	"E13N10", "E13N15", "E13X15", "E15N10", "E15N15", "E20F15", "E20X15"
]

const HopperTurton_C = String[
	"c$(c)-p$(p)" for c in 1:7 for p in 1:3
]

function run_G2OPP_exploratory_experiment(
	solver :: String,
	all_instances :: AbstractVector{String},
	mock_instance :: String,
	; instance_folder :: String = "../instances/"
	, output_folder   :: String = "./experiments_outputs/" * Dates.format(
		Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS"
	)
)
	isdir(output_folder) || mkpath(output_folder)
	instance_paths = instance_folder .* all_instances
	time_limit = 3600.0 * 3 # three hours

	common_options = [
		"--generic-time-limit", "$time_limit", "--PPG2KP-building-time-limit",
		"$time_limit", "--PPG2KP-verbose", "--PPG2KP-round2disc",
		"--PPG2KP-pricing", "none",	]
	option_sets = [
		String["--Gurobi-raw-parameter",
			"Pair{String, Any}[\"BarHomogeneous\" => 1]"],
		String["--Gurobi-raw-parameter",
			"Pair{String, Any}[\"BarHomogeneous\" => 1, \"Cutoff\" => 1.5]"],
	]
	@assert solver in ("Gurobi", "CPLEX")
	if solver == "Gurobi"
		append!.(option_sets, (["--Gurobi-LP-method", "2"],)) # barrier
	else
		append!.(option_sets, (["--CPLEX-LP-method", "4"],)) # barrier
	end
	solver_seeds = [1]
	for options in option_sets
		append!(options, common_options) # NOTE: changes `option_sets` elements
		run_batch(
			"G2CSP", "CPG_SSSCSP", "PPG2KP", solver,
			instance_folder * mock_instance, # This is the mock instance.
			instance_paths;
			options = options,
			solver_seeds = solver_seeds,
			output_folder = output_folder
		)
	end
	return
end

# The ordering of the instances is the ideal for prioritising solving
# a representative subset of them: Seed > Size > Category.
const CLASS = [
	"cl_$(lpad(c, 2, '0'))_$(lpad(n, 3, '0'))_$(lpad(s, 2, '0'))"
	for s in 1:10 for n in 20:20:100 for c in 1:10
]

const CLASS_50 = [
	"cl_$(lpad(c, 2, '0'))_$(lpad(n, 3, '0'))_01"
	for n in 20:20:100 for c in 1:10 
]

const A = [
	"A42",  "A4",  "A8",  "A3",  "A6", "A39", "A43",  "A1",
	"A22",  "A2", "A30", "A10", "A12", "A34", "A36",  "A7",
	"A26", "A31", "A41", "A37", "A28", "A17",  "A5", "A35",
	"A25", "A38", "A40", "A27", "A32", "A20", "A24", "A29",
	"A23", "A13", "A33", "A15", "A21", "A19", "A18", "A11",
	"A16",  "A9", "A14"
]

const A_MKP = [
	"A02_M2", "A02_M4", "A02_M8", "A03_M2", "A05_M2", "A05_M4", "A07_M2",
	"A07_M4", "A09_M2", "A09_M4", "A09_M8", "A11_M2", "A11_M4", "A11_M8",
	"A12_M2", "A12_M4", "A13_M2", "A13_M4", "A14_M2", "A14_M4", "A14_M8",
	"A15_M2", "A15_M4", "A15_M8", "A16_M2", "A16_M4", "A16_M8", "A17_M2",
	"A18_M2", "A18_M4", "A18_M8", "A19_M2", "A19_M4", "A19_M8", "A20_M2",
	"A20_M4", "A20_M8", "A21_M2", "A21_M4", "A21_M8", "A23_M2", "A23_M4",
	"A24_M2", "A24_M4", "A24_M8", "A25_M2", "A25_M4", "A26_M2", "A27_M2",
	"A27_M4", "A27_M8", "A28_M2", "A28_M4", "A29_M2", "A29_M4", "A29_M8",
	"A31_M2", "A32_M2", "A32_M4", "A32_M8", "A33_M2", "A33_M4", "A33_M8",
	"A34_M2", "A35_M2", "A35_M4", "A36_M2", "A37_M2", "A38_M2", "A38_M4",
	"A38_M8", "A40_M2", "A40_M4", "A41_M2", "A41_M4", "A42_M2", "A43_M2"
]

const CW_MKP = [
	"CW01_M2", "CW01_M4", "CW02_M2", "CW02_M4", "CW03_M2", "CW03_M4",
	"CW04_M2", "CW04_M4", "CW05_M2", "CW05_M4", "CW06_M2", "CW06_M4",
	"CW06_M8", "CW07_M2", "CW07_M4", "CW07_M8", "CW08_M2", "CW08_M4",
	"CW08_M8", "CW09_M2", "CW09_M4", "CW09_M8", "CW10_M2", "CW10_M4",
	"CW10_M8", "CW11_M2", "CW11_M4", "CW11_M8"
]

function run(
	solver          :: String,
	problem         :: String,
	format          :: String,
	all_instances   :: AbstractVector{String},
	mock_instance   :: String,
	instance_folder :: String
	; output_folder :: String = "./experiments_outputs/" * Dates.format(
		Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS"
	)
	, mock_instance_folder = instance_folder
	, MIPGap :: Float64 = 1e-6
	, time_limit :: Float64 = 3600.0
	, Gurobi_LP_method = 2
	, CPLEX_LP_method = 4
	, Gurobi_raw_parameters = String[]
	, CPLEX_raw_parameters = String[]
	, solver_seeds = (1,)
	, option_sets :: Vector{Vector{String}} = [["--PPG2KP-round2disc"]]
)
	isdir(output_folder) || mkpath(output_folder)
	instance_paths = instance_folder .* all_instances

	option_sets = deepcopy(option_sets)

	common_options = [
		"--generic-time-limit", "$time_limit",
		"--PPG2KP-building-time-limit", "$time_limit",
		"--PPG2KP-verbose",
		"--PPG2KP-pricing", "none",
	]
	@assert solver in ("Gurobi", "CPLEX")
	if solver == "Gurobi"
		default_gurobi_raw_params = [
			"\"NumericFocus\" => 3",
			"\"MIPGap\" => $(MIPGap)",
		]
		all_gurobi_raw_params = vcat(
			default_gurobi_raw_params, Gurobi_raw_parameters
		)
		gurobi_raw_params_str = join(string.(all_gurobi_raw_params), ',')
		append!.(option_sets, ([
			"--Gurobi-raw-parameters", "Pair{String, Any}[$(gurobi_raw_params_str)]"
		],))
		append!.(option_sets, (["--Gurobi-LP-method", string(Gurobi_LP_method)],))
	else
		append!.(option_sets, ([
			"--CPLEX-raw-parameters",
			"Pair{String, Any}[\"CPXPARAM_MIP_Tolerances_MIPGap\" => $(MIPGap)]"
		],))
		append!.(option_sets, (["--CPLEX-LP-method", string(CPLEX_LP_method)],))
	end

	for options in option_sets
		append!(options, common_options) # NOTE: changes `option_sets` elements
		run_batch(
			problem, format, "PPG2KP", solver,
			joinpath(mock_instance_folder, mock_instance), # This is the mock instance.
			instance_paths;
			options = options,
			solver_seeds = collect(solver_seeds),
			output_folder = output_folder
		)
	end
	return
end

#=
run(
	"Gurobi", "G2CSP", "CPG_SSSCSP", CLASS, CLASS[1],
	"../instances/G2CSP/CLASS/"; MIPGap = 1e-4
)

run(
	"Gurobi", "G2CSP", "CPG_SSSCSP", A, A[1], "../instances/G2CSP/A/";
	MIPGap = 1e-4
)

run(
	"Gurobi", "G2MKP", "CPG_MHLOPPW", A_MKP, A_MKP[1],
	"../instances/G2MKP/"; MIPGap = 1e-8
)

run(
	"Gurobi", "G2MKP", "CPG_MHLOPPW", CW_MKP, CW_MKP[1],
	"../instances/G2MKP/"; MIPGap = 1e-7
)

run(
	"Gurobi", "G2OPP", "CPG_SSSCSP", Clautiaux42, Clautiaux42[1],
	"../instances/G2OPP/Clautiaux42/"
)

run(
	"Gurobi", "G2OPP", "CPG_SSSCSP", HopperTurton_C, HopperTurton_C[1],
	"../instances/G2OPP/HopperTurton/C/"
)
=#

run(
	"Gurobi", "G2MKP", "CPG_MHLOPPW", A_MKP, A_MKP[1],
	"../instances/G2MKP/"; MIPGap = 1e-8,
	option_sets = [
		["--PPG2KP-round2disc"],
		["--PPG2KP-round2disc", "--PPG2KP-allow-rotation", "--PPG2KP-mirror-plates"]
	]
)

