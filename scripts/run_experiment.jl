#!/bin/bash
# -*- mode: julia -*-
#=
exec julia --project=@. --color=yes --startup-file=no -e "include(popfirst!(ARGS))" "${BASH_SOURCE[0]}" "$@"
=#

# Not directly used in this script, but needed to enable the glue code in
# GuillotineModels and make it able to use the solver with Requires.jl.
import Gurobi
import TimerOutputs
#import CPLEX

import Dates
import Dates: @dateformat_str

# If you want to just check if the calls are being assembled right
# comment the import below and uncomment the mock GuillotineModels
# module below.
import GuillotineModels
# Mock GuillotineModels.run used to check if the script is correct.
#=
module GuillotineModels
	import TimerOutputs
	const TIMER = TimerOutputs.TimerOutput()
	module CommandLine
		function run(args...; kwargs...)
			@show args
			@show kwargs
		end
	end
end
=#

function save_output_in_path(f, path :: AbstractString)
	open(path, "a+") do file
		save_output_in_path(f, file)
	end
end

function save_output_in_file(f, file :: IO)
	redirect_stdout(() -> redirect_stderr(f, file), file)
end

function safe_run(args, supported_solvers, implemented_models, verbose = true)
	try
		GuillotineModels.CommandLine.run(
			args; supported_solvers = supported_solvers
		)
	catch e
		# If `verbose == false` this is probably a mock run, and we do want to
		# have exception logging in a mock run because no exceptions should
		# happen inside a mock run (not even timeouts).
		if isa(e, GuillotineModels.TimeoutError)
			showerror(stderr, e)
		elseif isa(e, LoadError) && isa(e.error, GuillotineModels.TimeoutError)
			showerror(stderr, e.error)
		else
			showerror(stderr, e)
			open("UNEXPECTED_EXCEPTION_README.log", "a") do io
				println(io, args)
				println(io, Dates.format(
					Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS"
				))
				showerror(io, e)
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
		safe_run(args, supported_solvers, implemented_models, true)
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
	model            :: String
	, solver         :: String
	, mock_inst_path :: String
	, instance_paths :: AbstractVector{String}
	; options        :: Vector{String} = String[]
	, solver_seeds   :: Vector{Int} = Int[]
	, output_folder  :: String = "./"
)
	supported_solvers = [Symbol(solver)]
	implemented_models = [Symbol(model)]
	# The mock should not print anything, and even if it prints something,
	# it will not go into the file (it is not entirely supressed to make
	# the failure of these assumptions visible so we can correct them).
	mock_output_flags = [
		"--no-csv-output", "--$solver-no-output", "--$model-quiet"
	]
	tinkered_options = append!(["--warm-jit", "no"], options)
	mock_flags = vcat(tinkered_options, mock_output_flags)
	mock_args = append!([model, solver, mock_inst_path], mock_flags)
	safe_run(mock_args, supported_solvers, implemented_models, false)
	for inst_path in instance_paths
		non_seed_args = append!([model, solver, inst_path], tinkered_options)
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

function run_experiments(
	solver = "CPLEX"
	; instance_folder :: String = "../instances/"
	, output_folder   :: String = "./experiments_outputs/" * Dates.format(
		Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS"
	)
)
	isdir(output_folder) || mkpath(output_folder)
	easy_instance_names = vcat(
		"okp" .* split("1 4 5")
		, ["CU1", "STS4", "STS4s", "gcut9"]
	)
	#easy_instance_paths = instance_folder .* easy_instance_names
	# TODO: REMEBER TO CHANGE THIS BACK AFTER, WE JUST WANT TO GET
	# THE VALUES FOR A5 FOR A LITTLE CHECK
	easy_instance_paths = String[instance_folder * "A5"]

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
			"PPG2KP", solver, instance_folder * "gcut1", easy_instance_paths;
			options = options, solver_seeds = solver_seeds,
			output_folder = output_folder
		)
	end
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
	instance_paths = instance_folder .* THOMOPULOS_THESIS_INSTANCES
	time_limit = 3600.0

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
		run_batch(
			"PPG2KP", solver, instance_folder * "gcut1", instance_paths;
			options = options, solver_seeds = solver_seeds,
			output_folder = output_folder
		)
	end
end

#run_experiments("Gurobi")
run_faithful_reimplementation_experiment("Gurobi")
