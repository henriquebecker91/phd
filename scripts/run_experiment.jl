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

# Mock GuillotineModels.run used to check if the script is correct.
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
#import GuillotineModels

function save_output_in_path(f, path :: AbstractString)
	open(path, "a+") do file
		redirect_stdout(() -> redirect_stderr(f, file), file)
	end
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
	, output_path    :: String = Dates.format(
		Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS"
	) * ".log"
)
	supported_solvers = [Symbol(solver)]
	implemented_models = [Symbol(model)]
	# The mock should not print anything, and even if it prints something,
	# it will not go into the file (it is not entirely supressed to make
	# the failure of these assumptions visible).
	mock_output_flags = [
		"--no-csv-output", "--$solver-no-output", "--$model-quiet"
	]
	mock_flags = append!(["--warm-jit", "no"], mock_output_flags)
	mock_args = [model, solver, mock_inst_path, mock_flags..., options...]
	safe_run(mock_args, supported_solvers, implemented_models, false)
	save_output_in_path(output_path) do
		for inst_path in instance_paths
			non_seed_args = append!([model, solver, inst_path], options)
			if isempty(solver_seeds)
				safe_run(non_seed_args, supported_solvers, implemented_models)
			else
				for seed in solver_seeds
					args_with_seed = append!(["--$solver-seed", "$seed"], non_seed_args)
					safe_run(args_with_seed, supported_solvers, implemented_models)
				end
			end
		end
	end
	return
end

function run_experiments(
	solver = "CPLEX"
	; instance_folder :: String = "../instances/"
	, mock_inst_name :: String = "gcut1"
	, output_path    :: String = Dates.format(
		Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS"
	) * ".log"
)
	mock_inst_path = instance_folder * mock_inst_name
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
				"PPG2KP", solver, mock_inst_path, easy_instance_paths;
				options = options, #solver_seeds = solver_seeds,
				output_path = output_path
			)
		end
	end
end

#run_experiments("CPLEX")
run_experiments("Gurobi")

