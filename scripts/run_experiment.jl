#!/bin/julia

# Not directly used in this script, but needed to enable the glue code in
# GuillotineModels and make it able to use CPLEX with Requires.jl.
import CPLEX

import Dates
import Dates: @dateformat_str

#= Mock GuillotineModels.run used to print debug the script
module GuillotineModels
	module CommandLine
		run(s) = println(s)
	end
end
=#
import GuillotineModels

function save_output_in_path(f :: Function, path :: String)
	open(path, "a+") do file
		redirect_stdout(file) do
			redirect_stderr(file) do
				f()
			end
		end
	end
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
	, options        :: Vector{String} = String[]
	, solver_seeds   :: Vector{Int} = Int[]
	, output_path    :: String = Dates.format(
		Dates.now(), dateformat"Y-m-dTH:M:S"
	) * ".log"
)
	# The mock should not print anything, and even if it prints something,
	# it will not go into the file (it is not entirely supressed to make
	# the failure of these assumptions visible).
	mock_output_flags = ["--no-csv-output", "--$solver-no-output"]
	mock_flags = append!(["--do-not-mock-for-compilation"], mock_output_flags)
	mock_args = [model, solver, mock_inst_path, mock_flags..., options...]
	GuillotineModels.CommandLine.run(mock_args)
	save_output_in_path(output_path) do
		if isempty(solver_seeds)
			for inst_path in instance_paths
				args = append!([model, solver, inst_path], options)
				GuillotineModels.CommandLine.run(
					args; supported_solvers = [Symbol(solver)]
				)
			end
		else
			for inst_path in instance_paths
				args = append!([model, solver, inst_path], options)
				for seed in solver_seeds
					append!(["--$solver-seed $seed"], args)
					GuillotineModels.CommandLine.run(
						args; supported_solvers = [Symbol(solver)]
					)
				end
			end
		end
	end
end

function run_experiments(
	solver = "CPLEX"
	, instance_folder :: String = "../instances/"
	, mock_inst_name :: String = "gcut1"
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

	tl = 3600.0

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
			append!(options, common_flags)
			run_batch("PPG2KP", solver, mock_inst_path, all_instance_paths, options)
		end
	end

	let # just to avoid the variable names from leaking
		common_options = ["--$solver-time-limit $tl"]
		option_sets = [
			String[], ["--PPG2KP-round2disc"],
			["--PPG2KP-no-cut-position", "--PPG2KP-no-redundant-cut"]
		]
		solver_seeds = collect(1:10)
		for options in option_sets
			append!(options, common_options)
			run_batch(
				"PPG2KP", solver, mock_inst_path, all_instance_paths, options,
				solver_seeds
			)
			run_batch(
				"PPG2KP", solver, mock_inst_path, easy_instance_paths,
				append!(String["--PPG2KP-faithful2furini2016"], options), solver_seeds
			)
		end
	end
end

run_experiments()

