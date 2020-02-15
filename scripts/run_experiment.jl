#!/bin/julia

module GuillotineModels
	module CommandLine
		run(s) = println(s)
	end
end
#import GuillotineModels

# The concept of a batch is a set of runs that share the same model, solver,
# and most options but may vary the instances and the random seeds.
# The main utility of a batch is that we solve a mock instance before the
# other instances and then we have guarante that everything that depends
# on model, solver, and most options will be already compiled by JIT before
# we make serious time measurements.
function run_batch(
	model            :: String
	, solver         :: String
	, mock_inst_path :: String
	, instance_paths :: AbstractVector{String}
	, options        :: String = ""
	, solver_seeds   :: Vector{Int} = Int[]
)
	mock_output_flags = "--no-csv-output --$solver-no-output"
	mock_flags = "--do-not-mock-for-compilation $mock_output_flags"
	mock_args = "$model $solver $mock_inst_path $mock_flags $options"
	GuillotineModels.CommandLine.run(mock_args)
	if isempty(solver_seeds)
		for inst_path in instance_paths
			args = "$model $solver $inst_path $options"
			GuillotineModels.CommandLine.run(args)
		end
	else
		for inst_path in instance_paths
			args = "$model $solver $inst_path $options"
			for seed in solver_seeds
				args *= " --$solver-seed $seed"
				GuillotineModels.CommandLine.run(args)
			end
		end
	end
end

function run_experiments(
	solver = "CPLEX"
	, instance_folder :: String = "../intances/"
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
		common_flags = "--do-not-solve --faithful2furini2016"
		option_sets = [
			""
			, "--no-redundant-cut --no-cut-position"
			, "--no-redundant-cut"
			, "--no-cut-position"
			, "--round2disc"
		]
		for options in option_sets
			run_batch("PPG2KP", solver, mock_inst_path, all_instance_paths, options)
		end
	end

	let # just to avoid the variable names from leaking
		common_options = "--$solver-time-limit $tl"
		option_sets = [
			"", " --round2disc" , " --no-cut-position --no-redundant-cut"
		]
		option_sets = common_options .* option_sets
		solver_seeds = collect(1:10)
		for options in option_sets
			run_batch(
				"PPG2KP", solver, mock_inst_path, all_instance_paths, options,
				solver_seeds
			)
			run_batch(
				"PPG2KP", solver, mock_inst_path, easy_instance_paths,
				options * " --faithful2furini2016", solver_seeds
			)
		end
	end
end

run_experiments()

