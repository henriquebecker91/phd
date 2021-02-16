#!/bin/bash
# -*- mode: julia -*-
#=
exec julia -O3 --project=@. --color=yes --startup-file=no -e "include(popfirst!(ARGS))" "${BASH_SOURCE[0]}" "$@"
=#

import JuMP, CPLEX, Gurobi, MathOptInterface

const MOI = MathOptInterface

# Copyied from save_output_in_path and save_output_in_file
function save_output_in_path(f, path :: AbstractString)
  open(path, "a+") do file
    save_output_in_file(f, file)
  end
end
function save_output_in_file(f, file :: IO)
  redirect_stdout(file) do
		redirect_stderr(file) do
			f()
			Base.Libc.flush_cstdio()
		end
	end
end

# Copyied from GuillotineModels.Utilities
function num_all_constraints(m) :: Int64
	sum = 0 :: Int64
	for (ftype, stype) in JuMP.list_of_constraint_types(m)
		sum += JuMP.num_constraints(m, ftype, stype)
	end
	return sum
end

# Will solve a single instance, and return a list of pair-like objects
# with name of the information and its value.
function read_and_solve_file(
	filepath, # The path to the MPS file.
	optimizer, # CPLEX.Optimizer, Gurobi.Optimizer, etc...
	optimizer_conf, # A list of pairs name/value to be set as solver parameters.
	optimizer_out; # to where to redirect the solver output
	relax = false # If JuMP.relax_integrality is called before solving.
)
	model = JuMP.read_from_file(filepath)
	JuMP.set_optimizer(model, optimizer)
	for (name, value) in optimizer_conf
		JuMP.set_optimizer_attribute(model, name, value)
	end
	#silent && JuMP.set_optimizer_attribute(model, MOI.Silent(), true)
	save_output_in_file(optimizer_out) do
		relax && JuMP.relax_integrality(model)
		JuMP.optimize!(model)
	end

	# TODO: WE NEED TO SAVE TOTAL NUMBER OF VARIABLES AND CONSTRAINTS TOO
	info = Pair{Symbol,Any}[
		:termination_status => JuMP.termination_status(model),
		:raw_status => JuMP.raw_status(model),
		:primal_status => JuMP.primal_status(model),
		:dual_status => JuMP.dual_status(model),
		:simplex_iterations => MOI.get(model, MOI.SimplexIterations()),
		:number_of_variables => MOI.get(model, MOI.NumberOfVariables()),
		:number_of_constraints => num_all_constraints(model),
	]
	if JuMP.has_values(model)
		append!(info, [
			:objective_bound => JuMP.objective_bound(model),
			:objective_value => JuMP.objective_value(model),
			:dual_objective_value => JuMP.dual_objective_value(model),
			:solve_time => MOI.get(model, MOI.SolveTime()),
		])
	end

	if !relax
		append!(info, [
			:node_count => MOI.get(model, MOI.NodeCount()),
			:relative_gap => MOI.get(model, MOI.RelativeGap()),
		])
	end

	return info
end

const OPTIMIZER = Dict{Symbol, Any}([
	:CPLEX => CPLEX.Optimizer,
	:Gurobi => Gurobi.Optimizer,
])

function print_input(
	io, filepath, optimizer_id, optimizer_conf, relax
)
	filename = basename(filepath)
	without_ext, _ = splitext(filename)
	formulation, instance_name, rotation = split(without_ext, '-')
	println(io, "formulation = ", formulation)
	println(io, "instance_name = ", instance_name)
	println(io, "rotation = ", rotation)
	println(io, "filename = ", filename)
	println(io, "filepath = ", filepath)
	println(io, "solver = ", optimizer_id)
	println(io, "relax = ", convert(Int, relax))
	for (name, value) in optimizer_conf
		println(io, name, " = ", value)
	end
	return
end

function print_result(io, result)
	for (name, value) in result
		println(io, name, " = ", value)
	end
	return
end

function batch_read_solve_print(
	output_folder_path,
	filepaths, # A list with the path to each file.
	mock_filepath, # The path to the mock file (can be in filepaths or not).
	optimizer_id, # :CPLEX, :Gurobi, ...
	optimizer_conf # A list of pairs name/value to be set as solver parameters.
)
	# First do the two mock runs with the output disabled and do not
	# save the return.
	optimizer = OPTIMIZER[optimizer_id]
	open("/dev/null", "w+") do devnull
		println("Solving mock relaxation.")
		read_and_solve_file(
			mock_filepath, optimizer, optimizer_conf, devnull; relax = true
		)
		#println("Solving mock MIP.")
		#read_and_solve_file(
		#	mock_filepath, optimizer, optimizer_conf, devnull
		#)
	end
	# Now solve all MPS relaxed and then all MPS as MILP.
	LP_out_path = joinpath(
		output_folder_path, string(optimizer_id), "LP"
	)
	foreach(filepaths) do filepath
		println("Started solving relaxation of model: $(basename(filepath))")
		name, _ = splitext(basename(filepath))
		open(joinpath(LP_out_path, "$name.txt"), "w+") do io
			print_input(
				io, filepath, optimizer_id, optimizer_conf, 1
			)
			flush(io)
			result = read_and_solve_file(
				filepath, optimizer, optimizer_conf, io;
				relax = true, 
			)
			print_result(io, result)
		end
	end
	MIP_out_path = joinpath(
		output_folder_path, string(optimizer_id), "MIP"
	)
	foreach(filepaths) do filepath
		println("Started solving MIP model: $(basename(filepath))")
		name, _ = splitext(basename(filepath))
		open(joinpath(MIP_out_path, "$name.txt"), "w+") do io
			print_input(
				io, filepath, optimizer_id, optimizer_conf, 0
			)
			flush(io)
			result = read_and_solve_file(
				filepath, optimizer, optimizer_conf, io
			)
			print_result(io, result)
		end
	end

	return
end

const CWs = "CW" .* string.(1:11)
const CUs = "CU" .* string.(1:11)
const SET_B = vcat(CWs, CUs)
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

const PARAMETERS = Dict{Symbol, Vector{Pair{String, Any}}}([
	:CPLEX => [
		"CPXPARAM_TimeLimit" => 3600.0,
		"CPXPARAM_Threads" => 1,
	],
	:Gurobi => [
		"TimeLimit" => 3600.0,
		"Threads" => 1,
	],
])

function run_first_experiment(
	models_folder_path,
	output_folder_path,
	formulation,
	ext = "mps"
)
	setB_files = vcat(
		"$formulation-" .* SET_B .* "-0.$ext",
		"$formulation-" .* SET_B .* "-1.$ext",
	)
	setB_paths = joinpath.((models_folder_path,), setB_files)

	SCPG_files = vcat(
		"$formulation-" .* THOMOPULOS_THESIS_INSTANCES .* "-0.$ext",
		"$formulation-" .* THOMOPULOS_THESIS_INSTANCES .* "-1.$ext",
	)
	SCPG_paths = joinpath.((models_folder_path,), SCPG_files)

	filepaths = unique(vcat(SCPG_paths, setB_paths))

	for solver in (:Gurobi, :CPLEX)
		batch_read_solve_print(
			output_folder_path, filepaths, filepaths[1], # the gcut1 instance
			solver, PARAMETERS[solver]
		)
	end

	return
end

run_first_experiment(ARGS...)
