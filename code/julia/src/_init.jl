# script to be run in Julia REPL at start of the session making it 
# easier to test bits of code
using Revise # reloading modules now check for changes
using JSON
using JuMP
using Gurobi
push!(LOAD_PATH, "./") # modules below are defined in the current folder
using Chen1995
using Kurpel2018
using Check3DPackings

instances_folder = "../../../instances/hbd_basic_tests/"
instances_fnames = filter(x -> endswith(x, ".json"), readdir(instances_folder))
json_instances = []
for f in instances_fnames
  s = open(x -> read(x, String), instances_folder * f)
  push!(json_instances, JSON.parse(s))
end

function json2model(json;
  no_rotation = false,
  max_packed_volume = false,
  formulation = Chen1995
)
  model = Model(solver = GurobiSolver())
  pqrLWH = map(l -> json[l], split("pqrLWH", ""))

  formulation.build(model, pqrLWH...;
    no_rotation = no_rotation, max_packed_volume = max_packed_volume
  )
end

