using JSON
using JuMP
using GLPKMathProgInterface
#using Gurobi
using Chen1995

@testset "Chen1995" begin
  inst_dir = "../../../instances/myinstances/"
  test_instances = filter(x -> endswith(x, ".json"), readdir(inst_dir))
  @testset "$f" for f in test_instances
    json = JSON.parse(open(x -> read(x, String), inst_dir * f))
    pqr = map(l -> json[l], split("pqr", ""))
    LWH = map(l -> json[l], split("LWH", ""))
    model = Chen1995.build(Model(solver = GLPKSolverMIP()), pqr..., LWH...)
    status = solve(model, suppress_warnings = true)
    @test string(status) == json["expected_status"]
    if status == :Optimal
      d, xyz, pqr′ = Chen1995.extract_solution(model, pqr)
      @test !has_violations(d, xyz, pqr′, LWH)
    end
  end
end

