"""
Write the model to a file as a CSV.
"""
function generateCSV(modelName, sol)
  df1 = DataFrames.DataFrame(sol)
  DataFrames.rename!(df1, Dict(:timestamp => "time"))
  modelName = replace(modelName, "."=>"_")
  finalFileName = string("./Data/" * modelName,"julia_res.csv")
  CSV.write(finalFileName, df1)
end

"""
  Test runs ESCIMO without using special initialization values.
"""
function compileESCIMO(tspan = (1850., 2500.))
  @assert(first(tspan) == 1850., "Simulation should start at 1850.")
  (ESCIMOModel_problem, callbacks, ivs, ESCIMOModel_ReducedSystem, tspan, pars, vars, irreductable) = ESCIMOModel(tspan)
end

"""
  Runs the escimo file using precalculated initial conditions.
  See comments below.
"""
function simulate(;tspan = (1850., 2500.), solver = DifferentialEquations.Rodas5(autodiff=false))
  println(tspan)
  @assert(first(tspan) == 1850., "Simulation should start at 1850.")
  println("Creating model...")
  (ESCIMOModel_problem, callbacks, ivs, ESCIMOModel_ReducedSystem, tspan, pars, vars, irreductable) = ESCIMOModel(tspan)
  #=
  Currently (2023-05-18), MTK has some issues with calculating initial values.
  For now we read the initial value from an external file.
  =#
  println("Reading initial values from external file....")
  df = CSV.read("./initial.csv", DataFrames.DataFrame)
  initial_val_dict = esci_dict = Dict(pairs(eachcol(df)))
  ks = map([keys(initial_val_dict)...]) do x
    string(x)
  end
  vals = map(x->first(x), values(initial_val_dict))
  initialValDict2 = Dict(zip(ks, vals))
  newStartVals = Pair{Symbolics.Num}[]
  #= Get initial values from OpenModelica =#
  for (i,p) in enumerate(ivs)
    name = replace(string(first(p)), "(t)" => "")
    if name in keys(initialValDict2)
      valInAux = initialValDict2[name]
      println("New value for $name:" * string(valInAux))
      @assign p.second = valInAux
      ivs[i] = p
    end
  end
  println("Reinitialize the problem...")
  println("Value of tspan... $(tspan)")
  modifiedProblem = ModelingToolkit.ODEProblem(ESCIMOModel_ReducedSystem, ivs, tspan, pars, callback = callbacks)
  println("Simulating....")
  tstops = [i for i in 1850:10:2500]
  sol = DifferentialEquations.solve(modifiedProblem, solver; tstops = tstops)
  println("Simulation done...")
  println("sol.t: $(sol.t)")
  println("sol.retcode: $(sol.retcode)")
  generateCSV("ESCIMO", sol)
  sol
end
