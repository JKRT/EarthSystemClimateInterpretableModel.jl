#=
Installation script, author:johti17
=#
import Pkg

@info "Installing the Earth System Climate Interpretable Model"
@assert(Sys.iswindows(), "Windows is currently the only supported operating System.")

pkgsToAdd = [
  "ModelingToolkit",
  "DifferentialEquations",
  "Plots",
  "CSV",
  "DataFrames",
  "Symbolics",
  "Git",
]

@info "Adding dependencies..."
for p in pkgsToAdd
  Pkg.add(p)
end

using Git

#= Clone the external runtime... =#
try
  run(git(["clone", "https://github.com/OpenModelica/OMRuntimeExternalC.jl.git"]))
  @info("Building OMRuntimeExternalC...")
  Pkg.develop(path="./OMRuntimeExternalC.jl")
  @time Pkg.build("OMRuntimeExternalC.jl")
  @info "Developing The Earth System Climate Interpretable Model"
  @time Pkg.develop(path=pwd())
catch e
  println("Installation failed")
  throw(e)
end
