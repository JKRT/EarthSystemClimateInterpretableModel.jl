module EarthSystemClimateInterpretableModel


import CSV
import DataFrames
import Symbolics
import ModelingToolkit
import DifferentialEquations

using Revise
using MetaModelica

greet() = print("Running EarthSystemClimateInterpretableModel")

include("EscimoModel.jl")
include("util.jl")

end # module EarthSystemClimateInterpretableModel
