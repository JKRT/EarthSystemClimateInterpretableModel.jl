module EarthSystemClimateInterpretableModel


import CSV
import DataFrames
import Symbolics
import ModelingToolkit
import DifferentialEquations

greet() = print("Running EarthSystemClimateInterpretableModel")

include("EscimoModel.jl")
include("util.jl")

end # module EarthSystemClimateInterpretableModel
