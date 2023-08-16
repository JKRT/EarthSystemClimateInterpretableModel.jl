# EarthSystemClimateInterpretableModel.jl
The complete [Earth System Climate Interpretable Model (ESCIMO) model](https://esd.copernicus.org/articles/7/831/2016/) implemented in Julia, and made available for Julia scripting.
The model was presented in an article by: Jorgen Randers, Ulrich Golüke, Fred Wenstøp, and Søren Wenstøp in::
*A user-friendly earth system model of low complexity: the ESCIMO system dynamics model of global warming towards 2100*

The model provided here is based upon the model as presented in the Nature Article
*[A user-friendly earth system model of low complexity: the ESCIMO system dynamics model of global warming towards 2100](https://www.nature.com/articles/s41598-020-75481-z)*

## Dependencies
Dependencies include the following Julia packages:
- ModelingToolkit
- DifferentialEquations
- OMRuntimeExternalC
- Plots

### Note
This package makes use of `OMRuntimeExternalC`, which is an external package
that abstracts the C simulation runtime of the OpenModelica Simulation Environment.
Hence, to install and use this package, an installation script is provided see installation.

## Installation
There are two main options for installing and running the model.
Either via manual installation or via an installation script defined in install.jl

### Automatic Installation
To install this package, run install.jl

1. Install Julia
2. Clone/Download this repository
3. Navigate to the folder you downloaded/cloned the software.

```
julia> include("install.jl")
```

### Manual Installation (Advanced)
To manually install this package, make sure you have Julia > 1.9 installed.
Install all the dependent packages as usual.

```
<Navigate to a suitable location to clone OMRuntimeExternalC>
git clone https://github.com/OpenModelica/OMRuntimeExternalC.jl
cd OMRuntimeExternalC
```

Start Julia
```
julia> import Pkg
julia> Pkg.develop("path=<path to OMRuntimeExternalC>")
julia> Pkg.build("OMRuntimeExternalC")
```

## Usage

Example use case:
```
julia> import EarthSystemClimateInterpretableModel
julia> EarthSystemClimateInterpretableModel.simulate(;tspan = (1850., 2500.))
```
Optionally a different solver can be specified using the solver keyword argument:
```
julia> using EarthSystemClimateInterpretableModel
julia> simulate(;tspan = (1850., 2500.); solver = DifferentialEquations.Tsit5())
```

After a successful simulation, a CSV file is generated in the `./Data` folder, and a solution file is returned.
The CSV file can be used in external tools to plot the results. Alternatively, you can further process the returned solution.

```
julia> using EarthSystemClimateInterpretableModel
julia> using Plots
julia> sol = simulate(;tspan = (1850., 2500.); solver = DifferentialEquations.Tsit5())
julia> #= Plot the temperature in Celsius from 1850 to 2500 =#
julia> plot(sol; idxs = (EarthSystemClimateInterpretableModel.Temp_surface_C))
```

![image](https://github.com/JKRT/ESCIMO.jl/assets/8775827/57fcd86c-411d-4bc4-92be-68609acc5633)


## Collaboration & Contact
Please email me at the email located here [LiU-page](https://liu.se/en/employee/johti17)
