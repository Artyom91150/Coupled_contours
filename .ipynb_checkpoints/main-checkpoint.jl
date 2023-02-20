## Importing libraries ##
using YAML: load_file, write_file
using Configurations
using Dates: now, format
using JLD: save
println("Libraries uploaded")

## Importing utils ##
include("Utils\\Computation_Utils.jl")
include("Utils\\CallBack_Utils.jl")
include("Utils\\YAML_Utils.jl")
include("Utils\\Graphics_Utils.jl")

println("Utils uploaded")


## Loading configuration file ##
args = Base.ARGS
config_dict = load_file(args[1])
println("config uploaded")

StartUpTime = now()

TaskName = config_dict["Task"]
SaveFolder = config_dict["OutputProperties"]["DataSettings"]["SavePathFolder"] * "\\" * TaskName * " " * format(StartUpTime, "yyyy_mm_dd HH_MM_SS")
SaveFormat = config_dict["OutputProperties"]["DataSettings"]["DataFileFormat"]
mkpath(SaveFolder)

## Save configuration file
YAML.write_file(SaveFolder * "\\" * "config.yaml", config_dict)


# Required fileds
module_name = include(config_dict["ODE_File_Name"])
ODE_func = getfield(module_name, Meta.parse(config_dict["ODE_Func_Name"]))(YAML_Utils.MakeParametersDict(config_dict), eval(Meta.parse(config_dict["Couple_func"])))
InitCond = config_dict["ODESolver_setings"]["InitCond"]
time_span = config_dict["ODESolver_setings"]["time_span"]
println("Required fileds procceed")

# Optional fields
alg = haskey(config_dict["ODESolver_setings"], "alg") ? config_dict["ODESolver_setings"]["alg"] : nothing
ODESolverKwargs = haskey(config_dict["ODESolver_setings"], "kwargs") ? config_dict["ODESolver_setings"]["kwargs"] : nothing
trans_time = haskey(config_dict["ODESolver_setings"], "trans_time") ? config_dict["ODESolver_setings"]["trans_time"] : nothing
callback = haskey(config_dict["ODESolver_setings"], "callback") ? config_dict["ODESolver_setings"]["callback"] : nothing
println("Optional fileds procceed")

## Solve ODE system ##
sol = Computation_Utils.SolveODE(ODE_func, InitCond, time_span;
                                 alg = alg,
                                 trans_time = trans_time,
                                 kwargs = config_dict["ODESolver_setings"]["kwargs"],
                                 callback_func = callback);
println("SolveODE succeed")

## Save results ##
save(SaveFolder * "\\" * TaskName * "Data" * "." * SaveFormat, "sol", sol)
println("Results saved")

## Plot results ##
Graphics_Utils.PlotRequestedGraphics(sol, config_dict, SaveFolder)



