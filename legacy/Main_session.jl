using YAML: load_file, write_file
using Configurations

include("Utils\\IO_Utils.jl")

## Loading configuration file ##
args = Base.ARGS
Session_config = load_file(args[1])
println("### Config uploaded ###")

## Loading files with task module ##
for module_path in Session_config["Session_Input"]["Reuired_modules"]
    include(module_path)
    module_name = module_path[findlast("\\", module_path)[1] + 1 : findlast(".", module_path)[1] - 1]
    eval(Meta.parse("using ." * module_name))
end
println("### Task modules has been loaded ###")

## Making current session folder ##
Folder_path = IO_Utils.Make_results_folder(Session_config["Session_folder"], Session_config["Session_name"])
println("### Folder has been created ###")

## Saving config file ##
write_file(Folder_path * "\\" * "Session_config.yaml", Session_config)

## Tasks execution ##
for Task_config in Session_config["Session_Input"]["Tasks_List"]

    ## Importing task function ##
    Task_Func = eval(Meta.parse(Task_config["Task_Func"]))

    ## Execution task function ##
    Task_Func(Task_config["Task_Input"], Task_config["Task_Output"], Folder_path)
    println("## $(Task_Func) completed ##")
end