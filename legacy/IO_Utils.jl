module IO_Utils
    using Dates: now, format
    using JLD2: load_object
    
    include("Callback_Utils.jl")
    include("Computation_Utils.jl")

    function ParseConfig(config)
        ## Parsing parameters of ODE system
        for parameter in config["Parameters"]
            parameter["value"] = EvalStringExpr(parameter["value"])
        end

        ## Parsing kwargs of ODE solver
        config["ODESolver_setings"]["kwargs"] = ParseKwargs(config["ODESolver_setings"]["kwargs"])

        ## Parsing kwargs of graphics plots if they exist
        for i in 1 : length(config["OutputProperties"]["RequestedPlots"])
            if hasfield(config["OutputProperties"]["RequestedPlots"][i], "kwargs")
                config["OutputProperties"]["RequestedPlots"][i] = ParseKwargs(config["OutputProperties"]["RequestedPlots"])
            end
        end

        ## Parsing string couple function into anon function
        config["Couple_func"] = EvalStringExpr(config["Couple_func"])

        ## Parsing string callback function into anon function
        if hasfield(config["ODESolver_setings"], "callback")
            config["ODESolver_setings"]["callback"] = Callback_Utils.ParseAnonCondition(config["ODESolver_setings"]["callback"])
        end

        ## Parsing 
    end

    function Make_results_folder(Output_path_folder, Task_name)
        Folder_path = Output_path_folder * "\\" * Task_name * " " * format(now(), "yyyy_mm_dd HH_MM_SS")
        mkpath(Folder_path)
        return Folder_path
    end

    # Transform system parameters list from config. file into dictionary with pairs "value" => value like:
    # [Dict("name" => a, "value" => 123), Dict(...), ...] => Dict("a" => 123, ...)
    function Make_parameters_dict(Parameters_list::Vector)
        Parameters_dict = Dict()
        for p in Parameters_list
            get!(Parameters_dict, p["name"], Parse_str_expr(p["value"]))
        end
        return Parameters_dict
    end

    function Parse_str_expr(str)
        if(typeof(str) == String)
            str = eval(Meta.parse(str))
        end
        return str
    end

    function Parse_kwargs(kwargs::Dict)
        Parsed_kwargs = Dict()
        for key in keys(kwargs)
            Parsed_kwargs[Meta.parse(key)] = Parse_str_expr(kwargs[key])
        end
        return Parsed_kwargs
    end

    function Load_JLD2_file(file_path)
        #JLD2.rconvert(::Type{Computation_Utils.py_sol}, nt::NamedTuple) = Computation_Utils.py_sol(nt.t, nt.y, nt.t_events, nt.y_events)
        return load_object(file_path)
        #return load(file_path, "sol", typemap=Dict("Main.BS_ODE_Solving.Computation_Utils.py_sol" => Computation_Utils.py_sol))
    end
end