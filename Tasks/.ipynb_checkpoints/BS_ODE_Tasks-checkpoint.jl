module BS_ODE_Tasks
    import PyCall: pyimport, PyVector
    import JLD2: save_object, load_object
    using YAML: load_file, write_file
    import Statistics: mean, std
    using LinearAlgebra: eigen
    using SciMLBase
    using OrdinaryDiffEq

    export SolveODE, plotActivationDiagram, BS_Sync, GetPeriod, plotPeriod, SolveTanODE
    
    include("..\\Utils\\Computation_Utils.jl")
    include("..\\Utils\\PeriodSync_Utils.jl")
    include("..\\Utils\\CallBack_Utils.jl")
    include("..\\Utils\\IO_Utils.jl")

    pushfirst!(PyVector(pyimport("sys")."path"), "")
    MU = pyimport("Utils.Matplotlib_Utils")

    # По хорошему, должен подтягиваться из конфига, но пока что увы
    ODE_module = include("..\\ODEs\\BS_Uni_System.jl")
    #ODE_module = include("..\\ODEs\\SL_System.jl")
    #Couple_func = x -> 1/(1 + exp(10*(cos(x))))
    Couple_func = x -> 1 - cos(x)
    #callback = CallBack_Utils.ParseAnonCondition(x -> sin(x[1] - 0.25))
    callback = nothing

    function SolveODE(Task_Input, Task_Output, Current_folder)
    ## Make subdictionaries ##
        ODESystem_setings = Task_Input["ODESystem_setings"]
        ODESolver_settings = Task_Input["ODESolver_settings"]

    ## Required fileds ##
        # Module with ODE system implementations
        #ODE_module = include(ODESystem_setings["ODE_file_path"])

        # ODE system function for solving
        ODE_system = getfield(ODE_module, Meta.parse(ODESystem_setings["ODE_func_name"]))
        #cf = Parse_str_expr(ODESystem_setings["Support_properties"]["Couple_func"])
        ODE_system = ODE_system(IO_Utils.Make_parameters_dict(ODESystem_setings["ODE_Parameters"]), Couple_func)

        # Initial conditions for ODE solving
        Init_cond = ODESolver_settings["init_cond"]

        # Time span for ODE solver
        time_span = ODESolver_settings["time_span"]
        println("### Required fileds procceed ###")

    ## Optional fields ##
        # ODE solver algorithm
        alg = haskey(ODESolver_settings, "alg") ? ODESolver_settings["alg"] : nothing
        alg = typeof(alg) == String ? eval(Meta.parse(alg)) : alg

        # ODE solver kwargs
        ODESolver_kwargs = haskey(ODESolver_settings, "kwargs") ? ODESolver_settings["kwargs"] : nothing
        ODESolver_kwargs = IO_Utils.Parse_kwargs(ODESolver_kwargs)

        # Transit time for ODE solver
        trans_time = haskey(ODESolver_settings, "trans_time") ? ODESolver_settings["trans_time"] : nothing

        # Callback (event) function for ODE solver
        #callback = haskey(ODESolver_settings, "callback") ? ODESolver_settings["callback"] : nothing
        #callback = typeof(callback) == String ? eval(Meta.parse(callback)) : callback
        println("### Optional fileds procceed ###")

    ## Solve ODE system ##
        sol = Computation_Utils.SolveODE(ODE_system, Init_cond, time_span;
                                        alg = alg,
                                        trans_time = trans_time,
                                        kwargs = ODESolver_kwargs,
                                        callback = callback);
        println("### SolveODE succeed ###")

    ## Save results ##
        save_object(Current_folder * "\\" * Task_Output["Data_file_name"], sol)
        println("### SolveODE results saved ###")
    end

    function plotActivationDiagram(Task_Input, Task_Output, Current_folder)
        sol = IO_Utils.Load_JLD2_file(Current_folder * "\\" *Task_Input["Data_file_name"])
        title = haskey(Task_Input, "title") ? Task_Input["title"] : nothing
        savePath = Current_folder === nothing ? nothing : Current_folder * "\\" * Task_Output
        varNames = haskey(Task_Input, "varNames") ? Task_Input["varNames"] : nothing

        MU.plotActivationDiagram(sol, varNames = varNames, savePath = savePath, title = title)
    end

    function BS_Sync(Task_Input, Task_Output, Current_folder)
        sol = IO_Utils.Load_JLD2_file(Current_folder * "\\" *Task_Input["Data_file_name"])
        savePath = Current_folder === nothing ? nothing : Current_folder * "\\" * Task_Output
        Sync_config = getSync(sol)

        if savePath !== nothing 
            write_file(Current_folder * "\\" * Task_Output, Sync_config)
        end
    end

    function GetPeriod(Task_Input, Task_Output, Current_folder)
        sol = IO_Utils.Load_JLD2_file(Current_folder * "\\" *Task_Input["Data_file_name"])
        savePath = Current_folder === nothing ? nothing : Current_folder * "\\" * Task_Output
        Period = PeriodSync_Utils.getPeriodTime(sol)

        if savePath !== nothing
            save_object(Current_folder * "\\" * Task_Output, Period)
        end
    end

    function plotPeriod(Task_Input, Task_Output, Current_folder)
        sol = IO_Utils.Load_JLD2_file(Current_folder * "\\" *Task_Input["Data_file_name"])
        savePath = Current_folder === nothing ? nothing : Current_folder * "\\" * Task_Output

        if savePath !== nothing
            PeriodSync_Utils.plotNorms(sol; savePath = savePath)
        end
    end

    function getSync(sol)
        delay = 0
        Stds = [[std(map(cos, scndSol) - map(cos, frstSol)) for scndSol in sol.y[4:6]] for frstSol in sol.y[1:3]]
        Stds_pi = [[std(map(cos, scndSol) - map(x -> cos(x + pi), frstSol)) for scndSol in sol.y[4:6]] for frstSol in sol.y[1:3]]
        if minimum(Stds[1]) < minimum(Stds_pi[1]) 
            delay = 0
        else 
            delay = Float64(pi)
            Stds = Stds_pi
        end
        Syncs = Dict()
        for i = 1 : 3
            nearestSync = argmin(Stds[i])
            Syncs["phi$(i)"] = Dict("psi" => nearestSync, "std" => Stds[i][nearestSync], "delay" => delay)
        end
        return Syncs
    end

    function SolveTanODE(Task_Input, Task_Output, Current_folder)
        ## Make subdictionaries ##
            ODESystem_setings = Task_Input["ODESystem_setings"]
            ODESolver_settings = Task_Input["ODESolver_settings"]
    
            # ODE system function for solving
            Syncs = load_file(Current_folder * "\\" * ODESystem_setings["Support_properties"]["Syncs"])

            ODE_system = getfield(ODE_module, Meta.parse(ODESystem_setings["ODE_func_name"]))
            ODE_system = ODE_system(IO_Utils.Make_parameters_dict(ODESystem_setings["ODE_Parameters"]), Couple_func, Syncs)
    
            # Tangential ODE system function for solving
            Tan_system = getfield(ODE_module, Meta.parse(ODESystem_setings["Tan_func_name"]))
            Tan_system = Tan_system(IO_Utils.Make_parameters_dict(ODESystem_setings["ODE_Parameters"]), Couple_func, Syncs)
    
            # Initial conditions for ODE solving
            u0 = ODESolver_settings["u0"]
            if typeof(u0) == String
                u0 = IO_Utils.Load_JLD2_file(Current_folder * "\\" * u0)
                u0 = [x[end] for x in u0.y[1 : 3]]
            end

            # Initial conditions for Tangential ODE solving
            Q0 = ODESolver_settings["Q0"]
    
            # Time span for ODE solver
            time_span = ODESolver_settings["time_span"]
            if typeof(time_span) == String
                time_span = load_object(Current_folder * "\\" * time_span)
            end
            println("Required fileds procceed")
    
        ## Optional fields ##
            # ODE solver algorithm
            alg = haskey(ODESolver_settings, "alg") ? ODESolver_settings["alg"] : nothing
            alg = typeof(alg) == String ? eval(Meta.parse(alg)) : alg
    
            # ODE solver kwargs
            ODESolver_kwargs = haskey(ODESolver_settings, "kwargs") ? ODESolver_settings["kwargs"] : nothing
            ODESolver_kwargs = IO_Utils.Parse_kwargs(ODESolver_kwargs)
    
            # Transit time for ODE solver
            trans_time = haskey(ODESolver_settings, "trans_time") ? ODESolver_settings["trans_time"] : nothing
    
            println("Optional fileds procceed")
    
        ## Solve ODE system ##
            x, W = Computation_Utils.SolveTanODE(ODE_system, Tan_system, u0, Q0, time_span;
                                                alg = alg,
                                                trans_time = trans_time,
                                                kwargs = ODESolver_kwargs);
            println("SolveODE succeed")

        ## Eigen values of W ##
            eig = eigen(W)
    
        ## Save results ##
            write_file(Current_folder * "\\" * Task_Output, Dict("x" => x, "W" => W, "Eig" => eig))
            println("Results saved")

        end

end