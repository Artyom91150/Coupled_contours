mutable struct Session
    StartUp_Time::DateTime
    Folder_path::String
    Log_Name::String
    Tasks::Vector{AbstractTask}
    Results::Vector{Any}

    function Session(SessionName::String; Result_Folder_Name::String = "Results", Log_Name::String = "OutputLog.txt")
        StartUp_Time = now()
        Folder_path = "$(Result_Folder_Name)//$(SessionName) [$(Dates.format(StartUp_Time, "yyyy_mm_dd HH_MM_SS"))]"
        mkpath(Folder_path)

        Tasks = Vector{AbstractTask}(undef, 0)
        Results = Vector{Any}(undef, 0)
        
        WriteLog(Folder_path, Log_Name, "Session <$(Folder_path)> output log. ")

        new(StartUp_Time, Folder_path, Log_Name, Tasks, Results)
    end
end

function WriteLog(Folder_path::String, Log_Name::String, Message::String)
    Log_path = "$(Folder_path)//$(Log_Name)"
    open(Log_path, "w") do log
        write(log, Message)
        write(log, "\n\n")
    end
end

function PushTask!(s::Session, task::T) where T <: AbstractTask
    push!(s.Tasks, task)
    println("Task $T successfully added")
end

function ClearTasks!(s::Session)
    n_tasks = length(s.Tasks)
    s.Tasks = Vector{AbstractTask}(undef, 0)
    println("$(n_tasks) tasks successfully cleared")
end

function Base.show(io::IO, s::Session)
    for property in propertynames(s)
        value = getproperty(s, property)
        println(io, "[$property::$(typeof(value))]: $(value)")
    end
end

function Base.show(io::IO, ts::Vector{AbstractTask})
    for (i, t) in enumerate(ts)
        println(io, "Task #$(i) ::$(typeof(t))")
        println(io, t)
    end
end