struct MySavedValues{tType, savevalType}
    t::Vector{tType}
    saveval::Vector{savevalType}
end

function MySavedValues(::Type{tType}, ::Type{savevalType}) where {tType,savevalType}
    MySavedValues{tType, savevalType}(Vector{tType}(), Vector{savevalType}())
end

function Base.show(io::IO, saved_values::MySavedValues)
    tType = eltype(saved_values.t)
    savevalType = eltype(saved_values.saveval)
    print(io, "SavedValues{tType=", tType, ", savevalType=", savevalType, "}",
                "\nt:\n", saved_values.t, "\nsaveval:\n", saved_values.saveval)
end

mutable struct MySavingAffect{tType, savevalType}
    saved_values::MySavedValues{tType, savevalType}
end

function (affect!::MySavingAffect)(integrator)
    push!(affect!.saved_values.t, integrator.t)
    push!(affect!.saved_values.saveval, Tuple(Float64(v) for v in integrator.u))
end

function GetCallBack(Condition)
    observer_array = MySavedValues(Float64, Tuple{Float64, Float64, Float64, Float64, Float64, Float64})
    affectThis! = MySavingAffect(observer_array)
    return [ContinuousCallback(Condition, affectThis!, nothing), observer_array]
end