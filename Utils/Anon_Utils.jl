mutable struct AnonFunc <: Function
    func::Function
    formula::String

    function AnonFunc(formula::String)
        new(eval(Meta.parse(formula)), formula)
    end
end

function Base.print(io::IO, λ::AnonFunc) 
    println(io, "$(λ.func): $(λ.formula)")
end

function (λ::AnonFunc)(x)
    return λ.func(x)
end