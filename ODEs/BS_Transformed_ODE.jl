module BS_Transformed_ODE

export Get_Fast_BS_F, Get_Fast_BS_R, Get_Fast_BS_R1

export Phi2Z, Z2Phi

using LaTeXStrings

# Φ -> Z Transform
Phi2Z(phi::Float64) = 2log(exp(1), tan(phi / 2.0))
Phi2Z(phi::String) = L"2log(exp(tan(\frac{%$phi}{2})))"

# Z -> Φ Transform
Z2Phi(z::Float64) = 2atan(exp(z / 2.0))
Z2Phi(z::String) = L"2atan(exp(\frac{%$z}{2}))"





### Transformed double BS system with "NECHESTNAYA" couple < Sin(ϕ - ψ)/2 >
function Get_Fast_BS_F()

    function Unit_Rhs(cur, prev, next, A, B, C)
        return 2.0 * (A * (tanh(next / 2.0) - tanh(prev / 2.0)) - B * tanh(cur / 2.0) + C)
    end

    function Population_Rhs(X, A, B, C)
        x1, x2, x3 = X
        return [Unit_Rhs(x1, x3, x2, A, B, C), 
                Unit_Rhs(x2, x1, x3, A, B, C),
                Unit_Rhs(x3, x2, x1, A, B, C)]
    end

    function Couple(x, y) 
        return sech(y / 2.0) * sinh(x / 2.0) - tanh(y / 2.0)
    end

    @inbounds function BS_ODE_Duo_Fast(X::T, p, t = 0.0) where T
        z1, z2, z3, w1, w2, w3 = X
        A, B, C, Eps = p
    
        dZ = Population_Rhs(X[1:3], A, B, C)
        dZ .-= Eps .* Couple.(X[1:3], X[4:6])

        dW = Population_Rhs(X[4:6], -A, B, C)
        dW .-= Eps .* Couple.(X[4:6], X[1:3])

        return T([dZ; dW])
    end
    return BS_ODE_Duo_Fast
end


function f(u::T, p ,t = 0.0) where T
    x1, x2, x3, y1, y2, y3 = u
    A, B, C, Eps = p

    dx1 = 2.0 * (A * (tanh(x2 / 2.0) - tanh(x3 / 2.0)) - B * tanh(x1 / 2.0) + C) - Eps * (sech(y1 / 2.0) * sinh(x1 / 2.0) - tanh(y1 / 2.0))
    dx2 = 2.0 * (A * (tanh(x3 / 2.0) - tanh(x1 / 2.0)) - B * tanh(x2 / 2.0) + C) - Eps * (sech(y2 / 2.0) * sinh(x2 / 2.0) - tanh(y2 / 2.0)) 
    dx3 = 2.0 * (A * (tanh(x1 / 2.0) - tanh(x2 / 2.0)) - B * tanh(x3 / 2.0) + C) - Eps * (sech(y3 / 2.0) * sinh(x3 / 2.0) - tanh(y3 / 2.0))
    dy1 = 2.0 * (-A * (tanh(y2 / 2.0) - tanh(y3 / 2.0)) - B * tanh(y1 / 2.0) + C) - Eps * (sech(x1 / 2.0) * sinh(y1 / 2.0) - tanh(x1 / 2.0))
    dy2 = 2.0 * (-A * (tanh(y3 / 2.0) - tanh(y1 / 2.0)) - B * tanh(y2 / 2.0) + C) - Eps * (sech(x2 / 2.0) * sinh(y2 / 2.0) - tanh(x2 / 2.0))
    dy3 = 2.0 * (-A * (tanh(y1 / 2.0) - tanh(y2 / 2.0)) - B * tanh(y3 / 2.0) + C) - Eps * (sech(x3 / 2.0) * sinh(y3 / 2.0) - tanh(x3 / 2.0))

    return T([dx1, dx2, dx3, dy1, dy2, dy3])
end




### Transformed double BS system with "CHESTNAYA" couple < Cos(ψ + β) >
function Get_Fast_BS_R()

    function Unit_Rhs(cur, prev, next, A, B, C)
        return 2.0 * (A * (tanh(next / 2.0) - tanh(prev / 2.0)) - B * tanh(cur / 2.0) + C)
    end

    function Population_Rhs(X, A, B, C)
        x1, x2, x3 = X
        return [Unit_Rhs(x1, x3, x2, A, B, C), 
                Unit_Rhs(x2, x1, x3, A, B, C),
                Unit_Rhs(x3, x2, x1, A, B, C)]
    end

    function Couple(x, β) 
        return 2.0 * (sin(β) / cosh(x / 2.0) + cos(β) * tanh(x / 2.0))
    end

    @inbounds function BS_ODE_Duo_Fast(X::T, p, t = 0.0) where T
        z1, z2, z3, w1, w2, w3 = X
        A, B, C, Eps, β = p
    
        dZ = Population_Rhs(X[1:3], A, B, C)
        dZ .+= Eps .* Couple.(X[4:6], β)

        dW = Population_Rhs(X[4:6], -A, B, C)
        dW .+= Eps .* Couple.(X[1:3], β)

        return T([dZ; dW])
    end
    return BS_ODE_Duo_Fast
end


function Get_Fast_BS_R1()
    @inbounds function BS_ODE_Duo_Fast(X::T, p, t = 0.0) where T
        phi1, phi2, phi3, psi1, psi2, psi3 = X
        A, B, C, Eps, β = p
    
        dPhi1 = sin(phi1) * (A * (cos(phi3) - cos(phi2)) + B * cos(phi1) + C - Eps * cos(β + psi1))
        dPhi2 = sin(phi2) * (A * (cos(phi1) - cos(phi3)) + B * cos(phi2) + C - Eps * cos(β + psi2))
        dPhi3 = sin(phi3) * (A * (cos(phi2) - cos(phi1)) + B * cos(phi3) + C - Eps * cos(β + psi3))
    
        dPsi1 = sin(psi1) * (-A * (cos(psi3) - cos(psi2)) + B * cos(psi1) + C - Eps * cos(β + phi1))
        dPsi2 = sin(psi2) * (-A * (cos(psi1) - cos(psi3)) + B * cos(psi2) + C - Eps * cos(β + phi2))
        dPsi3 = sin(psi3) * (-A * (cos(psi2) - cos(psi1)) + B * cos(psi3) + C - Eps * cos(β + phi3))
        return T([dPhi1, dPhi2, dPhi3, dPsi1, dPsi2, dPsi3])
    end
end

end