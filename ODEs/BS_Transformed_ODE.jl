module BS_Transformed_ODE

export BS_Log, BS_Log_2, BS_Log_β

export Phi2Z, Z2Phi

using LaTeXStrings

# Φ -> Z Transform
Phi2Z(phi::Float64) = 2log(exp(1), tan(phi / 2.0))
Phi2Z(phi::String) = L"2log(exp(tan(\frac{%$phi}{2})))"

# Z -> Φ Transform
Z2Phi(z::Float64) = 2atan(exp(z / 2.0))
Z2Phi(z::String) = L"2atan(exp(\frac{%$z}{2}))"



### Transformed double BS system with "NECHESTNAYA" couple < Sin(ϕ - ψ)/2 >
@inbounds function BS_Log(u::T, p, t=0.0) where T
    x1, x2, x3, y1, y2, y3 = u
    A, B, C, Eps = p

    th1, th2, th3, th4, th5, th6 = tanh.(u)
    sh1, sh2, sh3, sh4, sh5, sh6 = sinh.(u)

    dx1 = (A * (th2 - th3) - B * th1 + C) - Eps * (sech(y1) * sh1 - th4)
    dx2 = (A * (th3 - th1) - B * th2 + C) - Eps * (sech(y2) * sh2 - th5) 
    dx3 = (A * (th1 - th2) - B * th3 + C) - Eps * (sech(y3) * sh3 - th6)
    dx4 = (-A * (th5 - th6) - B * th4 + C) - Eps * (sech(x1) * sh4 - th1)
    dx5 = (-A * (th6 - th4) - B * th5 + C) - Eps * (sech(x2) * sh5 - th2)
    dx6 = (-A * (th4 - th5) - B * th6 + C) - Eps * (sech(x3) * sh6 - th3)

    return T([dx1, dx2, dx3, dx4, dx5, dx6])
end

### The same system but with x/2
@inbounds function BS_Log_2(u::T, p, t=0.0) where T
    u = u ./ 2.0
    x1, x2, x3, y1, y2, y3 = u
    A, B, C, Eps = p

    th1, th2, th3, th4, th5, th6 = tanh.(u)
    sh1, sh2, sh3, sh4, sh5, sh6 = sinh.(u)

    dx1 = 2.0 * (A * (th2 - th3) - B * th1 + C) - Eps * (sech(y1) * sh1 - th4)
    dx2 = 2.0 * (A * (th3 - th1) - B * th2 + C) - Eps * (sech(y2) * sh2 - th5) 
    dx3 = 2.0 * (A * (th1 - th2) - B * th3 + C) - Eps * (sech(y3) * sh3 - th6)
    dx4 = 2.0 * (-A * (th5 - th6) - B * th4 + C) - Eps * (sech(x1) * sh4 - th1)
    dx5 = 2.0 * (-A * (th6 - th4) - B * th5 + C) - Eps * (sech(x2) * sh5 - th2)
    dx6 = 2.0 * (-A * (th4 - th5) - B * th6 + C) - Eps * (sech(x3) * sh6 - th3)

    return T([dx1, dx2, dx3, dx4, dx5, dx6])
end


### Transformed double BS system with "CHESTNAYA" couple < Cos(ψ + β) >
@inbounds function BS_Log_β(u::T, p, t=0.0) where T
    phi1, phi2, phi3, psi1, psi2, psi3 = u
    A, B, C, Eps, β = p

    dx1 = sin(phi1) * (A * (cos(phi3) - cos(phi2)) + B * cos(phi1) + C - Eps * cos(β + psi1))
    dx2 = sin(phi2) * (A * (cos(phi1) - cos(phi3)) + B * cos(phi2) + C - Eps * cos(β + psi2))
    dx3 = sin(phi3) * (A * (cos(phi2) - cos(phi1)) + B * cos(phi3) + C - Eps * cos(β + psi3))
    dx4 = sin(psi1) * (-A * (cos(psi3) - cos(psi2)) + B * cos(psi1) + C - Eps * cos(β + phi1))
    dx5 = sin(psi2) * (-A * (cos(psi1) - cos(psi3)) + B * cos(psi2) + C - Eps * cos(β + phi2))
    dx6 = sin(psi3) * (-A * (cos(psi2) - cos(psi1)) + B * cos(psi3) + C - Eps * cos(β + phi3))

    return T([dx1, dx2, dx3, dx4, dx5, dx6])
end

end