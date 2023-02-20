module BS_Uni_Sys
using StaticArrays
using DifferentialEquations
using DynamicalSystems

    function g2_1(x, a2, r)
        return sin(x + a2) - r * sin(2 * (x + a2))
    end

    function g2_2(x, a2, r)
        return (g2_1(-x, a2, r) - g2_1(x, a2, r))/2
    end

    function g4_1(x, a4)
        return sin(x + a4)
    end

    function g4_2(x, a4)
        return (g4_1(-x, a4) - g4_1(x, a4))/2
    end

    function Get_Sngl_Sys(dict::Dict)
        a2 = dict["a2"]
        a4 = dict["a4"]
        K = dict["K"]
        r = dict["r"]
        function BS_Sngl_Sys!(dX, X, p, t)
            #K, a2, a4, r = p
            dX[1] = 2 * g2_2(X[1], a2, r) - (K/2) * (g4_2( X[1] + X[3], a4) + g4_2( X[1] - X[3], a4)) + (K/2) * (g4_2( X[1] + X[2], a4) + g4_2( X[1] - X[2], a4))
            dX[2] = 2 * g2_2(X[2], a2, r) - (K/2) * (g4_2( X[2] + X[1], a4) + g4_2( X[2] - X[1], a4)) + (K/2) * (g4_2( X[2] + X[3], a4) + g4_2( X[2] - X[3], a4))
            dX[3] = 2 * g2_2(X[3], a2, r) - (K/2) * (g4_2( X[3] + X[2], a4) + g4_2( X[3] - X[2], a4)) + (K/2) * (g4_2( X[3] + X[1], a4) + g4_2( X[3] - X[1], a4))
            return SVector{3}(dX)
        end
        return BS_Sngl_Sys!
    end

    function Get_Uni_Sys(dict::Dict, Couple_Func)
        # System parameters from config. file
        a2 = dict["a2"]
        a4 = dict["a4"]
        K = dict["K"]
        r = dict["r"]
        Eps = dict["Eps"]

        Sngl_Sys_1! = Get_Sngl_Sys(Dict("a2" => a2, "a4" => a4, "K" => K, "r" => r))
        Sngl_Sys_2! = Get_Sngl_Sys(Dict("a2" => a2, "a4" => a4, "K" => -K, "r" => r))

        function Uni_Sys!(dX, X, p, t)
            #K, a2, a4, r, Eps = p
            dX[1:3] = Sngl_Sys_1!(dX[1:3], X[1:3], [], t) + Eps * [Couple_Func(p) for p in  X[4:6] - X[1:3]]
            dX[4:6] = Sngl_Sys_2!(dX[4:6], X[4:6], [], t) + Eps * [Couple_Func(p) for p in  X[1:3] - X[4:6]]
            return SVector{6}(dX)
        end
    return Uni_Sys!
    end

    function Get_Reduced_Sys(dict::Dict, Couple_Func, syncWith)
        a2 = dict["a2"]
        a4 = dict["a4"]
        K = dict["K"]
        r = dict["r"]
        Eps = dict["Eps"]
        function BS_Red_Sys!(dX, X, p, t)
            psi1 = X[syncWith["phi1"]["psi"]] + syncWith["phi1"]["delay"]
            psi2 = X[syncWith["phi2"]["psi"]] + syncWith["phi2"]["delay"]
            psi3 = X[syncWith["phi3"]["psi"]] + syncWith["phi3"]["delay"]

            dX[1] = 2 * g2_2(X[1], a2, r) - (K/2) * (g4_2( X[1] + X[3], a4) + g4_2( X[1] - X[3], a4)) + (K/2) * (g4_2( X[1] + X[2], a4) + g4_2( X[1] - X[2], a4)) + Eps * Couple_Func(X[1] - psi1)
            dX[2] = 2 * g2_2(X[2], a2, r) - (K/2) * (g4_2( X[2] + X[1], a4) + g4_2( X[2] - X[1], a4)) + (K/2) * (g4_2( X[2] + X[3], a4) + g4_2( X[2] - X[3], a4)) + Eps * Couple_Func(X[2] - psi2)
            dX[3] = 2 * g2_2(X[3], a2, r) - (K/2) * (g4_2( X[3] + X[2], a4) + g4_2( X[3] - X[2], a4)) + (K/2) * (g4_2( X[3] + X[1], a4) + g4_2( X[3] - X[1], a4)) + Eps * Couple_Func(X[3] - psi3)
            return SVector{3}(dX)
        end
        return BS_Red_Sys!
    end

    function Get_Tan_Red_Sys(dict::Dict, Couple_Func, syncWith)
        function BS_Tan_Red_Sys!(dW, W, X)
            phi1, phi2, phi3 = X
            psi1 = X[syncWith["phi1"]["psi"]] + syncWith["phi1"]["delay"]
            psi2 = X[syncWith["phi2"]["psi"]] + syncWith["phi2"]["delay"]
            psi3 = X[syncWith["phi3"]["psi"]] + syncWith["phi3"]["delay"]
            dW .= jacobian(ContinuousDynamicalSystem(Get_Uni_Sys(dict, Couple_Func), zeros(6), []), [phi1, phi2, phi3, psi1, psi2, psi3]) * W
            return SVector{6}(dW)
        end
        return BS_Tan_Red_Sys!
    end

    function Get_Tan_Full_Sys(dict::Dict, Couple_Func)
        function BS_Tan_Full_Sys!(dW, W, X)
            dW .= jacobian(ContinuousDynamicalSystem(Get_Uni_Sys(dict, Couple_Func), zeros(6), []), X) * W
            return SVector{6}(dW)
        end
        return BS_Tan_Full_Sys!
    end
end