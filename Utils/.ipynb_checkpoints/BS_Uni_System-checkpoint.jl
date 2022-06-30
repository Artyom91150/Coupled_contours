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


function single_Bick(X, p, t)
    K, a2, a4, r = p
    dX = zeros(3)
    for i = 1:3
        dX[i] = 2 * g2_2(X[i], a2, r) - (K/2) * (g4_2( X[i] + X[(i + 1) % 3 + 1], a4) + g4_2( X[i] - X[(i + 1) % 3 + 1], a4)) + (K/2) * (g4_2( X[i] + X[(i) % 3 + 1], a4) + g4_2( X[i] - X[(i) % 3 + 1], a4))
    end
    return dX
end

function couple_sin(x)
    return sin(x)
end

function couple_cos(x)
    return 1 - cos(x)
end

function couple_exp(x)
    return 1/(1 + exp(50*cos(x)))
end

function BS_Uni_Sys(couple_func = couple_cos)
    function rhs(dX, X, p, t)
        K, a2, a4, r, Eps = p
        #dX = zeros(6)
        dX[1:3] = single_Bick(X[1:3], [K, a2, a4, r], t) + Eps * [couple_func(p) for p in  X[4:6] - X[1:3]]
        dX[4:6] = single_Bick(X[4:6], [-K, a2, a4, r], t) + Eps * [couple_func(p) for p in  X[1:3] - X[4:6]]
        return SVector{6}(dX)
    end
return rhs
end