function BS_Sin_Sys(X, p, t)
    K, a2, a4, r, Eps = p
    x1, x2, x3, y1, y2, y3 = X
    dX = zeros(6)
    dX[1] = sin(a2 - x1) - sin(a2 + x1) - (K*(sin(a4 + x1 + x2)/2 + sin(a4 + x1 - x2)/2 - sin(a4 - x1 + x2)/2 + sin(x1 - a4 + x2)/2))/2 + (K*(sin(a4 + x1 + x3)/2 + sin(a4 + x1 - x3)/2 - sin(a4 - x1 + x3)/2 + sin(x1 - a4 + x3)/2))/2 - r*sin(2*a2 - 2*x1) + r*sin(2*a2 + 2*x1) - Eps*sin(x1 - y1)
    dX[2] = sin(a2 - x2) - sin(a2 + x2) + (K*(sin(a4 + x1 + x2)/2 - sin(a4 + x1 - x2)/2 + sin(a4 - x1 + x2)/2 + sin(x1 - a4 + x2)/2))/2 - (K*(sin(a4 + x2 + x3)/2 + sin(a4 + x2 - x3)/2 - sin(a4 - x2 + x3)/2 + sin(x2 - a4 + x3)/2))/2 - r*sin(2*a2 - 2*x2) + r*sin(2*a2 + 2*x2) - Eps*sin(x2 - y2)
    dX[3] = sin(a2 - x3) - sin(a2 + x3) - (K*(sin(a4 + x1 + x3)/2 - sin(a4 + x1 - x3)/2 + sin(a4 - x1 + x3)/2 + sin(x1 - a4 + x3)/2))/2 + (K*(sin(a4 + x2 + x3)/2 - sin(a4 + x2 - x3)/2 + sin(a4 - x2 + x3)/2 + sin(x2 - a4 + x3)/2))/2 - r*sin(2*a2 - 2*x3) + r*sin(2*a2 + 2*x3) - Eps*sin(x3 - y3)
    dX[4] = sin(a2 - y1) - sin(a2 + y1) - r*sin(2*a2 - 2*y1) + r*sin(2*a2 + 2*y1) + (K*(sin(a4 + y1 + y2)/2 + sin(a4 + y1 - y2)/2 - sin(a4 - y1 + y2)/2 + sin(y1 - a4 + y2)/2))/2 - (K*(sin(a4 + y1 + y3)/2 + sin(a4 + y1 - y3)/2 - sin(a4 - y1 + y3)/2 + sin(y1 - a4 + y3)/2))/2 + Eps*sin(x1 - y1)
    dX[5] = sin(a2 - y2) - sin(a2 + y2) - r*sin(2*a2 - 2*y2) + r*sin(2*a2 + 2*y2) - (K*(sin(a4 + y1 + y2)/2 - sin(a4 + y1 - y2)/2 + sin(a4 - y1 + y2)/2 + sin(y1 - a4 + y2)/2))/2 + (K*(sin(a4 + y2 + y3)/2 + sin(a4 + y2 - y3)/2 - sin(a4 - y2 + y3)/2 + sin(y2 - a4 + y3)/2))/2 + Eps*sin(x2 - y2)
    dX[6] = sin(a2 - y3) - sin(a2 + y3) - r*sin(2*a2 - 2*y3) + r*sin(2*a2 + 2*y3) + (K*(sin(a4 + y1 + y3)/2 - sin(a4 + y1 - y3)/2 + sin(a4 - y1 + y3)/2 + sin(y1 - a4 + y3)/2))/2 - (K*(sin(a4 + y2 + y3)/2 - sin(a4 + y2 - y3)/2 + sin(a4 - y2 + y3)/2 + sin(y2 - a4 + y3)/2))/2 + Eps*sin(x3 - y3)
    return SVector{6}(dX)
end

function BS_Sin_Jac(X, p, t)
    x1, x2, x3, y1, y2, y3 = X
    K, a2, a4, r, Eps = p
    Jac = zeros(Float64, 6, 6)
    Jac[1, :] = [2*r*cos(2*a2 - 2*x1) - cos(a2 - x1) - cos(a2 + x1) + 2*r*cos(2*a2 + 2*x1) + Eps*cos(x1 - y1) - (K*(cos(a4 + x1 + x2)/2 + cos(a4 + x1 - x2)/2 + cos(a4 - x1 + x2)/2 + cos(x1 - a4 + x2)/2))/2 + (K*(cos(a4 + x1 + x3)/2 + cos(a4 + x1 - x3)/2 + cos(a4 - x1 + x3)/2 + cos(x1 - a4 + x3)/2))/2, -(K*(cos(a4 + x1 + x2)/2 - cos(a4 + x1 - x2)/2 - cos(a4 - x1 + x2)/2 + cos(x1 - a4 + x2)/2))/2, (K*(cos(a4 + x1 + x3)/2 - cos(a4 + x1 - x3)/2 - cos(a4 - x1 + x3)/2 + cos(x1 - a4 + x3)/2))/2, -Eps*cos(x1 - y1), 0, 0]
    Jac[2, :] = [(K*(cos(a4 + x1 + x2)/2 - cos(a4 + x1 - x2)/2 - cos(a4 - x1 + x2)/2 + cos(x1 - a4 + x2)/2))/2, 2*r*cos(2*a2 - 2*x2) - cos(a2 - x2) - cos(a2 + x2) + 2*r*cos(2*a2 + 2*x2) + Eps*cos(x2 - y2) + (K*(cos(a4 + x1 + x2)/2 + cos(a4 + x1 - x2)/2 + cos(a4 - x1 + x2)/2 + cos(x1 - a4 + x2)/2))/2 - (K*(cos(a4 + x2 + x3)/2 + cos(a4 + x2 - x3)/2 + cos(a4 - x2 + x3)/2 + cos(x2 - a4 + x3)/2))/2, -(K*(cos(a4 + x2 + x3)/2 - cos(a4 + x2 - x3)/2 - cos(a4 - x2 + x3)/2 + cos(x2 - a4 + x3)/2))/2, 0, -Eps*cos(x2 - y2), 0]
    Jac[3, :] = [-(K*(cos(a4 + x1 + x3)/2 - cos(a4 + x1 - x3)/2 - cos(a4 - x1 + x3)/2 + cos(x1 - a4 + x3)/2))/2, (K*(cos(a4 + x2 + x3)/2 - cos(a4 + x2 - x3)/2 - cos(a4 - x2 + x3)/2 + cos(x2 - a4 + x3)/2))/2, 2*r*cos(2*a2 - 2*x3) - cos(a2 - x3) - cos(a2 + x3) + 2*r*cos(2*a2 + 2*x3) + Eps*cos(x3 - y3) - (K*(cos(a4 + x1 + x3)/2 + cos(a4 + x1 - x3)/2 + cos(a4 - x1 + x3)/2 + cos(x1 - a4 + x3)/2))/2 + (K*(cos(a4 + x2 + x3)/2 + cos(a4 + x2 - x3)/2 + cos(a4 - x2 + x3)/2 + cos(x2 - a4 + x3)/2))/2, 0, 0, -Eps*cos(x3 - y3)]
    Jac[4, :] = [-Eps*cos(x1 - y1), 0, 0, 2*r*cos(2*a2 - 2*y1) - cos(a2 - y1) - cos(a2 + y1) + 2*r*cos(2*a2 + 2*y1) + Eps*cos(x1 - y1) + (K*(cos(a4 + y1 + y2)/2 + cos(a4 + y1 - y2)/2 + cos(a4 - y1 + y2)/2 + cos(y1 - a4 + y2)/2))/2 - (K*(cos(a4 + y1 + y3)/2 + cos(a4 + y1 - y3)/2 + cos(a4 - y1 + y3)/2 + cos(y1 - a4 + y3)/2))/2, (K*(cos(a4 + y1 + y2)/2 - cos(a4 + y1 - y2)/2 - cos(a4 - y1 + y2)/2 + cos(y1 - a4 + y2)/2))/2, -(K*(cos(a4 + y1 + y3)/2 - cos(a4 + y1 - y3)/2 - cos(a4 - y1 + y3)/2 + cos(y1 - a4 + y3)/2))/2]
    Jac[5, :] = [0, -Eps*cos(x2 - y2), 0, -(K*(cos(a4 + y1 + y2)/2 - cos(a4 + y1 - y2)/2 - cos(a4 - y1 + y2)/2 + cos(y1 - a4 + y2)/2))/2, 2*r*cos(2*a2 - 2*y2) - cos(a2 - y2) - cos(a2 + y2) + 2*r*cos(2*a2 + 2*y2) + Eps*cos(x2 - y2) - (K*(cos(a4 + y1 + y2)/2 + cos(a4 + y1 - y2)/2 + cos(a4 - y1 + y2)/2 + cos(y1 - a4 + y2)/2))/2 + (K*(cos(a4 + y2 + y3)/2 + cos(a4 + y2 - y3)/2 + cos(a4 - y2 + y3)/2 + cos(y2 - a4 + y3)/2))/2, (K*(cos(a4 + y2 + y3)/2 - cos(a4 + y2 - y3)/2 - cos(a4 - y2 + y3)/2 + cos(y2 - a4 + y3)/2))/2]
    Jac[6, :] = [0, 0, -Eps*cos(x3 - y3), (K*(cos(a4 + y1 + y3)/2 - cos(a4 + y1 - y3)/2 - cos(a4 - y1 + y3)/2 + cos(y1 - a4 + y3)/2))/2, -(K*(cos(a4 + y2 + y3)/2 - cos(a4 + y2 - y3)/2 - cos(a4 - y2 + y3)/2 + cos(y2 - a4 + y3)/2))/2, 2*r*cos(2*a2 - 2*y3) - cos(a2 - y3) - cos(a2 + y3) + 2*r*cos(2*a2 + 2*y3) + Eps*cos(x3 - y3) + (K*(cos(a4 + y1 + y3)/2 + cos(a4 + y1 - y3)/2 + cos(a4 - y1 + y3)/2 + cos(y1 - a4 + y3)/2))/2 - (K*(cos(a4 + y2 + y3)/2 + cos(a4 + y2 - y3)/2 + cos(a4 - y2 + y3)/2 + cos(y2 - a4 + y3)/2))/2]
    return Jac
end