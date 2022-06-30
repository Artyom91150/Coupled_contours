function BS_Exp_Sys(dX, X, p, t)
    K, a2, a4, r, Eps = p
    x1, x2, x3, y1, y2, y3 = X
    #dX = zeros(6)
    dX[1] = sin(a2 - x1) - sin(a2 + x1) + Eps/(0.5 + exp(5*cos(x1 - y1))) - (K*(sin(a4 + x1 + x2)/2 + sin(a4 + x1 - x2)/2 - sin(a4 - x1 + x2)/2 + sin(x1 - a4 + x2)/2))/2 + (K*(sin(a4 + x1 + x3)/2 + sin(a4 + x1 - x3)/2 - sin(a4 - x1 + x3)/2 + sin(x1 - a4 + x3)/2))/2 - r*sin(2*a2 - 2*x1) + r*sin(2*a2 + 2*x1)
    dX[2] = sin(a2 - x2) - sin(a2 + x2) + Eps/(0.5 + exp(5*cos(x2 - y2))) + (K*(sin(a4 + x1 + x2)/2 - sin(a4 + x1 - x2)/2 + sin(a4 - x1 + x2)/2 + sin(x1 - a4 + x2)/2))/2 - (K*(sin(a4 + x2 + x3)/2 + sin(a4 + x2 - x3)/2 - sin(a4 - x2 + x3)/2 + sin(x2 - a4 + x3)/2))/2 - r*sin(2*a2 - 2*x2) + r*sin(2*a2 + 2*x2)
    dX[3] = sin(a2 - x3) - sin(a2 + x3) + Eps/(0.5 + exp(5*cos(x3 - y3))) - (K*(sin(a4 + x1 + x3)/2 - sin(a4 + x1 - x3)/2 + sin(a4 - x1 + x3)/2 + sin(x1 - a4 + x3)/2))/2 + (K*(sin(a4 + x2 + x3)/2 - sin(a4 + x2 - x3)/2 + sin(a4 - x2 + x3)/2 + sin(x2 - a4 + x3)/2))/2 - r*sin(2*a2 - 2*x3) + r*sin(2*a2 + 2*x3)
    dX[4] = sin(a2 - y1) - sin(a2 + y1) + Eps/(0.5 + exp(5*cos(x1 - y1))) - r*sin(2*a2 - 2*y1) + r*sin(2*a2 + 2*y1) + (K*(sin(a4 + y1 + y2)/2 + sin(a4 + y1 - y2)/2 - sin(a4 - y1 + y2)/2 + sin(y1 - a4 + y2)/2))/2 - (K*(sin(a4 + y1 + y3)/2 + sin(a4 + y1 - y3)/2 - sin(a4 - y1 + y3)/2 + sin(y1 - a4 + y3)/2))/2
    dX[5] = sin(a2 - y2) - sin(a2 + y2) + Eps/(0.5 + exp(5*cos(x2 - y2))) - r*sin(2*a2 - 2*y2) + r*sin(2*a2 + 2*y2) - (K*(sin(a4 + y1 + y2)/2 - sin(a4 + y1 - y2)/2 + sin(a4 - y1 + y2)/2 + sin(y1 - a4 + y2)/2))/2 + (K*(sin(a4 + y2 + y3)/2 + sin(a4 + y2 - y3)/2 - sin(a4 - y2 + y3)/2 + sin(y2 - a4 + y3)/2))/2
    dX[6] = sin(a2 - y3) - sin(a2 + y3) + Eps/(0.5 + exp(5*cos(x3 - y3))) - r*sin(2*a2 - 2*y3) + r*sin(2*a2 + 2*y3) + (K*(sin(a4 + y1 + y3)/2 - sin(a4 + y1 - y3)/2 + sin(a4 - y1 + y3)/2 + sin(y1 - a4 + y3)/2))/2 - (K*(sin(a4 + y2 + y3)/2 - sin(a4 + y2 - y3)/2 + sin(a4 - y2 + y3)/2 + sin(y2 - a4 + y3)/2))/2
    return SVector{6}(dX)
end


