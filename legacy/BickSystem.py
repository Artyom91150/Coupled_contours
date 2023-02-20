import math as np

def singeBickSystem(K, a2, a4, r):
    def rhs(t, X):
        dXdt = [None]*3
        for i in range(3):
            dXdt[i] = 2 * g2Func2(X[i], a2, r) - K/2 * (g4Func2(X[(i + 2) % 3] + X[i], a4) + g4Func2(X[i] - X[(i + 2) % 3], a4)) + K/2 * (g4Func2(X[(i + 1) % 3] + X[i], a4) + g4Func2(X[i] - X[(i + 1) % 3], a4))
        return dXdt
    return rhs

def g4Func1(x, a4):
    g4 = np.sin(x + a4)
    return g4

def g4Func2(x, a4):
    g4 = (g4Func1(-x, a4) - g4Func1(x, a4)) / 2
    return g4

def g2Func1(x, a2, r):
    g2 = np.sin(x + a2) - r * np.sin(2 * (x + a2))
    return g2

def g2Func2(x, a2, r):
    g2 = (g2Func1(-x, a2, r) - g2Func1(x, a2, r)) / 2
    return g2

def CoupleFunc(x):
    return np.sin(x)
    #return (1 - np.cos(x))
    #return (x + np.pi) % (2 * np.pi) - np.pi

def coupleBickSystem(K, Eps, a2, a4, r):
    # bidirectional coupling in all pairs x_i - y_i
    def rhs(t, X):
        # we assume that X is arranged as
        # x1, x2, x3, y1, y2, y3
        Y1s = X[0:3]
        Y2s = X[3:6]
        Y1ders = singeBickSystem(K, a2, a4, r)(t, Y1s)
        Y2ders = singeBickSystem(-K, a2, a4, r)(t, Y2s)
        Y1ders[0] += Eps*(CoupleFunc(Y2s[0]-Y1s[0]))
        Y1ders[1] += Eps*(CoupleFunc(Y2s[1]-Y1s[1]))
        Y1ders[2] += Eps*(CoupleFunc(Y2s[2]-Y1s[2]))
        Y2ders[0] += Eps*(CoupleFunc(Y1s[0]-Y2s[0]))
        Y2ders[1] += Eps*(CoupleFunc(Y1s[1]-Y2s[1]))
        Y2ders[2] += Eps*(CoupleFunc(Y1s[2]-Y2s[2]))
        ders = Y1ders + Y2ders
        return ders
    return rhs

def coupleBickSystem_Simplifyed(K, Eps, a2, a4, r):
    def rhs(t, X):
        x1, x2, x3, y1, y2, y3 = X
        dX = [Eps - 2*np.cos(a2)*np.sin(x1) - Eps*np.cos(x1)*np.cos(y1) - Eps*np.sin(x1)*np.sin(y1) - K*np.cos(a4)*np.cos(x2)*np.sin(x1) + K*np.cos(a4)*np.cos(x3)*np.sin(x1) + 4*r*np.cos(x1)*np.sin(x1)*(2*np.cos(a2)**2 - 1),
    Eps - 2*np.cos(a2)*np.sin(x2) - Eps*np.cos(x2)*np.cos(y2) - Eps*np.sin(x2)*np.sin(y2) + K*np.cos(a4)*np.cos(x1)*np.sin(x2) - K*np.cos(a4)*np.cos(x3)*np.sin(x2) + 4*r*np.cos(x2)*np.sin(x2)*(2*np.cos(a2)**2 - 1),
    Eps - 2*np.cos(a2)*np.sin(x3) - Eps*np.cos(x3)*np.cos(y3) - Eps*np.sin(x3)*np.sin(y3) - K*np.cos(a4)*np.cos(x1)*np.sin(x3) + K*np.cos(a4)*np.cos(x2)*np.sin(x3) + 4*r*np.cos(x3)*np.sin(x3)*(2*np.cos(a2)**2 - 1),
    Eps - 2*np.cos(a2)*np.sin(y1) - Eps*np.cos(x1)*np.cos(y1) - Eps*np.sin(x1)*np.sin(y1) + K*np.cos(a4)*np.cos(y2)*np.sin(y1) - K*np.cos(a4)*np.cos(y3)*np.sin(y1) + 4*r*np.cos(y1)*np.sin(y1)*(2*np.cos(a2)**2 - 1),
    Eps - 2*np.cos(a2)*np.sin(y2) - Eps*np.cos(x2)*np.cos(y2) - Eps*np.sin(x2)*np.sin(y2) - K*np.cos(a4)*np.cos(y1)*np.sin(y2) + K*np.cos(a4)*np.cos(y3)*np.sin(y2) + 4*r*np.cos(y2)*np.sin(y2)*(2*np.cos(a2)**2 - 1),
    Eps - 2*np.cos(a2)*np.sin(y3) - Eps*np.cos(x3)*np.cos(y3) - Eps*np.sin(x3)*np.sin(y3) + K*np.cos(a4)*np.cos(y1)*np.sin(y3) - K*np.cos(a4)*np.cos(y2)*np.sin(y3) + 4*r*np.cos(y3)*np.sin(y3)*(2*np.cos(a2)**2 - 1)]
        return dX
    return rhs

def coupleBickSystem_RAW(K, Eps, a2, a4, r):
    def rhs(t, X):
        x1, x2, x3, y1, y2, y3 = X
        dX = [np.sin(a2 - x1) - np.sin(a2 + x1) - Eps*(np.sin(x1 - y1)) - (K*(np.sin(a4 + x1 + x2)/2 + np.sin(a4 + x1 - x2)/2 - np.sin(a4 - x1 + x2)/2 + np.sin(x1 - a4 + x2)/2))/2 + (K*(np.sin(a4 + x1 + x3)/2 + np.sin(a4 + x1 - x3)/2 - np.sin(a4 - x1 + x3)/2 + np.sin(x1 - a4 + x3)/2))/2 - r*np.sin(2*a2 - 2*x1) + r*np.sin(2*a2 + 2*x1),
              np.sin(a2 - x2) - np.sin(a2 + x2) - Eps*(np.sin(x2 - y2)) + (K*(np.sin(a4 + x1 + x2)/2 - np.sin(a4 + x1 - x2)/2 + np.sin(a4 - x1 + x2)/2 + np.sin(x1 - a4 + x2)/2))/2 - (K*(np.sin(a4 + x2 + x3)/2 + np.sin(a4 + x2 - x3)/2 - np.sin(a4 - x2 + x3)/2 + np.sin(x2 - a4 + x3)/2))/2 - r*np.sin(2*a2 - 2*x2) + r*np.sin(2*a2 + 2*x2),
              np.sin(a2 - x3) - np.sin(a2 + x3) - Eps*(np.sin(x3 - y3)) - (K*(np.sin(a4 + x1 + x3)/2 - np.sin(a4 + x1 - x3)/2 + np.sin(a4 - x1 + x3)/2 + np.sin(x1 - a4 + x3)/2))/2 + (K*(np.sin(a4 + x2 + x3)/2 - np.sin(a4 + x2 - x3)/2 + np.sin(a4 - x2 + x3)/2 + np.sin(x2 - a4 + x3)/2))/2 - r*np.sin(2*a2 - 2*x3) + r*np.sin(2*a2 + 2*x3),
              np.sin(a2 - y1) - np.sin(a2 + y1) - Eps*(np.sin(x1 - y1)) - r*np.sin(2*a2 - 2*y1) + r*np.sin(2*a2 + 2*y1) + (K*(np.sin(a4 + y1 + y2)/2 + np.sin(a4 + y1 - y2)/2 - np.sin(a4 - y1 + y2)/2 + np.sin(y1 - a4 + y2)/2))/2 - (K*(np.sin(a4 + y1 + y3)/2 + np.sin(a4 + y1 - y3)/2 - np.sin(a4 - y1 + y3)/2 + np.sin(y1 - a4 + y3)/2))/2,
              np.sin(a2 - y2) - np.sin(a2 + y2) - Eps*(np.sin(x2 - y2)) - r*np.sin(2*a2 - 2*y2) + r*np.sin(2*a2 + 2*y2) - (K*(np.sin(a4 + y1 + y2)/2 - np.sin(a4 + y1 - y2)/2 + np.sin(a4 - y1 + y2)/2 + np.sin(y1 - a4 + y2)/2))/2 + (K*(np.sin(a4 + y2 + y3)/2 + np.sin(a4 + y2 - y3)/2 - np.sin(a4 - y2 + y3)/2 + np.sin(y2 - a4 + y3)/2))/2,
              np.sin(a2 - y3) - np.sin(a2 + y3) - Eps*(np.sin(x3 - y3)) - r*np.sin(2*a2 - 2*y3) + r*np.sin(2*a2 + 2*y3) + (K*(np.sin(a4 + y1 + y3)/2 - np.sin(a4 + y1 - y3)/2 + np.sin(a4 - y1 + y3)/2 + np.sin(y1 - a4 + y3)/2))/2 - (K*(np.sin(a4 + y2 + y3)/2 - np.sin(a4 + y2 - y3)/2 + np.sin(a4 - y2 + y3)/2 + np.sin(y2 - a4 + y3)/2))/2]
        return dX
    return rhs

def coupleBickSystem_SympComb(K, Eps, a2, a4, r):
    def rhs(t, X):
        x1, x2, x3, y1, y2, y3 = X
        dX = [np.sin(a2 - x1) - np.sin(a2 + x1) - (K*np.sin(a4 + x1 - x2))/4 + (K*np.sin(a4 - x1 + x2))/4 - (K*np.sin(x1 - a4 + x2))/4 + (K*np.sin(a4 + x1 - x3))/4 - (K*np.sin(a4 - x1 + x3))/4 + (K*np.sin(x1 - a4 + x3))/4 - Eps*np.sin(x1 - y1) - (K*np.sin(a4 + x1 + x2))/4 + (K*np.sin(a4 + x1 + x3))/4 + 2*r*np.cos(2*a2)*np.sin(2*x1),
              np.sin(a2 - x2) - np.sin(a2 + x2) - (K*np.sin(a4 + x1 - x2))/4 + (K*np.sin(a4 - x1 + x2))/4 + (K*np.sin(x1 - a4 + x2))/4 - (K*np.sin(a4 + x2 - x3))/4 + (K*np.sin(a4 - x2 + x3))/4 - (K*np.sin(x2 - a4 + x3))/4 - Eps*np.sin(x2 - y2) + (K*np.sin(a4 + x1 + x2))/4 - (K*np.sin(a4 + x2 + x3))/4 + 2*r*np.cos(2*a2)*np.sin(2*x2),
              np.sin(a2 - x3) - np.sin(a2 + x3) + (K*np.sin(a4 + x1 - x3))/4 - (K*np.sin(a4 - x1 + x3))/4 - (K*np.sin(x1 - a4 + x3))/4 - (K*np.sin(a4 + x2 - x3))/4 + (K*np.sin(a4 - x2 + x3))/4 + (K*np.sin(x2 - a4 + x3))/4 - Eps*np.sin(x3 - y3) - (K*np.sin(a4 + x1 + x3))/4 + (K*np.sin(a4 + x2 + x3))/4 + 2*r*np.cos(2*a2)*np.sin(2*x3),
              np.sin(a2 - y1) - np.sin(a2 + y1) + (K*np.sin(a4 + y1 - y2))/4 - (K*np.sin(a4 - y1 + y2))/4 + (K*np.sin(y1 - a4 + y2))/4 - (K*np.sin(a4 + y1 - y3))/4 + (K*np.sin(a4 - y1 + y3))/4 - (K*np.sin(y1 - a4 + y3))/4 - Eps*np.sin(x1 - y1) + (K*np.sin(a4 + y1 + y2))/4 - (K*np.sin(a4 + y1 + y3))/4 + 2*r*np.cos(2*a2)*np.sin(2*y1),
              np.sin(a2 - y2) - np.sin(a2 + y2) + (K*np.sin(a4 + y1 - y2))/4 - (K*np.sin(a4 - y1 + y2))/4 - (K*np.sin(y1 - a4 + y2))/4 + (K*np.sin(a4 + y2 - y3))/4 - (K*np.sin(a4 - y2 + y3))/4 + (K*np.sin(y2 - a4 + y3))/4 - Eps*np.sin(x2 - y2) - (K*np.sin(a4 + y1 + y2))/4 + (K*np.sin(a4 + y2 + y3))/4 + 2*r*np.cos(2*a2)*np.sin(2*y2),
              np.sin(a2 - y3) - np.sin(a2 + y3) - (K*np.sin(a4 + y1 - y3))/4 + (K*np.sin(a4 - y1 + y3))/4 + (K*np.sin(y1 - a4 + y3))/4 + (K*np.sin(a4 + y2 - y3))/4 - (K*np.sin(a4 - y2 + y3))/4 - (K*np.sin(y2 - a4 + y3))/4 - Eps*np.sin(x3 - y3) + (K*np.sin(a4 + y1 + y3))/4 - (K*np.sin(a4 + y2 + y3))/4 + 2*r*np.cos(2*a2)*np.sin(2*y3)]
        return dX
    return rhs

def coupleBickSystemJac(K, Eps, a2, a4, r):
    def jac(t, X):
        x1, x2, x3, y1, y2, y3 = X
        tmp = [[Eps*np.cos(y1)*np.sin(x1) - Eps*np.cos(x1)*np.sin(y1) - 2*np.cos(a2)*np.cos(x1) + 4*r*(2*np.cos(a2)**2 - 1)*(2*np.cos(x1)**2 - 1) - K*np.cos(a4)*np.cos(x1)*np.cos(x2) + K*np.cos(a4)*np.cos(x1)*np.cos(x3), K*np.cos(a4)*np.sin(x1)*np.sin(x2), -K*np.cos(a4)*np.sin(x1)*np.sin(x3), -Eps*np.sin(x1 - y1), 0, 0],
 [-K*np.cos(a4)*np.sin(x1)*np.sin(x2), Eps*np.cos(y2)*np.sin(x2) - Eps*np.cos(x2)*np.sin(y2) - 2*np.cos(a2)*np.cos(x2) + 4*r*(2*np.cos(a2)**2 - 1)*(2*np.cos(x2)**2 - 1) + K*np.cos(a4)*np.cos(x1)*np.cos(x2) - K*np.cos(a4)*np.cos(x2)*np.cos(x3), K*np.cos(a4)*np.sin(x2)*np.sin(x3), 0, -Eps*np.sin(x2 - y2), 0],
 [K*np.cos(a4)*np.sin(x1)*np.sin(x3), -K*np.cos(a4)*np.sin(x2)*np.sin(x3), Eps*np.cos(y3)*np.sin(x3) - Eps*np.cos(x3)*np.sin(y3) - 2*np.cos(a2)*np.cos(x3) + 4*r*(2*np.cos(a2)**2 - 1)*(2*np.cos(x3)**2 - 1) - K*np.cos(a4)*np.cos(x1)*np.cos(x3) + K*np.cos(a4)*np.cos(x2)*np.cos(x3), 0, 0, -Eps*np.sin(x3 - y3)],
 [Eps*np.sin(x1 - y1), 0, 0, Eps*np.cos(x1)*np.sin(y1) - 2*np.cos(a2)*np.cos(y1) - Eps*np.cos(y1)*np.sin(x1) + 4*r*(2*np.cos(a2)**2 - 1)*(2*np.cos(y1)**2 - 1) + K*np.cos(a4)*np.cos(y1)*np.cos(y2) - K*np.cos(a4)*np.cos(y1)*np.cos(y3), -K*np.cos(a4)*np.sin(y1)*np.sin(y2), K*np.cos(a4)*np.sin(y1)*np.sin(y3)],
 [0, Eps*np.sin(x2 - y2), 0, K*np.cos(a4)*np.sin(y1)*np.sin(y2), Eps*np.cos(x2)*np.sin(y2) - 2*np.cos(a2)*np.cos(y2) - Eps*np.cos(y2)*np.sin(x2) + 4*r*(2*np.cos(a2)**2 - 1)*(2*np.cos(y2)**2 - 1) - K*np.cos(a4)*np.cos(y1)*np.cos(y2) + K*np.cos(a4)*np.cos(y2)*np.cos(y3), -K*np.cos(a4)*np.sin(y2)*np.sin(y3)],
 [0, 0, Eps*np.sin(x3 - y3), -K*np.cos(a4)*np.sin(y1)*np.sin(y3), K*np.cos(a4)*np.sin(y2)*np.sin(y3), Eps*np.cos(x3)*np.sin(y3) - 2*np.cos(a2)*np.cos(y3) - Eps*np.cos(y3)*np.sin(x3) + 4*r*(2*np.cos(a2)**2 - 1)*(2*np.cos(y3)**2 - 1) + K*np.cos(a4)*np.cos(y1)*np.cos(y3) - K*np.cos(a4)*np.cos(y2)*np.cos(y3)]]
        return tmp
    return jac









