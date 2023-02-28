import matplotlib.pyplot as plt
import matplotlib.colors as mpc
import numpy as np

#############################################################################

def projHalfCos(X):
    return np.cos(X/2)
projHalfCos.label= lambda varName: r'$\cos{{\frac{{{}}}{{2}}}}$'.format(varName)

def projHalfSin(X):
    return np.sin(X/2)
projHalfSin.label= lambda varName: r'$\sin{{\frac{{{}}}{{2}}}}$'.format(varName)

def projSin(X):
    return np.sin(X)
projSin.label= lambda varName: r'$\sin{{{}}}$'.format(varName)

def projCos(X):
    return np.cos(X)
projCos.label= lambda varName: r'$\cos{{{}}}$'.format(varName)

def projArctan2(X):
    return np.arctan2(np.sin(X), np.cos(X))
projArctan2.label= lambda varName: r'$\arctan2{{{}}}$'.format(varName)

def projPhase(X):
    return np.mod(X, 2*np.pi)
projPhase.label = lambda varName: r'${}\;\, {{\rm mod}} \; 2\pi$'.format(varName)

def projNone(X):
    return X
projNone.label = lambda varName: r'${}$'.format(varName)

def projLogSin(X):
    Y = [0] * len(X)
    for i in range(0, len(Y)):
        if np.sin(X[i]/2) > 0:
            Y[i] = -np.log(1e-15 + 1 - np.sin(X[i]/2))
        else:
            Y[i] = np.log(1e-15 + 1 + np.sin(X[i]/2))
    return Y
projLogSin.label = lambda varName: r'$\pm\log{{(1\pm\sin{{\frac{{{}}}{{2}}}})}}$'.format(varName)

#############################################################################

def plotTimeSeries(sol, savePath = None, varNames = None, projFunc = None, title = None, plotKwargs = None): 
    if varNames == None : varNames = [f"y_{p}" for p in range(1, len(sol.y) + 1)]
    if projFunc == None : projFunc = projNone
    if title == None : title = r'$\bf{Временные ~диаграммы}$'
    if plotKwargs == None : plotKwargs = {}

    fig, axes = plt.subplots(len(sol.y), 1, sharex = True, figsize = (6, len(sol.y)))
    for i in range(0, len(sol.y)):
        axes[i].plot(sol.t, projFunc(sol.y[i]), **plotKwargs)
        axes[i].set_ylabel(projFunc.label(varNames[i]), fontsize=16, rotation = 60)
    axes[-1].set_xlabel("t", fontsize=16)

    fig.suptitle(title)
    plt.tight_layout()
    if savePath != None : fig.savefig(savePath)

    return fig
    

#############################################################################

def plotPoincare(sol, varPairs, /, savePath = None, varNames = None, showEvents = False, projFunc = None, title = None, plotKwargs = None):
    if varNames == None : varNames = [f"y_{p}" for p in range(1, len(sol.y) + 1)]
    if projFunc == None : projFunc = projNone
    if title == None : title = r"$\bf{Отображения ~Пуанкаре}$"
    if plotKwargs == None : plotKwargs = {}
    
    pairsRows = np.shape(varPairs)[0]
    pairsCols = np.shape(varPairs)[1]
    figsizeScale = 2
    
    fig, axes = plt.subplots(pairsRows, pairsCols, sharex = True, sharey = True, figsize = (figsizeScale * pairsCols + 1, figsizeScale * pairsRows))
    if pairsRows == 1:
        axes = [axes]

    for i in range(0, pairsRows):
        for j in range(0, pairsCols):
            firstVarInd = varPairs[i][j][0]
            secondVarInd = varPairs[i][j][1]
            axes[i][j].plot(projFunc(sol.y[firstVarInd]), projFunc(sol.y[secondVarInd]), alpha = 0.5)
            if (len(sol.t_events[0]) > 0) & showEvents :
                axes[i][j].scatter(projFunc(sol.y_events[0][:, firstVarInd]), projFunc(sol.y_events[0][:, secondVarInd]), **plotKwargs)
            axes[i][j].set_xlabel(projFunc.label(varNames[firstVarInd]), fontsize=16)
            axes[i][j].set_ylabel(projFunc.label(varNames[secondVarInd]), fontsize=16)
            axes[i][j].set_aspect('equal')

    fig.suptitle(title)
    plt.tight_layout()
    if savePath != None : fig.savefig(savePath)

    return fig

#############################################################################

def normNone(X):
    return X
normNone.label = lambda varName: r'${}$'.format(varName)
normNone.title = "Времена возврата"

def normDefault(X):
    return (X - min(X)) / (max(X) - min(X))
normDefault.label = lambda varName: r'$||{}||$'.format(varName)
normDefault.title = "Нормированные времена возврата"


def plotReturnTime(sol, /, savePath = None, normFunc = None, title = None, plotKwargs = None):
    if normFunc == None : normFunc = normDefault
    if title == None : title = r"$\bf{Времена ~возврата ~на ~секущую ~Пуанкаре}$"
    if plotKwargs == None : plotKwargs = {}
    
    fig = plt.figure(figsize=(10, 4))
    assert len(sol.t_events[0]) > 1, f"Not enough event points in plotReturnTime"

    t_diff = normFunc(sol.t_events[0][1:] - sol.t_events[0][0:-1])
    plt.scatter(sol.t_events[0][: -1], t_diff, **plotKwargs)

    plt.xlabel("$t$", fontsize=14)
    plt.ylabel(normFunc.title, fontsize=14)

    plt.title(title)
    plt.tight_layout()
    if savePath != None : fig.savefig(savePath)

    return fig

#############################################################################

def circleDist(x, y):
    xStd = x % (2 * np.pi)
    yStd = y % (2 * np.pi)
    return min(np.abs(xStd - yStd), 2 * np.pi - np.abs(xStd - yStd))


def normValue(x):
    eps = 1e-2
    label = 1.5
    if circleDist(x, 0.) < eps:
        label = 0.5  # white
    elif circleDist(x, np.pi) < eps:
        label = 2.5  # black

    return label


class ColorActivation:
    def __init__(self, colorList, normalizeFunc):
        assert len(colorList) == 3, 'Must be three colors!'
        self.colorList = colorList
        # create colormap to use further
        myCmap = mpc.ListedColormap(colorList)
        boundaries = [0, 1, 2, 3]
        myNorm = mpc.BoundaryNorm(boundaries, myCmap.N, clip=True)
        # assign to fields
        self.cmap = myCmap
        self.norm = myNorm
        self.normFunc = normalizeFunc


def plotActivationDiagram(sol, varNames = None, colorizer = None, title = None, savePath = None):
    if varNames == None : varNames = [f"y_{{{p}}}" for p in range(1, len(sol.y) + 1)]
    if colorizer == None : colorizer = ColorActivation([(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)], normValue)
    if title == None : title = r"$\bf{Диаграмма~активности}$"

    solMat = np.array(sol.y)
    clrdPlt = [[colorizer.normFunc(v) for v in row] for row in solMat]
    plt.close()
    fig = plt.figure(figsize=(15, 6))
    N = len(varNames)
    pcm = plt.pcolormesh(sol.t, range(N), clrdPlt, cmap=colorizer.cmap, norm=colorizer.norm, shading='nearest')
    cbar = fig.colorbar(pcm)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label(label = r"$\phi_i \approx 0$ -> черный; $\phi_i \in (0, \pi)$ -> серый; $\phi_i \approx pi$ -> белый", size=12) 
    plt.gca().set_yticks(0.5 + np.arange(N))
    plt.gca().set_yticklabels([projNone.label(var) for var in varNames], fontsize=14)

    plt.gca().tick_params(which="major", width=1.0, labelsize=14)
    plt.gca().tick_params(which="major", length=5, labelsize=14)
    plt.xlabel(r'$t$', fontsize=14)

    plt.title(title)
    plt.tight_layout()
    if savePath != None : fig.savefig(savePath)

    return fig

def plotActivationDiagram_continuos(sol, varNames = None, cmap = None, title = None, savePath = None):
    if varNames == None : varNames = [f"y_{{{p}}}" for p in range(1, len(sol.y) + 1)]
    if cmap == None : cmap = plt.colormaps['Greys']
    if title == None : title = r"$\bf{Диаграмма~активности}$"

    solMat = np.array(sol.y)
    clrdPlt = [[circleDist(y, 0.0) for y in s]for s in solMat]
    plt.close()
    fig = plt.figure(figsize=(15, 6))
    N = len(varNames)
    pcm = plt.pcolormesh(sol.t, range(N), clrdPlt, cmap=cmap, shading='nearest', vmin=0.0, vmax=np.pi)
    cbar = fig.colorbar(pcm)
    cbar.ax.tick_params(labelsize=12)
    plt.gca().set_yticks(0.5 + np.arange(N))
    plt.gca().set_yticklabels([projNone.label(var) for var in varNames], fontsize=14)

    plt.gca().tick_params(which="major", width=1.0, labelsize=14)
    plt.gca().tick_params(which="major", length=5, labelsize=14)
    plt.xlabel(r'$t$', fontsize=14)

    plt.title(title)
    plt.tight_layout()
    if savePath != None : fig.savefig(savePath)

    return fig