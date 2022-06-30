import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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

def projPhase(X):
    return np.mod(X, 2*np.pi)
projPhase.label = lambda varName: r'${}\;\, {{\rm mod}} \; 2\pi$'.format(varName)

def projNone(X):
    return X
projNone.label = lambda varName: r'${}$'.format(varName)

#############################################################################

def plotTimeSeries(sol, projFunc, outFilePath, plotKwargs):
    fig = plt.figure(figsize=(20, 8))
    gs = gridspec.GridSpec(nrows=6, ncols=1, wspace=0.25, hspace=0.25)
    
    varNames = [r'\psi_1', r'\psi_2', r'\psi_3',
                r'\phi_1', r'\phi_2', r'\phi_3']
    
    for i in range(0, 6):
        ax = fig.add_subplot(gs[i, 0])
        ax.plot(sol.t, projFunc(sol.y[i]), **plotKwargs)
        plt.xlabel("$t$", fontsize=20)
        plt.ylabel(projFunc.label(varNames[i]), fontsize=20)

    # TIGHT LAYOUT!!!!!!!!!!!!
    plt.tight_layout()
    fig.savefig(outFilePath)
    return None

#############################################################################

def plotProjections(sol, projFunc, outFilePath, plotKwargs):
    fig = plt.figure(figsize=(20, 30))
    gs = gridspec.GridSpec(nrows=5, ncols=3, wspace=0.25, hspace=0.25)

    Pairs = [[0, 1], [0, 2], [1, 2],
             [3, 4], [3, 5], [4, 5],
             [0, 3], [0, 4], [0, 5],
             [1, 3], [1, 4], [1, 5],
             [2, 3], [2, 4], [2, 5],]

    varNames = [r'\psi_1', r'\psi_2', r'\psi_3',
                r'\phi_1', r'\phi_2', r'\phi_3']

    Labels = [(projFunc.label(varNames[i]), projFunc.label(varNames[j])) for i, j in Pairs]

    k = 0
    for i in range(0, 5):
        for j in range(0, 3):
            if k < 16 :
                ax = fig.add_subplot(gs[i, j])
                ax.plot(projFunc(sol.y[Pairs[k][0]]), projFunc(sol.y[Pairs[k][1]]), **plotKwargs)
                plt.xlabel(Labels[k][0], fontsize=20)
                plt.ylabel(Labels[k][1],fontsize=20)
                k= k + 1

    # TIGHT LAYOUT!!!!!!!!!!!!
    fig.savefig(outFilePath)
    return None


#############################################################################

def plotPoincare(sol, projFunc, outFilePath, plotKwargs):
    fig = plt.figure(figsize=(20, 30))
    gs = gridspec.GridSpec(nrows=5, ncols=3, wspace=0.25, hspace=0.25)

    Pairs = [[0, 1], [0, 2], [1, 2],
             [3, 4], [3, 5], [4, 5],
             [0, 3], [0, 4], [0, 5],
             [1, 3], [1, 4], [1, 5],
             [2, 3], [2, 4], [2, 5],]

    varNames = [r'\psi_1', r'\psi_2', r'\psi_3',
                r'\phi_1', r'\phi_2', r'\phi_3']

    Labels = [(projFunc.label(varNames[i]), projFunc.label(varNames[j])) for i, j in Pairs]

    k = 0
    for i in range(0, 5):
        for j in range(0, 3):
            if k < 16 :
                ax = fig.add_subplot(gs[i, j])
                ax.scatter(projFunc(sol.y_events[0][:, Pairs[k][0]]), projFunc(sol.y_events[0][:, Pairs[k][1]]), **plotKwargs)
                plt.xlabel(Labels[k][0], fontsize=20)
                plt.ylabel(Labels[k][1], fontsize=20)
                k= k + 1

    # TIGHT LAYOUT!!!!!!!!!!!!
    fig.savefig(outFilePath)
    return None


#############################################################################

def normNone(X):
    return X
normNone.label = lambda varName: r'${}$'.format(varName)

def normDefault(X):
    return (X - min(X)) / (max(X) - min(X))
normDefault.label = lambda varName: r'$||{}||$'.format(varName)


def plotreturnTime(sol, normFunc, outFilePath, plotKwargs):
    fig = plt.figure(figsize=(10, 5))

    if len(sol.t_events[0] > 1):
        t_diff = normFunc(sol.t_events[0][1:] - sol.t_events[0][0:-1])
        plt.scatter(sol.t_events[0][: -1], t_diff)
        
        if normFunc == normNone:
            plt.xlabel("$t$", fontsize=20)
            plt.ylabel("Время возврата", fontsize=20)
            
        else:
            plt.xlabel("$t$", fontsize=20)
            plt.ylabel("нормированное\n время возврата", fontsize=20)
            
        fig.savefig(outFilePath)
        
    else:
        print("Not enough points")
        
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


def plotActivationDiagram(sol, varNamesList, colorizer, outFileName):
    assert len(varNamesList) == len(sol.y), "Dimension of phase space and length of list with variable names do not coincide!"
    solMat = np.array(sol.y)
    clrdPlt = [[colorizer.normFunc(v) for v in row] for row in solMat]
    plt.close()
    fig = plt.figure(figsize=(15, 6))
    N = len(varNamesList)
    plt.pcolormesh(sol.t, range(N), clrdPlt, cmap=colorizer.cmap, norm=colorizer.norm, shading='nearest')
    plt.gca().set_yticks(0.5 + np.arange(N))
    plt.gca().set_yticklabels(varNamesList, fontsize=20)
    #plt.xlim([80000, 100000])
    plt.gca().tick_params(which="major", width=1.0, labelsize=20)
    plt.gca().tick_params(which="major", length=5, labelsize=20)
    plt.xlabel(r'$t$',fontsize=20)
    plt.tight_layout()
    plt.savefig(outFileName, facecolor='white')
    return fig











