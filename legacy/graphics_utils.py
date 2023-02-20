import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mpc
import matplotlib.cm as cm
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
            Y[i] = -np.log10(1e-15 + 1 - np.sin(X[i]/2))
        else:
            Y[i] = np.log10(1e-15 + 1 + np.sin(X[i]/2))
    return Y
projLogSin.label = lambda varName: r'$\pm\log{{(1\pm\sin{{\frac{{{}}}{{2}}}})}}$'.format(varName)

#############################################################################

def plotTimeSeries(sol, projFunc, plot_title, savePath, plotKwargs, plot_func = "plot"):
    if len(sol.y) == 6:
        fig, ax = plt.subplots(figsize=(16, 8), nrows=6, ncols = len(projFunc))

        varNames = [r'\psi_1', r'\psi_2', r'\psi_3',
                    r'\phi_1', r'\phi_2', r'\phi_3']

        for i in range(0, 6):
            for j in range(0, len(projFunc)):
                if plot_func == "scatter":
                    ax[i][j].scatter(sol.t, projFunc[j](sol.y[i]), **plotKwargs)
                else:
                    ax[i][j].plot(sol.t, projFunc[j](sol.y[i]), **plotKwargs)
                ax[i][j].set_xlabel("$t$", fontsize=14)
                ax[i][j].set_ylabel(projFunc[j].label(varNames[i]), fontsize=14, rotation = 60)

    elif len(sol.y) == 3:
        fig, ax = plt.subplots(figsize=(16, 6), nrows=3, ncols = len(projFunc))

        varNames = [r'\psi_1', r'\psi_2', r'\psi_3']

        for i in range(0, 3):
            for j in range(0, len(projFunc)):
                if plot_func == "scatter":
                    ax[i][j].scatter(sol.t, projFunc[j](sol.y[i]), **plotKwargs)
                else:
                    ax[i][j].plot(sol.t, projFunc[j](sol.y[i]), **plotKwargs)
                    
                ax[i][j].set_xlabel("$t$", fontsize=14)
                ax[i][j].set_ylabel(projFunc[j].label(varNames[i]), fontsize=14, rotation = 60)
            
    fig.suptitle(r"$\bf{Временные ~диаграммы}$" + '\n' + plot_title)  
    plt.tight_layout()
    fig.savefig(savePath)
    return (fig, ax)

#############################################################################

def plotProjections(sol, projFunc, plot_title, savePath, plotKwargs, plot_func = "plot"):
    if len(sol.y) == 6:
        fig, ax = plt.subplots(figsize=(8, 10), nrows=5, ncols=3)

        if Pairs == []:
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
                if k < len(Pairs) :
                    if plot_func == "scatter":
                        ax[i][j].scatter(projFunc(sol.y[Pairs[k][0]]), projFunc(sol.y[Pairs[k][1]]), **plotKwargs)
                    else:
                        ax[i][j].plot(projFunc(sol.y[Pairs[k][0]]), projFunc(sol.y[Pairs[k][1]]), **plotKwargs)
                        
                    ax[i][j].set_xlabel(Labels[k][0], fontsize=12)
                    ax[i][j].set_ylabel(Labels[k][1],fontsize=12)
                    ax[i][j].set_aspect('equal')
                    k= k + 1
    elif len(sol.y) == 3:
        fig, ax = plt.subplots(figsize=(8, 4), ncols=3)

        Pairs = [[0, 1], [0, 2], [1, 2]]
        varNames = [r'\psi_1', r'\psi_2', r'\psi_3']
        Labels = [(projFunc.label(varNames[i]), projFunc.label(varNames[j])) for i, j in Pairs]

        for k in [0, 1, 2]:
            ax[k].plot(projFunc(sol.y[Pairs[k][0]]), projFunc(sol.y[Pairs[k][1]]), **plotKwargs)
            ax[k].set_xlabel(Labels[k][0], fontsize=12)
            ax[k].set_ylabel(Labels[k][1],fontsize=12)
            ax[k].set_aspect('equal')

    fig.suptitle(r"$\bf{Проекции ~на ~плоскости}$" + '\n' + plot_title)
    plt.tight_layout()
    fig.savefig(savePath)
    return (fig, ax)


#############################################################################

def plotPoincare(sol, projFunc, plot_title, savePath, plotKwargs, trajectory_plot_func = "plot"):
    if len(sol.y) == 6:
        fig, ax = plt.subplots(figsize=(8, 10), nrows=5, ncols=3)

        Pairs = [[0, 1], [0, 2], [1, 2],
                 [3, 4], [3, 5], [4, 5],
                 [0, 3], [0, 4], [0, 5],
                 [1, 3], [1, 4], [1, 5],
                 [2, 3], [2, 4], [2, 5]]

        varNames = [r'\psi_1', r'\psi_2', r'\psi_3',
                    r'\phi_1', r'\phi_2', r'\phi_3']

        Labels = [(projFunc.label(varNames[i]), projFunc.label(varNames[j])) for i, j in Pairs]

        k = 0
        for i in range(0, 5):
            for j in range(0, 3):
                if k < len(Pairs):
                    ax[i][j].scatter(projFunc(sol.y_events[0][:, Pairs[k][0]]), projFunc(sol.y_events[0][:, Pairs[k][1]]), **plotKwargs)
                    if trajectory_plot_func == "plot":
                        ax[i][j].plot(projFunc(sol.y[Pairs[k][0]]), projFunc(sol.y[Pairs[k][1]]), alpha = 0.5)
                    elif trajectory_plot_func == "scatter":
                        ax[i][j].scatter(projFunc(sol.y[Pairs[k][0]]), projFunc(sol.y[Pairs[k][1]]), alpha = 0.5)
                    ax[i][j].set_xlabel(Labels[k][0], fontsize=12)
                    ax[i][j].set_ylabel(Labels[k][1], fontsize=12)
                    ax[i][j].set_aspect('equal')
                    k= k + 1
                
    elif len(sol.y) == 3:
        fig, ax = plt.subplots(figsize=(8, 4), nrows=1, ncols=3)

        Pairs = [[0, 1], [0, 2], [1, 2]]
        varNames = [r'\psi_1', r'\psi_2', r'\psi_3']
        Labels = [(projFunc.label(varNames[i]), projFunc.label(varNames[j])) for i, j in Pairs]

        for k in [0, 1, 2]:
            ax[k].scatter(projFunc(sol.y_events[0][:, Pairs[k][0]]), projFunc(sol.y_events[0][:, Pairs[k][1]]), **plotKwargs)
            ax[k].plot(projFunc(sol.y[Pairs[k][0]]), projFunc(sol.y[Pairs[k][1]]), alpha = 0.5)
            ax[k].set_xlabel(Labels[k][0], fontsize=12)
            ax[k].set_ylabel(Labels[k][1], fontsize=12)
            ax[k].set_aspect('equal')

    fig.suptitle(r"$\bf{Отображения ~Пуанкаре}$" + '\n' + plot_title)
    plt.tight_layout()
    fig.savefig(savePath)
    return (fig, ax)


#############################################################################

def normNone(X):
    return X
normNone.label = lambda varName: r'${}$'.format(varName)

def normDefault(X):
    return (X - min(X)) / (max(X) - min(X))
normDefault.label = lambda varName: r'$||{}||$'.format(varName)


def plotreturnTime(sol, normFunc, plot_title, savePath, plotKwargs):
    fig = plt.figure(figsize=(10, 4))

    if len(sol.t_events[0]) > 1:
        t_diff = normFunc(sol.t_events[0][1:] - sol.t_events[0][0:-1])
        plt.scatter(sol.t_events[0][: -1], t_diff)
        
        if normFunc.label == normNone.label:
            avg = np.mean(t_diff)
            max_min_diff = abs(max(t_diff) - min(t_diff))
            fig.gca().set_ylim([avg - max(max_min_diff/5, 10), avg + max(max_min_diff/5, 10)])
            plt.xlabel("$t$", fontsize=14)
            plt.ylabel("Время возврата", fontsize=14)
            
        else:
            plt.xlabel("$t$", fontsize=14)
            plt.ylabel("нормированное\n время возврата", fontsize=14)
        
    else:
        print("Not enough points")
        
    plt.title(r"$\bf{Времена ~возврата ~на ~секущую ~Пуанкаре}$" + '\n' + plot_title)
    plt.tight_layout()
    fig.savefig(savePath)
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


def plotActivationDiagram(sol, varNamesList, colorizer, plot_title, savePath):
    assert len(varNamesList) == len(sol.y), "Dimension of phase space and length of list with variable names do not coincide!"
    solMat = np.array(sol.y)
    clrdPlt = [[colorizer.normFunc(v) for v in row] for row in solMat]
    plt.close()
    fig = plt.figure(figsize=(15, 6))
    N = len(varNamesList)
    pcm = plt.pcolormesh(sol.t, range(N), clrdPlt, cmap=colorizer.cmap, norm=colorizer.norm, shading='nearest')
    cbar = fig.colorbar(pcm)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label(label = r"$\phi_i \approx 0$ -> черный; $\phi_i \in (0, \pi)$ -> серый; $\phi_i \approx pi$ -> белый", size=12) 
    plt.gca().set_yticks(0.5 + np.arange(N))
    plt.gca().set_yticklabels(varNamesList, fontsize=14)
    #plt.xlim([80000, 100000])
    plt.gca().tick_params(which="major", width=1.0, labelsize=14)
    plt.gca().tick_params(which="major", length=5, labelsize=14)
    plt.xlabel(r'$t$', fontsize=14)
    plt.title(r"$\bf{Диаграмма~активности}$" + '\n' + plot_title)
    plt.tight_layout()
    fig.savefig(savePath)
    return fig











