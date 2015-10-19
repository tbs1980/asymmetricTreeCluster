import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib.cm as cm
import matplotlib as mpl

def plotTree(treeDumpFileName,dim1=0,dim2=1,radiusOfSphere=None):
    """
    A function for plotting the tree.

    @param treeDumpFileName file name of the tree-dump
    @param dim1 plot dimension 1
    @param dim2 plot dimension 2
    """
    tree = np.loadtxt(treeDumpFileName,delimiter=",")

    if len(tree.shape) == 1:
        tree = tree.reshape([1,tree.shape[0]])

    startCol = 3;
    numDims = (len(tree[0]) - startCol - 2)/2

    minWeight = np.min(tree[:,tree.shape[1]-1])
    maxWeight = np.max(tree[:,tree.shape[1]-1])

    norm = mpl.colors.Normalize(vmin=minWeight, vmax=maxWeight)
    cmap = cm.hot

    for nodeInfo in tree:
        b1Min = nodeInfo[startCol+dim1]
        b2Min = nodeInfo[startCol+dim2]

        b1Max = nodeInfo[startCol+numDims+dim1]
        b2Max = nodeInfo[startCol+numDims+dim2]

        plt.plot([b1Min,b1Max],[b2Min,b2Min],color='k')
        plt.plot([b1Min,b1Max],[b2Max,b2Max],color='k')
        plt.plot([b1Min,b1Min],[b2Min,b2Max],color='k')
        plt.plot([b1Max,b1Max],[b2Min,b2Max],color='k')

        numPoints = nodeInfo[2]
        if numPoints > 0 :
            weight = nodeInfo[tree.shape[1]-1]
            plt.fill_between([b1Min,b1Max],[b2Min,b2Min],[b2Max,b2Max],
                color=cm.ScalarMappable(norm=norm, cmap=cmap).to_rgba(weight))

    if radiusOfSphere != None:
        circle1=plt.Circle((0,0),radius=radiusOfSphere,color='m',fill=False)
        fig = plt.gcf()
        fig.gca().add_artist(circle1)


    plt.show()

if __name__ == "__main__" :
    if len(sys.argv) == 2:
        plotTree(sys.argv[1],0,1)
    elif len(sys.argv) == 3:
        plotTree(sys.argv[1],0,1,float(sys.argv[2]))
    else:
        print "usage: python ", sys.argv[0], "<tree-file-name>", "(radius)"
