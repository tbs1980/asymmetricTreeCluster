from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib.cm as cm
import matplotlib as mpl
from math import exp

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

    # col 0 tree-index
    # col 1 split-dimension
    # col 2 points-size
    # col 3 weight min
    # col 4 weight max
    # col 5 weight mean
    # col 6 weight std-dvn
    # col 7 volume of the node
    # col 8 node-characterstic
    # col 9 acceptance ratio of the node
    # ----------------------------------
    # col 10 -> end node min and node max

    startCol = 10;

    numDims = (len(tree[0]) - startCol)/2

    print("numDims ",numDims)

    numPts = tree[:,2]
    wtMin = tree[:,3]
    wtMax = tree[:,4]

    # what column should be used for heat map
    heatMapPropertyCol = 3

    minWeight = exp(np.min(tree[np.where(numPts>0) ,heatMapPropertyCol]))
    maxWeight = exp(np.max(tree[np.where(numPts>0) ,heatMapPropertyCol]))

    norm = mpl.colors.Normalize(vmin=minWeight, vmax=maxWeight)
    cmap = cm.Spectral_r

    for nodeInfo in tree:
        #if int(nodeInfo[1]) == int(dim1) or int(nodeInfo[1]) == int(dim2):
        b1Min = nodeInfo[startCol+dim1]
        b2Min = nodeInfo[startCol+dim2]

        b1Max = nodeInfo[startCol+numDims+dim1]
        b2Max = nodeInfo[startCol+numDims+dim2]

        plt.plot([b1Min,b1Max],[b2Min,b2Min],color='k')
        plt.plot([b1Min,b1Max],[b2Max,b2Max],color='k')
        plt.plot([b1Min,b1Min],[b2Min,b2Max],color='k')
        plt.plot([b1Max,b1Max],[b2Min,b2Max],color='k')

        numPoints = nodeInfo[2] # col 2 points-size
        if numPoints > 0 :
            treeInd = int(nodeInfo[0]) # col 0 tree-index

            weight = exp(nodeInfo[heatMapPropertyCol]) # plotting the mean weight

            plt.fill_between([b1Min,b1Max],[b2Min,b2Min],[b2Max,b2Max],
                color=cm.ScalarMappable(norm=norm, cmap=cmap).to_rgba(weight))

    if radiusOfSphere != None:
        circle1=plt.Circle((0,0),radius=radiusOfSphere,color='m',fill=False)
        fig = plt.gcf()
        fig.gca().add_artist(circle1)

    fig = plt.gca()
    fig.set_axis_bgcolor('#B0B1B2')

    #plt.show()
    #plt.savefig('plotTree.svg')
    plt.axes().set_aspect('equal', 'datalim')
    plt.axes().set_xlim(-5,5)
    plt.axes().set_ylim(-5,5)
    plt.savefig(treeDumpFileName.replace("dat","png"),transparent=True)

if __name__ == "__main__" :
    if len(sys.argv) == 2:
        plotTree(sys.argv[1],0,1)
    elif len(sys.argv) == 3:
        plotTree(sys.argv[1],0,1,float(sys.argv[2]))
    else:
        print("usage: python ", sys.argv[0], "<tree-file-name>", "(radius)")
