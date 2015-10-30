from __future__ import print_function
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
    numDims = (len(tree[0]) - startCol - 4)/2

    print("numDims ",numDims)

    wtMin = tree[:,tree.shape[1]-4]
    wtMax = tree[:,tree.shape[1]-3]
    numPts = tree[:,3]

    #minWeight = np.min(wtMin[numPts>0])
    #maxWeight = np.max(wtMax[numPts>0])

    minWeight = np.min(wtMin)
    maxWeight = np.max(wtMax)

    print(minWeight, maxWeight)

    norm = mpl.colors.Normalize(vmin=minWeight, vmax=maxWeight)
    cmap = cm.hot

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

        numPoints = nodeInfo[2]
        if numPoints > 0 :
            treeInd = int(nodeInfo[0])
            #plt.text(0.5*(b1Min+b1Max),0.5*(b2Min+b2Max),str(treeInd))
            weight = nodeInfo[tree.shape[1]-2] # plotting the mean weight
            plt.fill_between([b1Min,b1Max],[b2Min,b2Min],[b2Max,b2Max],
                color=cm.ScalarMappable(norm=norm, cmap=cmap).to_rgba(weight))

    if radiusOfSphere != None:
        circle1=plt.Circle((0,0),radius=radiusOfSphere,color='m',fill=False)
        fig = plt.gcf()
        fig.gca().add_artist(circle1)

    #ax = plt.gca()
    #ax.set_axis_bgcolor('red')
    plt.show()

if __name__ == "__main__" :
    if len(sys.argv) == 2:
        plotTree(sys.argv[1],0,1)
    elif len(sys.argv) == 3:
        plotTree(sys.argv[1],0,1,float(sys.argv[2]))
    else:
        print("usage: python ", sys.argv[0], "<tree-file-name>", "(radius)")
