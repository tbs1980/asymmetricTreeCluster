import numpy as np
import matplotlib.pyplot as plt

def plot_acc(file_list,label_list,title,outputFileName):
    """
    A function that plots the acceptance rate
    """

    if len(file_list) != len(label_list):
        raise ValueError("file_list and label_list should have the same length")

    for file_name,label_name in zip(file_list,label_list):
        acc = np.loadtxt(file_name,delimiter=",")
        plt.plot(acc[:,0],acc[:,1],label=label_name)


    plt.title(title)
    plt.legend(loc=0)
    plt.savefig(outputFileName)


#file_list = ['./acceptance2.dat','./acceptance10.dat',
#    './acceptance50.dat','./acceptance100.dat']
#label_list = ['2','10','50','100']
file_list = ['./acceptance2.dat']
label_list = ['2']
#title = 'lv = 1000, thr = 10'
title = 'lv = 1000, thr = 20'
#outputFileName = 'accVsDimLv1kThr10.png'
outputFileName = 'accVsDimLv1kThr20.png'

plot_acc(file_list,label_list,title,outputFileName)