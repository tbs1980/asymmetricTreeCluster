import matplotlib.pyplot as plt
import numpy as np
from math import exp
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.mlab import bivariate_normal
import matplotlib.cm as cm

def unit_gauss_func_2d(x,y):
    return exp(-0.5*np.sum(x*y))

x = np.linspace(-5, 5, 300)
y = x
X,Y = np.meshgrid(x, y)
Z = bivariate_normal(X,Y)

"""
plt.contour(X,Y,Z,[1.9641280346397437e-05,0.0028007247208697364,0.058549831524319168,0.12394999430965298])
plt.imshow(Z, interpolation='bilinear', origin='lower',cmap=cm.Spectral_r, extent=(-5, 5, -5, 5))
plt.plot(-3,-3,'bo')
plt.plot(2.2,-1.8,'bo')
plt.plot(-1,1,'bo')
plt.plot(-0.5,-0.5,'bo')
#plt.show()
plt.savefig("biGaussMapCont.png",transparent=True)
"""

plt.imshow(Z, interpolation='bilinear', origin='lower',cmap=cm.Spectral_r, extent=(-5, 5, -5, 5))
pts = np.loadtxt("./points.dat")
plt.plot(pts[:,0],pts[:,1],'o',color='blue')
plt.xlim(-5,5)
plt.ylim(-5,5)
plt.savefig("biGaussMapPoints.png",transparent=True)


"""
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z,cmap=cm.Spectral_r)
ax.set_alpha(0.1)
plt.savefig("biGauss.png",transparent=True)
"""