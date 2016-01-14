import numpy as np
import matplotlib.pyplot as plt

pts = np.loadtxt("./points.dat")
plt.plot(pts[:,0],pts[:,1],'o',color='blue')
plt.xlim(-5,5)
plt.ylim(-5,5)
fig = plt.gca()
fig.set_axis_bgcolor('#B0B1B2')
plt.show()