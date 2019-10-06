import pylab as p
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LightSource
from matplotlib import cm
import imageio
import glob, os, sys

import time

path = "./plots/"

# Remove old plotfiles
for filename in glob.glob(path+"wave*.png"):
    os.remove(filename)

# Load arrays
u = np.load("u.npy")
h = np.load("h.npy")
y = np.load("y.npy")
x = np.load("x.npy")

T = np.shape(u)[2]

X, Y = np.meshgrid(x, y)

# Plot for each timestep
for t in range(T):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    #Plot the wave
    surf = ax.plot_surface(X, Y, u[:,:,t], rstride=1, cstride=1, cmap=cm.get_cmap("winter_r"),#color="blue",
                            linewidth=0, antialiased=True, shade=True, 
                            vmin=p.amin(u),vmax=p.amax(u))
    
    #Plot the bottom
    surf2 = ax.plot_surface(X, Y, h, rstride=1, cstride=1, color="brown",
                                linewidth=0, antialiased=True, shade=True,
                                alpha=0.6)
    ax.set_zlim3d(p.amin(h), p.amax(u))
    fig.colorbar(surf)
    plt.savefig(path + "wave%05d"  % (t) + ".png")
    

# Create gif
filenames = glob.glob('./plots/*.png')
filenames.sort()
with imageio.get_writer('./wave.gif', mode='I') as writer:
    for filename in filenames:
        writer.append_data(imageio.imread((filename)))
