import argparse
import sys

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
np.set_printoptions(threshold=sys.maxsize)


parser = argparse.ArgumentParser()
parser.add_argument('-a', '--animate', action='store_true')

args = parser.parse_args()

fig, (ax, cax) = plt.subplots(nrows=2, figsize=(6,5), gridspec_kw={"height_ratios":[1, 0.02], "hspace":0.3})


def get_data():
    z = np.load('en.npy')
    print(z)
    x, y = np.load('grid.npy')
    # z = (z-np.min(z))*627.51
    z = (z--171.612669228135)*627.51
    z = z[:-1, :-1]
    return x, y, z



def animate(i):
    x, y, z = get_data()
    ax.clear()
    plot(x, y, z)


def plot(x, y, z):
    z_min, z_max = -20, 20 # 
    # z_min, z_max = np.abs(z).min(), (np.abs(z).min()+60 if z[z>100].any() else np.max(z))

    c = ax.pcolormesh(x, y, z, cmap='coolwarm', vmin=z_min, vmax= z_max)
    ax.set_title('2D PES around N-H and C-N bond formation')
    ax.set_xlabel('C-C distance (67-87) [Å]')
    ax.set_ylabel('C-C distance (68-86) [Å]')

    # set the limits of the plot to the limits of the data
    ax.axis([x.min(), x.max(), y.min(), y.max()])
    fig.colorbar(c, cax=cax, orientation="horizontal")

    
if args.animate:
    ani = animation.FuncAnimation(fig, animate, interval=10000)
else:
    x, y, z = get_data()
    plot(x, y, z)
    

plt.show()