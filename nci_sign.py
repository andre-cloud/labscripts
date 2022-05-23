import argparse
from collections import Counter
import os
import sys
import shutil

from alive_progress import alive_bar
import cclib
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import numpy as np
import pickle

parser = argparse.ArgumentParser()
parser.add_argument('file', help='Scan log file(s)', nargs='+')
parser.add_argument('-ns', '--negative_sign', default=-0.1, type=float)
parser.add_argument('-ps', '--positive_sign', default=0.1, type=float)
parser.add_argument('-th', '--threshold', default=.4, type=float)
parser.add_argument('-tit', '--title', help='define the title of the graph. Default %(default)s', default=r's($\rho$) Plot')
parser.add_argument('-c', '--colors', nargs='+', help='Declaire the colors of the graph')



args = parser.parse_args()




colors = ['#382EFF', '#2CDB1D', '#F51118']
cmap_name='BGR'
cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=1000)

def read_sign(filename):
    data = np.loadtxt(filename)
    x = data[:, 0]
    y = data[:, 1]

    d = {k:v for k,v in zip(x,y)}
    colors = {x:d[x] for x in d if (args.negative_sign <= x <= args.positive_sign and d[x] <= args.threshold)}
    not_colors = {x:d[x] for x in d if x not in colors}
    
    return colors, not_colors

def show_graph(colors, not_colors):
    fig, (ax, cax)  = plt.subplots(nrows=2,figsize=(4,4), 
                  gridspec_kw={"height_ratios":[1, 0.02]})

    z = ax.scatter(list(colors.keys()),list(colors.values()),  c=list(colors.keys()), cmap=cmap)
    ax.scatter(list(not_colors.keys()),list(not_colors.values()))
    fig.colorbar(z, cax=cax, orientation="horizontal")



    ax.set_title(args.title)
    ax.set_xlabel(r'$sign(\lambda_2)\rho$ (a.u.)')
    ax.set_ylabel('Reduced Gradient')




    plt.tight_layout()


    plt.show()



if __name__=='__main__':
    for i in args.file:
        colors, not_colors = read_sign(i)
        show_graph(colors, not_colors)