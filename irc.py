import argparse
import os
import sys
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pickle


parser = argparse.ArgumentParser()

parser.add_argument('file', help='Scan log file(s)', nargs='+')
parser.add_argument('-tit', '--title', help='define the title of the graph. Default %(default)s', default='IRC')

parser.add_argument('-pr','--prod_right', help='True is the product is on the right of the graph', action='store_false')
parser.add_argument('-a', '--arrow', help='Toggle the arrow on the graph', action='store_true')
parser.add_argument('-yp', '--y_pad', type=int, default=0)

parser.add_argument('--save', help='Save pickle and csvs of the graph', action='store_true')
parser.add_argument('-gd','--graph_directory', help='Define the directory in which you want to save the files of the graph. Default: %(default)s', default='irc_graph')

parser.add_argument('-c', '--color', help='Define the color of the graph (see matplotlib docs for the table of colors). Default %(default)s', default='blue')

args = parser.parse_args()


def get_irc(filename):
    
    with open(filename) as f:
        fl = f.read()

    if "O   R   C   A" in fl:
        x, y, factor= orca(fl)
        # factor=1
    else:
        x, y = gaussian(fl)
        factor=1

    x = np.array(x)*factor
    plt.plot(x, y, color=args.color, alpha=0.4)
    plt.scatter(x, y, color=args.color)

    x = list(x)
    text_on_graphs(x, y, prod_left=args.prod_right, arrow=args.arrow)
    show_graph(y)


def orca(fl):
    fac = float([i for i in fl.splitlines() if "Length of steepest" in i][0].split()[6])
    splitter = ' IRC PATH SUMMARY '
    fl = fl.split(splitter)[1]
    fl = fl.split('Timings for individual modules')[0]
    data = fl.splitlines()[5:-1]
    x = [float(i.split()[0]) for i in data if (i and i.split())]
    y = [float(i.split()[2]) for i in data if (i and i.split())]
    idx = y.index(max(y))
    x_max = x[idx]
    x = list(np.array(x) - x_max)
    
    return x,y,fac

def gaussian(fl):
    splitter = 'Summary of reaction path following'
    fl = fl.split(splitter)[1]
    fl = fl.split('--------------------------------------------------------------------------')[1]
    data = fl.split('\n')[2:]
    x = [float(i.split()[2]) for i in data if (i and i.split())]
    y = [float(i.split()[1]) for i in data if (i and i.split())]
    y = (np.array(y) - max(y))*627.51
    return x,y 


def text_on_graphs(x, y, xpad=0.2, ypad=2, idx=3, prod_left=True, arrow=False):
    idx_ = y.index(max(y))
    l, r = x[:idx_], x[idx_:]
    yl1, yl2 = y[x.index(l[idx])] + ypad, y[x.index(l[-idx])] + ypad
    yr1, yr2 = y[x.index(r[idx])] + ypad, y[x.index(r[-idx])] + ypad

    if arrow:
        plt.arrow(l[-idx]-xpad, yl2, l[idx]-l[-idx], yl1-yl2, head_width=0.09)
        plt.arrow(r[idx]+xpad, yr1, r[-idx]-r[idx], yr2-yr1, head_width=0.09)

    left_text = 'To products' if prod_left else 'To reagents'
    right_text = 'To reagents' if prod_left else 'To products'

    plt.text(x=l[len(l)//2], y=max(y)-2+args.y_pad, s=left_text, horizontalalignment='center', verticalalignment='center')
    plt.text(x=r[len(r)//2], y=max(y)-2+args.y_pad, s=right_text, horizontalalignment='center', verticalalignment='center')

    plt.vlines(0, min(x)-20, max(x)+10, linestyles='dashed', alpha=0.4, color="gray")


def show_graph(y):
    plt.title(args.title)
    plt.xlabel('Reaction coordinates')
    plt.ylabel('$\Delta$E [kcal/mol]')
    plt.ylim((min(y)-3, max(y)+3))

    plt.tight_layout()
    fig = plt.gcf()
    
    if args.save:
        save_graph(fig)
    try:
        plt.show()
    except Exception:
        print('It was impossible to visulaize the graph of these file because you\'re running a non interactive terminal session. Graph is dumped in pickle and png format.')
        if not args.save:
            save_graph(fig)

def save_graph(fig):
    directory = args.graph_directory

    if os.path.exists(directory):
        if 'y' in input(f'A directory named {directory} already exists. Existing directory  will be deleted, wanna procede? [y/n]').lower():
            shutil.rmtree(directory)   
        else:
            sys.exit()
    os.mkdir(directory)
  
    with open(os.path.join(directory, 'irc.pickle'), 'wb') as f:
        pickle.dump(fig, f)
    plt.savefig(os.path.join(directory, 'irc.png'), dpi=700)


if __name__=='__main__':
    for i in args.file:
        get_irc(i)