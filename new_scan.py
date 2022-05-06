import argparse
import os
import sys
import shutil

from alive_progress import alive_bar
import cclib
import matplotlib.pyplot as plt
import numpy as np
import pickle

parser = argparse.ArgumentParser()

parser.add_argument('file', help='Scan log file(s)', nargs='+')
parser.add_argument('-thr', '--threshold', help='Energy value below which no maximum can be found (in kcal/mol). Default: %(default)s', default=5, type=float)
parser.add_argument('-tit', '--title', help='define the title of the graph. Default %(default)s', default='Scan')
parser.add_argument('-p', '--padding', help='Define the padding of two or more graphs. Default %(default)s', default=20, type=int)
parser.add_argument('-xpad', help='Define the x shift of the text for the labels of the peaks. Defaults %(default)s', default=2, type=int)
parser.add_argument('-ypad', help='Define the y shift of the text for the labels of the peaks. Defaults %(default)s', default=0, type=int)
parser.add_argument('-tick', '--tick_height', help='Set the height of the ticks for the maximum values. Default %(default)s', default=1.5, type=float)
parser.add_argument('-nl', '--no_legend', help='Disable the legend on the graph', action='store_true')
parser.add_argument('-c', '--colors', nargs='+', help='Declaire the colors of the graph')
parser.add_argument('-xl', '--x_label', help='Set the x label for the graph. Default %(default)s', default='Dihedral angle')

parser.add_argument('--save', help='Save pickle and csvs of the graph', action='store_true')
parser.add_argument('-gd','--graph_directory', help='Define the directory in which you want to save the files of the graph. Default: %(default)s', default='scan_graph')


args = parser.parse_args()


def get_scan(filename, i):
    padding = i*args.padding
    data = cclib.io.ccread(filename)
    x = data.scanparm[0]
    y = list((np.array(data.scanenergies)- min(np.array(data.scanenergies)))*23.060541945329)
    dict_sort = {}
    for i in sorted(list(x)):
        if i not in dict_sort:
            dict_sort[i] = y[list(x).index(i)]

    x = np.array(list(dict_sort.keys()))
    y = np.array(list(dict_sort.values()))
    maxs = find_max(x, y)
    y_max = np.array([y[list(x).index(i)] for i in maxs])

    write_max(maxs, y_max, pad=padding, xpad=args.xpad, ypad=args.ypad)
    if not args.colors:
        plt.plot(x,y+padding, linewidth=0.5, alpha=0.4)
        plt.scatter(x,y+padding, label=os.path.split(filename)[1].strip('.log') if len(args.file)>1 else None)
        plt.scatter(maxs, y_max+padding, alpha=0.4, color='red')
    else: 
        plt.plot(x,y+padding, linewidth=0.5, alpha=0.4, color=args.colors[args.file.index(filename)])
        plt.scatter(x,y+padding, label=os.path.split(filename)[1].strip('.log') if len(args.file)>1 else None, color=args.colors[args.file.index(filename)])
        plt.scatter(maxs, y_max+padding, alpha=0.4, color='red')

    # plt.plot(x,y+padding, '--', alpha=0.4, )

    


def find_max(x, y):
    maxs = []
    for idx, i in enumerate(x):
        back, fowr = idx-1, idx+1
        if back<0: back = -1
        if fowr>=len(x): fowr=0
        if y[idx] > y[back] and y[idx] > y[fowr] and y[idx] > args.threshold: maxs.append(i)
    return maxs


def write_max(x,y, pad, xpad, ypad):
    for idx, i in enumerate(x):
        plt.vlines(x=i, ymin=0, ymax=1, color='gray', alpha=0.6)
        plt.vlines(x=i, ymin=y[idx]-args.tick_height/2+pad, ymax=y[idx]+args.tick_height/2+pad, color='gray', alpha=0.6)
        ha = 'left' if xpad>=0 else 'right'
        plt.text(x=i+xpad, y=y[idx]+ypad+pad, s=f'∆E = {round(float(y[idx]), 2)}\n{round(float(i), 2)}', color='gray', horizontalalignment=ha, verticalalignment='center',)


def show_graph():
    if not args.no_legend: plt.legend(
        loc='upper center', bbox_to_anchor=(0.5, -.125), fancybox=True, shadow=True, ncol=3
    )
    plt.ylim(bottom=0)
    plt.title(args.title)
    plt.ylabel('Energy $\Delta E\quad[kcal/mol]$')
    plt.xlabel(args.x_label)
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

    ln = 2

    with alive_bar(ln, title='Saving plot') as bar:      
        with open(os.path.join(directory, 'scan.pickle'), 'wb') as f:
            pickle.dump(fig, f)
        bar()
        plt.savefig(os.path.join(directory, 'scan.png'), dpi=700)
        bar()




if __name__=='__main__':
    with alive_bar(len(args.file), title='Parsing files') as bar:
        for idx, file in enumerate(args.file):
            get_scan(file, idx)
            bar()

    show_graph()