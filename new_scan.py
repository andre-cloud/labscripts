import argparse
from collections import Counter
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
parser.add_argument('-xpad', help='Define the x shift of the text for the labels of the peaks. Defaults %(default)s', default=2, type=float)
parser.add_argument('-ypad', help='Define the y shift of the text for the labels of the peaks. Defaults %(default)s', default=0, type=float)
parser.add_argument('-tick', '--tick_height', help='Set the height of the ticks for the maximum values. Default %(default)s', default=1.5, type=float)
parser.add_argument('-nl', '--no_legend', help='Disable the legend on the graph', action='store_true')
parser.add_argument('-c', '--colors', nargs='+', help='Declaire the colors of the graph')
parser.add_argument('-xl', '--x_label', help='Set the x label for the graph. Default %(default)s', default='Dihedral angle')
parser.add_argument('-ff', '--from_file', help='Get the data from an xy file. For orca calucation use .relaxscanact.dat file',action='store_true')
parser.add_argument('--no_label', help='Don\'t write the label in the peak', action='store_true')

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
    x, y = check_x(x, y)
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
        back_2, fowr_2 = idx-2, idx+2
        if back<0: back_2=-2; back=-1
        if back_2<0 and back>=0: back_2=-1
        if fowr>=len(x): fowr=0; fowr_2=1
        if fowr_2==len(x): fowr_2=0 
        if y[idx] > y[back] and y[idx] > y[fowr] and y[idx] > args.threshold and y[idx] > y[back_2] and y[idx] > y[fowr_2]: maxs.append(i)
    return maxs


def write_max(x,y, pad, xpad, ypad):
    if not args.no_label:
        for idx, i in enumerate(x):
            plt.vlines(x=i, ymin=0, ymax=1, color='gray', alpha=0.6)
            plt.vlines(x=i, ymin=y[idx]-args.tick_height/2+pad, ymax=y[idx]+args.tick_height/2+pad, color='gray', alpha=0.6)
            ha = 'left' if xpad>=0 else 'right'
            plt.text(x=i+xpad, y=y[idx]+ypad+pad, s=f'∆E = {float(y[idx]):.2f}\n{float(i):.2f}', color='gray', horizontalalignment=ha, verticalalignment='center',)


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

def check_x(x, y):
    x_ = list(x)
    for idx, i in enumerate(x_):
        if i<-180 or i>180: 
            new_x = (360-abs(i)) * (1 if i<-180 else -1)
            x_[idx] = new_x
    dict_ = {k:v for k,v in zip(x_, y)}
    dict_ = {k: v for k, v in sorted(list(dict_.items()))}
    round_x = np.round(np.array(list(dict_.keys())), 3)
    if len(round_x) != len(set(round_x)):
        keys=[]
        dup = [item for item, count in Counter(round_x).items() if count > 1]
        for i in dup:
            for j in dict_:
                if i-j < 0.002:
                    keys.append(j)
            min_value=min([dict_[key] for key in keys])
            dict_[keys[0]] = min_value
            dict_.pop(keys[1])
            keys=[]

    return np.array(list(dict_.keys())), np.array(list(dict_.values()))


def from_file(filename, i):
    padding = i*args.padding

    with open(filename) as f:
        fl = f.read()

    data = [[float(i.strip().split()[0]), float(i.strip().split()[1])] for i in fl.split('\n') if i]
    data = np.array(data)
    x, y = data[:, 0], data[:, 1]
    y -= min(y)
    y *= 627.51
    x, y = check_x(x, y)
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




if __name__=='__main__':
    with alive_bar(len(args.file), title='Parsing files') as bar:
        for idx, file in enumerate(args.file):
            if not args.from_file:
                get_scan(file, idx)
            else:
                from_file(file, idx)
            bar()

            
    show_graph()