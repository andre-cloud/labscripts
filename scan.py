import argparse
from collections import Counter
import os
import random
import sys
import shutil

from alive_progress import alive_bar
import cclib
import matplotlib as mpl 
import matplotlib.pyplot as plt
import numpy as np
import pickle


graphs={}

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
parser.add_argument('-a', '--absolute', action='store_true', help='Toggle if absolute energy is desired')
parser.add_argument('-r', '--random', help='Activate random color generation', action='store_true')

parser.add_argument('--save', help='Save pickle and csvs of the graph', action='store_true')
parser.add_argument('-gd','--graph_directory', help='Define the directory in which you want to save the files of the graph. Default: %(default)s', default='scan_graph')


args = parser.parse_args()

plots = []
lines = []
labels = []

# mpl.rc('font', size=12)

def get_scan(filename, i):
    padding = i*args.padding
    data = cclib.io.ccread(filename)
    x = data.scanparm[0]
    if not args.absolute : 
        y = list((np.array(data.scanenergies)- min(np.array(data.scanenergies)))*23.060541945329)
    else:
        y = np.array(data.scanenergies)*0.03674930495120813
    dict_sort = {}
    x, y = check_x(x, y)
    maxs = find_max(x, y)
    y_max = np.array([y[list(x).index(i)] for i in maxs])

    write_max(maxs, y_max, pad=padding, xpad=args.xpad, ypad=args.ypad)

    plot(x, y, padding, filename, maxs, y_max)

    # plt.plot(x,y+padding, '--', alpha=0.4, )


def find_max(x, y):
    m = min(y)
    maxs = []
    for idx, i in enumerate(x):
        back, fowr = idx-1, idx+1
        back_2, fowr_2 = idx-2, idx+2
        if back<0: back_2=-2; back=-1
        if back_2<0 and back>=0: back_2=-1
        if fowr>=len(x): fowr=0; fowr_2=1
        if fowr_2==len(x): fowr_2=0 
        if y[idx] > y[back] and y[idx] > y[fowr] and y[idx] > m+args.threshold and y[idx] > y[back_2] and y[idx] > y[fowr_2]: maxs.append(i)
    return maxs


def write_max(x,y, pad, xpad, ypad):
    if not args.no_label:
        for idx, i in enumerate(x):
            plt.vlines(x=i, ymin=0, ymax=1, color='gray', alpha=0.6)
            plt.vlines(x=i, ymin=y[idx]-args.tick_height/2+pad, ymax=y[idx]+args.tick_height/2+pad, color='gray', alpha=0.6)
            ha = 'left' if xpad>=0 else 'right'
            plt.text(x=i+xpad, y=y[idx]+ypad+pad, s=f'∆E = {float(y[idx]):.2f}\n{float(i):.2f}', color='gray', horizontalalignment=ha, verticalalignment='center',)


def show_graph():
    ax = plt.gca()

    if not args.no_legend: legend = ax.legend(
        loc='upper center', bbox_to_anchor=(0.5, -.125), fancybox=True, shadow=True, ncol=3
    )
    if not args.absolute: plt.ylim(bottom=0)
    plt.title(args.title)
    plt.ylabel('Energy $\Delta E\quad[kcal/mol]$')
    plt.xlabel(args.x_label)
    plt.tight_layout()

    # fig = plt.gcf().set_dpi(1000)
    if args.save:
        save_graph(fig)
    
    if not args.no_legend:
        for obj in legend.__dict__['legendHandles']:
            obj.set_picker(5)

        def on_pick(event):
            flag = False
            groups = {key: [val, lin] for key, val, lin in zip(labels, plots, lines)}
            label = event.artist.get_label()
            for key in groups:
                if key == label:
                    if not groups[key][0].__dict__['_alpha'] or groups[key][0].__dict__['_alpha'] != 1:
                        groups[key][0].set_alpha(1.0)
                        groups[key][1][0].set_alpha(0.4)
                    else:
                        flag = True
                else:
                    groups[key][0].set_alpha(0.1)
                    groups[key][1][0].set_alpha(0)
            if flag:
                for key in groups:
                    groups[key][0].set_alpha(None)
                    groups[key][1][0].set_alpha(0.4)
                flag = False

            fig.canvas.draw()
    
        plt.connect('pick_event', on_pick)
    
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
    x, y = data[:, 0][::2], data[:, 1][::2]
    if not args.absolute : y -= min(y)
    if not args.absolute : y *= 627.51
    x, y = check_x(x, y)
    maxs = find_max(x, y)
    y_max = np.array([y[list(x).index(i)] for i in maxs])

    write_max(maxs, y_max, pad=padding, xpad=args.xpad, ypad=args.ypad)
    plot(x, y, padding, filename, maxs, y_max)


def plot(x, y, padding, filename, maxs, y_max):
    if not args.colors:
        if args.random:
            r = lambda: random.randint(0,255)
            c = '#%02X%02X%02X' % (r(),r(),r())
            l = plt.plot(x,y+padding, linewidth=0.5, alpha=0.4, color=c)
            a = plt.scatter(x,y+padding, label=os.path.split(filename)[1].strip('.log') if len(args.file)>1 else None, color=c)
            plt.scatter(maxs, y_max+padding, alpha=0.4, color='red')
            plots.append(a)
            lines.append(l)
            labels.append(os.path.split(filename)[1].strip('.log') if len(args.file)>1 else None)
            return
        l = plt.plot(x,y+padding, linewidth=0.5, alpha=0.4)
        a = plt.scatter(x,y+padding, label=os.path.split(filename)[1].strip('.log') if len(args.file)>1 else None)
        plt.scatter(maxs, y_max+padding, alpha=0.4, color='red')
        plots.append(a)
        lines.append(l)
        labels.append(os.path.split(filename)[1].strip('.log') if len(args.file)>1 else None)
        return

    l = plt.plot(x,y+padding, linewidth=1, alpha=0.4, color=args.colors[args.file.index(filename)])
    a = plt.scatter(x,y+padding, label=os.path.split(filename)[1].strip('.log') if len(args.file)>1 else None, color=args.colors[args.file.index(filename)])
    plt.scatter(maxs, y_max+padding, alpha=0.4, color='red')
    plots.append(a)
    lines.append(l)
    labels.append(os.path.split(filename)[1].strip('.log') if len(args.file)>1 else None)
    return

if __name__=='__main__':

    with alive_bar(len(args.file), title='Parsing files') as bar:
        for idx, file in enumerate(args.file):
            if not args.from_file:
                get_scan(file, idx)
            else:
                from_file(file, idx)
            bar()

            
    show_graph()
