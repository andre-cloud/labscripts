import argparse
import os
from functions import io

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import numpy as np
import pickle



Path = mpath.Path
fig, ax = plt.subplots()
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']="\\usepackage{amsmath}"


newax_text_list = []
leg_text = []











class InputError(Exception):

    def __init__(self, prm, len):
        self.prm = prm
        self.len = len

    def __str__(self):
        return f"The parameter {self.prm} and the number of the reactions ({self.len}) are not consitents. Check the input file"


class Interpreter:
    graph_prm = {
        'label' : [],
        'colors' : [None],
        'legend' : ['false'],
        'save' : ['false'],
        'relative' : ['false'],
        'labelpoint' : [None],
        'zero' : ['0'],
        'title' : ['Reaction Energy Path'], 
        'style' : ['-'],
        'linewidth' : ['1.5'],
        'points' : [],
        'pointscolor' : [],
        'pointsstyle' : [],
        'ypad' : ['0'],
        'horizontalalignment' : ['left'],
        'y_limit' : [None],
        'decimal_digits' : ['1'],
        'directory' : ['reaction']
    }

    exp_prm = {
        'label' : 'Str: Legend label for each path',
        'colors' : 'Str: Matplotlib colors (so either CSS or HEX codes)',
        'legend' : 'Bool: to display legend',
        'save' : 'Bool: to save png and pickle files',
        'relative' : 'Bool: to indicate that energies in PES section are relative (i+(i-1))',
        'labelpoint' : 'List[Int]: starting from 0, to indicate near which point is asked to show the ∆G‡ relative to the zero',
        'zero' : 'List[Int]: starting from zero to indicate the zero of each path',
        'title' : 'Str: Title of the graph', 
        'style' : 'List[Str]: matplotlib stiles for plots',
        'linewidth' : 'List[Float]: linewidth of the plot\'s stroke',
        'points' : 'List[Tuple]: (x,y) position of points',
        'pointscolor' : 'List[Str]: similar to color but for points',
        'pointsstyle' : 'List[Str]: similar to style but for points',
        'ypad' : 'List[Float]: y-direction pad for the point label',
        'horizontalalignment' : '',
        'y_limit' : 'Float: define the top y limit',
        'decimal_digits' : 'Int: define the number of decimals digits for the labeled points', 
        'directory' : 'Str: name of the direcotry in which to save all files'
    }

    def __init__(self, filename):
        self.filename = filename
        self.read_file()


        for idx, data in enumerate(list(self.pes)):
            create_path(
                data=data,
                color=self.graph_prm['colors'][idx], 
                label=self.graph_prm['label'][idx], 
                labelpoint=np.array(self.graph_prm['labelpoint'], dtype="int8") if self.graph_prm['labelpoint'] != [None] else None, 
                zero=int(self.graph_prm['zero'][idx]),
                linestyle=self.graph_prm['style'][idx],
                linewidth=self.graph_prm['linewidth'][idx],
                ypad=self.graph_prm['ypad'][idx],
                horizontalalignment=self.graph_prm['horizontalalignment'],
                decimals = self.graph_prm['decimal_digits'][0]
                )
        
        x_labels(
            self.labels[0], 
            )
        if self.graph_prm['points']:
            for idx, i in enumerate(self.graph_prm['points']):
                data = np.array(i.strip('(').strip(')').split(','), dtype=np.float64)
                display_points(
                data=data,
                color=self.graph_prm['pointscolor'][idx], 
                linestyle=self.graph_prm['pointsstyle'][idx],
                )
                

        show_graph(
            legend=True if self.graph_prm['legend'][0].lower() == 'true' else False,
            save=True if self.graph_prm['save'][0].lower() == 'true' else False, 
            directory=self.graph_prm['directory'][0],
            title=self.graph_prm['title'][0], 
            data=data, 
            y_lim=self.graph_prm['y_limit'][0],
            file=self.file
            )


    

    def read_file(self):
        with open(self.filename) as f:
            self.file = f.read()

        pes, labels, graph_prm = [i.strip().split('\n')[1:] for i in self.file.split('--# ')[1:] if i]
        pes = [[float(j) for j in i.split(',') if j] for i in pes if i]
        self.labels = [[j.strip() for j in i.split(';')] for i in labels]
        print(self.labels)
        self.pes = np.array(pes, dtype='float64')


        # self.graph_prm = {i.split(' : ')[0].lower():[j for j in i.split(' : ')[1].split(', ')] for i in graph_prm if i}
        for i in graph_prm:
            if i and not i.startswith('#'):
                self.graph_prm[i.split(' : ')[0].lower()] = [j.strip() for j in i.split(' : ')[1].split(';')]

        self.control_prm()

    def __repr__(self):
        return self.pes, self.labels, self.graph_prm


    def control_prm(self):
        if len(self.graph_prm['colors']) !=  len(self.pes): self.graph_prm['colors'] = [None] * len(self.pes)

        if len(self.graph_prm['label']) !=  len(self.pes): self.graph_prm['label'] = [None] * len(self.pes)
        if set(self.graph_prm['label']) == set([None]): self.graph_prm['legend'] = 'False'
        

        if len(self.graph_prm['zero']) !=  len(self.pes):
            if set(self.graph_prm['zero']) != set(['0']):
                raise InputError('zero', len(self.pes))
            self.graph_prm['zero'] = ['0'] * len(self.pes)
        
        if len(self.graph_prm['style']) !=  len(self.pes): self.graph_prm['style'] = ['-'] * len(self.pes)

        if len(self.graph_prm['linewidth']) !=  len(self.pes): self.graph_prm['linewidth'] = ['1.5'] * len(self.pes)
        
        if len(self.graph_prm['ypad']) !=  len(self.pes): self.graph_prm['ypad'] = ['0'] * len(self.pes)
        if len(self.graph_prm['horizontalalignment']) !=  len(self.graph_prm['labelpoint']): self.graph_prm['horizontalalignment'] = ['left'] * len(self.graph_prm['labelpoint'])

        if self.graph_prm['points']:
            if len(self.graph_prm['pointscolor']) !=  len(self.pes): self.graph_prm['pointscolor'] = [None] * len(self.pes)
            if len(self.graph_prm['pointsstyle']) !=  len(self.pes): self.graph_prm['pointsstyle'] = ['-'] * len(self.pes)

        
        if self.graph_prm['relative'][0].lower() == 'true':
            for i in range(len(self.pes)):
                for j in range(len(self.pes[i])):
                    if j != 0:
                        self.pes[i][j] = self.pes[i][j] + self.pes[i][j-1]





def display_points(data, color, linestyle):
    plt.scatter(data[0], data[1], color=color, alpha=0.6, linewidth=0.3)
    plt.hlines(data[1], data[0] - .1, data[0] + 0.1, color=color, alpha=0.7, linestyle=linestyle)


def create_path(data, color, label, labelpoint, zero, linestyle, linewidth:float, ypad, horizontalalignment, decimals):
    for j, d in enumerate(data):
        if j == len(data)-1:
            break
        if j == 0:
            path_patch = mpatches.PathPatch(
                Path([(j, data[j]), (j + 0.5, data[j]), (j + 0.5, data[j + 1]),
                        (j + 1, data[j + 1])],
                        [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]),
                fc="none", transform=ax.transData, color=color, label=u'{0}'.format(label), linestyle=linestyle, linewidth=float(linewidth))
        else:
            path_patch = mpatches.PathPatch(
                Path([(j, data[j]), (j + 0.5, data[j]), (j + 0.5, data[j + 1]),
                        (j + 1, data[j + 1])],
                        [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]),
                fc="none", transform=ax.transData, color=color, label=u'{0}'.format(label), linestyle=linestyle, linewidth=float(linewidth))

        ax.add_patch(path_patch)
        plt.hlines(data[j], j - 0.1, j + 0.1, color=color, alpha=0.7, linestyle=linestyle)
    plt.hlines(data[-1], len(data) - 1.1, len(data) -0.9, color=color, alpha=0.7, linestyle=linestyle)
    leg_text.append(path_patch)
    for x, y in enumerate(data):
        plt.scatter(x, y, color=color, alpha=0.6, linewidth=0.3)

    if labelpoint is not None:
        print(labelpoint)
        for idx, i in enumerate(labelpoint):
            va = 'center'
            ha = horizontalalignment[idx] if i != len(data)-1 else 'right'
            x = i+0.2 if ha != 'right' else i-.2
            print(va, ha)
            plt.text(x,data[i]+float(ypad), f"$\Delta G^\ddagger$ = {np.round(data[i]-data[zero], int(decimals))} kcal/mol", verticalalignment=va, horizontalalignment=ha, color=color)


def x_labels(x_label):
    plt.xticks(range(len(x_label)), x_label)
    


def show_graph(legend, title, data, save, directory, y_lim:float, file):
    ax.set_title(title)
    ax.set_ylabel(r"$\Delta G$ (kcal / mol)")
    plt.minorticks_on()

    ax.tick_params(axis='x', which='minor', bottom=False)
    ax.tick_params(which='minor', labelright=True, right=True)
    ax.tick_params(labelright=True, right=True)

    if y_lim:
        ax.set_ylim(top=float(y_lim))
    
    if legend:
        plt.legend(handles=leg_text,  loc='upper center', bbox_to_anchor=(0.5, -.075), fancybox=True, shadow=True, ncol=4)
    plt.tight_layout()
    if save:
        io.mkdir(directory)
        with open(os.path.join(directory, 'Rxn_profile.pickle'), 'wb') as f:
            pickle.dump(fig, f)
        plt.savefig(os.path.join(directory, 'Rxn_profile.png'), dpi=700)
        with open(os.path.join(directory, 'input.txt'), 'w') as f:
            f.write(file)

    plt.show()

ax_label = []

locs, labels = plt.xticks()


prms = '\n'.join([f"- {i} : {Interpreter.exp_prm[i]}" for i in list(Interpreter.graph_prm.keys())])
parser = argparse.ArgumentParser(description=f'''
Input file MUST contains, in this order, and in this formatting:

--# PES
[...]
--# LABELS
[...]
--# GRAPHS
[...]

--- Explanations ---

PES: Each line is a "reaction path". A line can contains many energy comma separated.

LABELS: Contains the labels of the points of the PES Paths. MUST BE THE SAME NUMBER OF THE PESs LINES

GRAPHS: Contains all the parameters that can be tweched in the graph. All setting are listed here below
{prms}
''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('file', help='Input file')

args = parser.parse_args()






if __name__=='__main__':


    Interpreter(args.file)

    # for idx, i in enumerate(datas):
    #     create_path(i, colors[idx])
    #     x_labels(xaxis_text, labels)

    # show_graph(save=True)

