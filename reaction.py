import argparse
import sys


import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import numpy as np
import pickle



Path = mpath.Path
fig, ax = plt.subplots()

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
        'labelpoint' : [None],
        'zero' : ['0'],
        'title' : ['Reaction Energy Path'], 
        'style' : ['-'],
        'linewidth' : ['1.5'],
        'points' : None,
        'pointscolor' : None,
        'pointsstyle' : None,
        'ypad' : ['0'],
        'horizontalalignment' : ['0'],
    }

    def __init__(self, filename):
        self.filename = filename
        self.read_file()

        print(self.graph_prm)


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
                horizontalalignment=self.graph_prm['horizontalalignment']
                )
            x_labels(
                self.labels[idx], 
                self.graph_prm['colors'][idx]
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
            title=self.graph_prm['title'][0], 
            data=data
            )


    

    def read_file(self):
        with open(self.filename) as f:
            self.file = f.read()

        pes, labels, graph_prm = [i.strip().split('\n')[1:] for i in self.file.split('--# ')[1:] if i]
        pes = [[float(j) for j in i.split(',') if j] for i in pes if i]
        self.labels = [[j for j in i.split(',')] for i in labels]
        self.pes = np.array(pes, dtype='float64')


        # self.graph_prm = {i.split(' : ')[0].lower():[j for j in i.split(' : ')[1].split(', ')] for i in graph_prm if i}
        for i in graph_prm:
            if i and not i.startswith('#'):
                self.graph_prm[i.split(' : ')[0].lower()] = [j.strip() for j in i.split(' : ')[1].split(';')]

        self.control_prm()

        assert len(self.pes) == len(self.labels)
        for i in range(len(self.pes)-1):
            assert len(self.pes[i]) == len(self.labels[i])


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




def display_points(data, color, linestyle):
    plt.scatter(data[0], data[1], color=color, alpha=0.6, linewidth=0.3)
    plt.hlines(data[1], data[0] - .1, data[0] + 0.1, color=color, alpha=0.7, linestyle=linestyle)


def create_path(data, color, label, labelpoint, zero, linestyle, linewidth:float, ypad, horizontalalignment):
    for j, d in enumerate(data):
        if j == len(data)-1:
            break
        if j == 0:
            path_patch = mpatches.PathPatch(
                Path([(j, data[j]), (j + 0.5, data[j]), (j + 0.5, data[j + 1]),
                        (j + 1, data[j + 1])],
                        [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]),
                fc="none", transform=ax.transData, color=color, label=label, linestyle=linestyle, linewidth=linewidth)
        else:
            path_patch = mpatches.PathPatch(
                Path([(j, data[j]), (j + 0.5, data[j]), (j + 0.5, data[j + 1]),
                        (j + 1, data[j + 1])],
                        [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]),
                fc="none", transform=ax.transData, color=color, label=label, linestyle=linestyle, linewidth=linewidth)

        ax.add_patch(path_patch)
        plt.hlines(data[j], j - 0.1, j + 0.1, color=color, alpha=0.7, linestyle=linestyle)
    plt.hlines(data[-1], len(data) - 1.1, len(data) -0.9, color=color, alpha=0.7, linestyle=linestyle)
    leg_text.append(path_patch)
    for x, y in enumerate(data):
        plt.scatter(x, y, color=color, alpha=0.6, linewidth=0.3)

    if labelpoint is not None:
        for idx, i in enumerate(labelpoint):
            va = 'center'
            ha = horizontalalignment[idx] if i != len(data)-1 else 'right'
            x = i+0.2 if ha != 'right' else i-.2
            plt.text(x,data[i]+float(ypad), f"$\Delta G^‡$ = {np.round(data[i]-data[zero], 2)} kcal/mol", verticalalignment=va, horizontalalignment=ha, color=color)


def x_labels(x_label, color):
    plt.xticks(range(len(x_label)), x_label)

    newax = []
    for i in range(len(ax_label)):
        if i > 0:
            y = ax.twiny()
            newax.append(y)
    for i in range(len(newax)):
        newax[i].set_xticks(locs)
        newax[i].set_xlim(ax.get_xlim())
        newax[i].tick_params(axis='x', color=color)
        newax[i].set_xticklabels(newax_text_list[i + 1])
        newax_text_list.append(newax[i])
        newax[i].xaxis.set_ticks_position('bottom')
        newax[i].xaxis.set_label_position('bottom')
        newax[i].xaxis.set_ticks_position('none')
        newax[i].spines['bottom'].set_position(('outward', 15 * (i + 1)))
        newax[i].spines['bottom'].set_visible(False)


def show_graph(legend, title, data, save=False):
    ax.set_title(title)
    ax.set_ylabel(r"$G_{rel}$ (kcal / mol)")
    plt.minorticks_on()

    ax.tick_params(axis='x', which='minor', bottom=False)
    ax.tick_params(which='minor', labelright=True, right=True)
    ax.tick_params(labelright=True, right=True)
    
    
    # frame1 = plt.gca()
    # frame1.axes.get_yaxis().set_ticks([])

    # plt.hlines(data[1], 1, 3+.3, linestyles='--', color='salmon', alpha=0.7)
    # plt.hlines(data[5], 5, 3-.3, linestyles='--', color='salmon', alpha=0.7)
    # plt.vlines(3, data[1], data[5], linestyles='--', color='salmon', alpha=0.7)
    # plt.text(3, np.mean((data[1], data[5])), s='$\Delta\Delta G^‡$',verticalalignment='center', horizontalalignment='center', bbox=dict(facecolor='white', edgecolor='salmon', boxstyle='round,pad=1'))
    
    
    if legend:
        plt.legend(handles=leg_text)
    plt.tight_layout()
    if save:
        with open('tests/Rxn_profile.pickle', 'wb') as f:
            pickle.dump(fig, f)
        plt.savefig('tests/Rxn_profile.png', dpi=700)
    plt.show()

ax_label = []

locs, labels = plt.xticks()


prms = '\n'.join([f"- {i} : " for i in list(Interpreter.graph_prm.keys())])
parser = argparse.ArgumentParser(description=f'''
Input file MUST contains, in this order, and in this formatting:

--# PES
[...]
--# LABELS
[...]
--# GRAPHS
[...]

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

