

import argparse 
import matplotlib.pyplot as plt
import time

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file', help='ORCA output file')
parser.add_argument('-s', '--starting_index', default=0, type=int)
parser.add_argument('-e', '--ending_index', default=-1, type=int)
parser.add_argument('-o', '--output', type=str, default='graphs', help='Specify the filename (w\out the extention) in which save the output graph')
parser.add_argument('-sl', '--silent', action='store_true', default=False)

args = parser.parse_args()


GRAPHS = ('Energy change', 'RMS gradient', 'MAX gradient', 'RMS step', 'MAX step')

def get_energies(file):
    with open(file) as f:
        lines = [i for i in f.readlines()]
    flag = False
    conv = []
    for i in range(len(lines)):
        if "|Geometry convergence|" in lines[i]:
            flag = True
        if flag:
            a = [lines[i+j] for j in range(2, 8) if not ('--' in lines[i+j] or '.....' in lines[i+j])]
            if len(a) == 4:
                a.insert(0, '          Energy change       0.0000000000            0.0000000000      NO\n')
            conv.append(a)
            flag = False

    print(f'============ LAST ITERATION ({len(conv)-1}) ============\n')
    print(''.join(conv[-1]))
    print('===============================================')

    tmp_dict = {idx: [] for idx, _ in enumerate(conv[0])}
    tmp_dict['it'] = [i for i in range(len(conv))]
    for i in conv:
        for idx, j in enumerate(i):
            tmp_dict[idx].append(float(j.strip().split()[2]))

    return tmp_dict


def get_graph(conv: dict):
    fig, ax = plt.subplots(3, 2, sharex=True)

    iterations = conv.pop('it')

    for idx, i in enumerate(GRAPHS, start=1):
        if idx == 1:
            continue
        ax.flat[idx].set(title=i)
        ax.flat[idx].plot(iterations[args.starting_index : args.ending_index], conv[idx-1][args.starting_index: args.ending_index])
    for i in ax[0, :]:
        i.remove()
    gs = ax[0, 0].get_gridspec()
    axbig = fig.add_subplot(gs[0, :])
    if min(conv[0][args.starting_index: args.ending_index]) < 0:
        axbig.hlines(0, iterations[args.starting_index], iterations[args.ending_index], color='r', alpha=0.2, linestyle='dashed')
    axbig.set(title=GRAPHS[0])
    axbig.plot(iterations[args.starting_index : args.ending_index], conv[0][args.starting_index: args.ending_index])


    plt.tight_layout()
    if not args.silent:
        plt.savefig(f'{args.output if args.output == "graphs" else args.output}.png', dpi=300)
    else:
        plt.show()
    



if __name__=='__main__':
    conv = get_energies(args.file)
    get_graph(conv)