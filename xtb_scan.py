import argparse
import re

import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument('file', help='Scan log file with 1D scan', nargs='+')
parser.add_argument('start', help='The initial point of the scan', type=float)
parser.add_argument('end', help='The final point of the scan', type=float)
parser.add_argument('step', help='Insert the number of step you\'ve scanned', type=int)

parser.add_argument('--output', default=' ')

args = parser.parse_args()

def get_energy(filename):
    with open(filename, 'r') as f:
        fl = f.read()
    
    first = fl.split('\n')[0]
    regex = r'energy:'
    structurs = re.split(regex, fl, maxsplit=args.step)
    data = []
    for i in structurs[1:]:
        energy = i.split('\n')[0]
        data.append((energy.split()[0]))
    return data


def write_file(energies, idx):
    if args.output != ' ':
        out = args.output+str(idx)+'.dat' if len(args.file) > 1 else args.output+'.dat'
    else:
        out = f'scan_{args.file[idx].split("/")[0]}.dat'
    with open(out, 'w') as f:
        for x, y in zip(np.linspace(args.start, args.end, args.step), energies):
            f.write(f'{x}   {y}\n')
    return None



if __name__=='__main__':
    for idx, i in enumerate(args.file):
        en = get_energy(i)
        write_file(en, idx)
    # show_graph(colors, not_colors)