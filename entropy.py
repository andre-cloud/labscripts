import argparse

import numpy as np
from scipy.constants import R

parser = argparse.ArgumentParser()

parser.add_argument('file', nargs='+')
parser.add_argument('-cre', '--crest_ensemble', action='store_true')
parser.add_argument('-cen', '--censo_ensemble', action='store_true')
parser.add_argument('-T', '--temp', type=float, default=298.15)
parser.add_argument('-ltx', action='store_true')


args = parser.parse_args()


if not (args.crest_ensemble or args.censo_ensemble):
    raise Exception('You have to declaire which ensemble you\'re parsing')

def split_file(filename):
    with open(filename) as f:
        fl = f.read()
    fl = fl.splitlines()
    points = []
    prev_i = 0
    for i in range(0, len(fl)+1, int(fl[0])+2):
        if fl[prev_i:i]: points.append('\n'.join(fl[prev_i:i])) 
        prev_i=i

    if args.crest_ensemble:
        en = [float(i.splitlines()[1]) for i in points]
        return np.array(en)

    if args.censo_ensemble:
        en = [float(i.splitlines()[1].split()[1]) for i in points]
        return np.array(en)

def boltzmann_perc(en):
    return np.exp(-(en-np.min(en))/R*args.temp)/np.sum(np.exp(-(en-np.min(en))/R*args.temp))



if __name__ == '__main__':
    print('Filename\tGav(kcal/mol)\tN Conf\tSCRE(J/mol K)')
    for file in args.file:
        en = split_file(file)
        p = boltzmann_perc(en*627.51)
        entropy = R*np.sum(p*np.log10(p))
        mean_energy = np.sum(en*p)
        confs = len(list(en))
        if not args.ltx : print(f"{file}\t{mean_energy:.4f}\t{confs}\t{entropy:.3f}")
        else: print(f"{file}\t{mean_energy:.4f} & {entropy:.3f} ({confs})")
