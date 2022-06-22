import argparse
import os, sys

import numpy as np 



parser = argparse.ArgumentParser()

parser.add_argument('file', nargs='+')
parser.add_argument('-i', '--index', help='Define the index (starting at 1) of the atoms', nargs='+', required=True)


args = parser.parse_args()


def dihedral(p) -> float: 
    """
    p : (3, 4) size numpy array of positions

    return : float
    """
    p0 = p[0]
    p1 = p[1]
    p2 = p[2]
    p3 = p[3]

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    b1 /= np.linalg.norm(b1)

    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


def get_position(file, indexes):
    fl = file.split('\n')[2:]
    return np.array([fl[int(i)-1].split()[1:] for i in indexes], dtype=np.float64)


def read_input(filename):
    """
    filename : str

    return : array with n geometrys
    """
    with open(filename) as f:
        fl = f.read()
    fl = fl.splitlines()
    points = []
    prev_i = 0
    for i in range(0, len(fl), int(fl[0])+2):
        if fl[prev_i:i]: points.append('\n'.join(fl[prev_i:i])) 
        prev_i=i
    return points


def sign(indexes, geometry):
    atoms = get_position(geometry, indexes)
    dh = dihedral(atoms)
    print(np.sign(dh))
    return np.sign(dh)


if __name__=='__main__':
    for file in args.file:
        geoms = read_input(file)
        p, n = [], []
        for idx, geom in enumerate(geoms):
            if sign(args.index, geom) > 0:
                p.append(idx)
            else:
                n.append(idx)
        with open(f'{file}_p.xyz', 'w') as f:
            f.write('\n'.join([geoms[i] for i in p]))
        with open(f'{file}_n.xyz', 'w') as f:
            f.write('\n'.join([geoms[i] for i in n]))





