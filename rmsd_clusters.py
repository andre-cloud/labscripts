from functions.io import parse_ensemble
from functions.chem import get_positions
import argparse
# from rmsd import rmsd
import numpy as np 
from itertools import product


parser = argparse.ArgumentParser()

parser.add_argument('files', nargs='+')
# parser.add_argument('thr', type=float)

args = parser.parse_args()


points = {file: parse_ensemble(file) for file in args.files}
geoms = {}
matches = {}



def rmsd(A, B):
    Coord = len(A[0])
    NAtom = len(A)
    cum = 0.0
    for i in range(NAtom):
        for j in range(Coord):
            cum += (A[i][j] - B[i][j])**2.0
    return np.sqrt(cum / NAtom)


def centroid(A):
    A = A.mean(axis=0)
    return A


def kabsch(coord_var, coord_ref):

    # covariance matrix
    covar = np.dot(coord_var.T, coord_ref)
    # SVD  http://en.wikipedia.org/wiki/Kabsch_algorithm
    v, s, wt = np.linalg.svd(covar)
    # proper/improper rotation, JCC 2004, 25, 1894.
    d = (np.linalg.det(v) * np.linalg.det(wt)) < 0.0

    if d: # antialigns of the last singular vector
        s[-1] = -s[-1]
        v[:, -1] = -v[:, -1]
    R = np.dot(v, wt)


    #####
    # use np.dot(coord_var_cen,R) + trans

    return R



def calc_rmsd(v, w):
    trans = centroid(v)
    coord_var_cen = w - centroid(w)
    coord_ref_cen = v - centroid(v)
    # 4. Generate rotation matrix by Kabsch algorithm
    R = kabsch(coord_var_cen, coord_ref_cen)

    # 5. Rotate and translate
    coord_var_shifted = np.dot(coord_var_cen,R) + trans

    # 7. Std-out initial and final RMSD
    return rmsd(coord_var_shifted, v)


for i in points:
    for idx, g in enumerate(points[i]):
        geoms[idx] = get_positions(g)
    points[i] = geoms.copy()
    geoms = {}

del geoms

lists = [[f'{i}?{j}' for j in points[i]] for i in points]

res = list(product(*lists))
for g1, g2 in res:
    d1, idx1 = g1.split('?')
    d2, idx2 = g2.split('?')
    v, w = points[d1][int(idx1)], points[d2][int(idx2)]

    rmsd2 = calc_rmsd(v, w)    

    if g1.replace('?', '_') not in matches: 
        matches[g1.replace('?', '_')] = [rmsd2] 
    else: 
        matches[g1.replace('?', '_')].append(rmsd2)


for i in matches:
    print(matches[i])