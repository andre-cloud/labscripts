from functions.io import parse_ensemble
from functions.chem import get_positions
import argparse, sys
import numpy as np 
from itertools import product


"""

RMSD is not a good parameter to consider if two structure are the same conformers.

CREST if the ∆E is higher than 0.05 kcal/mol are two different conformers if 
rotational constant are similar. 

"""


parser = argparse.ArgumentParser()

parser.add_argument('files', nargs='+')
parser.add_argument('--rmsd_thr', type=float, default=1.0)
parser.add_argument('--energy_thr', type=float, default=0.05)

args = parser.parse_args()


points = {file: parse_ensemble(file) for file in args.files}
ens_points = {file: parse_ensemble(file) for file in args.files}
geoms = {}
ens = {}
matches = {}
ens_matches = {}



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
        ens[idx] = float(g.splitlines()[1])
        geoms[idx] = get_positions(g)
    points[i] = geoms.copy()
    ens_points[i] = ens.copy()
    geoms = {}
    ens = {}

del geoms
del ens




lists = [[f'{i}?{j}' for j in points[i]] for i in points]


res = list(product(*lists))
for g1, g2 in res:
    d1, idx1 = g1.split('?')
    d2, idx2 = g2.split('?')
    v, w = points[d1][int(idx1)], points[d2][int(idx2)]
    e_v, e_w = ens_points[d1][int(idx1)], ens_points[d2][int(idx2)]

    rmsd2 = calc_rmsd(v, w)    

    if g1.replace('?', '_') not in matches: 
        matches[g1.replace('?', '_')] = [rmsd2] 
        ens_matches[g1.replace('?', '_')] = [np.abs((e_v-e_w)*627.51)] 
    else: 
        matches[g1.replace('?', '_')].append(rmsd2)
        ens_matches[g1.replace('?', '_')].append(np.abs((e_v-e_w)*627.51))


    

ma_rmsd = np.empty((len(matches), len(matches[list(matches.keys())[0]])))
ma_ens = np.empty((len(matches), len(matches[list(matches.keys())[0]])))
for idx, i in enumerate(list(matches.keys())):
    ma_rmsd[idx] = matches[i]
    ma_ens[idx] = ens_matches[i]

    

c = 0
for idx, i in enumerate(ma_rmsd):
    if True in set(ma_ens[idx]<=args.energy_thr):
        c += 1

print(f'Geometries of {args.files[0]} in {args.files[1]}: {c} out of {ma_rmsd.shape[0]} ({c/ma_rmsd.shape[0]*100:.2f}%)')


c = 0
for idx, i in enumerate(ma_rmsd.T):
    if True in set(ma_ens.T[idx]<=args.energy_thr):
        c += 1

print(f'Geometries of {args.files[1]} in {args.files[0]}: {c} out of {ma_rmsd.shape[1]} ({c/ma_rmsd.shape[1]*100:.2f}%)\n')

s = ' Tensor Statistics '
print('-'*15+s+'-'*15)
print(f'Min RMSD: {np.min(ma_rmsd):.3f} - Max RMSD: {np.max(ma_rmsd):.3f} - Mean RMSD: {np.mean(ma_rmsd):.3f}')
print(f'Min ∆E: {np.min(ma_ens):.5f} - Max ∆E: {np.max(ma_ens):.5f} - Mean ∆E: {np.mean(ma_ens):.5f}')
print('-'*(15+len(s)+15))
print()