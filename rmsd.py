import numpy as np
import argparse
from periodictable import core, covalent_radius, mass
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument('file', help="Ensemble file in xyz coordinates")
parser.add_argument('output', help="Filename for the output")
parser.add_argument('-r', '--reference', help="Conformer for reference. It starts at 0. Default 0", default=0)

args = parser.parse_args()

pt = core.PeriodicTable(table="H=1")
covalent_radius.init(pt)
mass.init(pt)


def _parse_xyz_str(fl: str):
    """
    Parse an xyz geom descriptor

    :param fl: string of the file splitted in a single geometry
    :type fl: str

    :return: list of atoms, XYZ position of the atoms
    :rtype: tuple
    """
    en = float(fl[1])
    fl = fl[2:]
    atoms, geom = [], []
    for line in fl:
        a, *g = line.split()
        atoms.append(a)
        geom.append(g)
    return np.array(atoms), np.array(geom, dtype=float), en

def read_ensemble(file) -> list:
    """
    Read the initial ensemble and return the ensemble list
    Not only XYZ file is supported. OBABEL is required

    :param file: initial ensemble file
    :type file: str

    :return: whole ensemble list
    :rtype: list
    """

    confs = []
    ens = []

    with open(file) as f:
        fl = f.readlines()

    n_atoms = int(fl[0])
    old_idx = 0
    for i in range(0, len(fl) + 1, n_atoms + 2):
        if i == old_idx:
            continue
        atoms, geom, en = _parse_xyz_str(fl[old_idx:i])
        confs.append(geom)
        ens.append(en)
        old_idx = i

    return confs, atoms, ens

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


confs, atoms, ens = read_ensemble(args.file)

ref = confs.pop(args.reference)
rmsds = [0, ]

for i in confs: 
    rmsds.append(np.round(calc_rmsd(ref, i), 3))

ens = np.array(ens)
ens_rel = np.round((ens - np.min(ens))*627.51, 2)

data = {'RMSD': rmsds, "∆E [kcal/mol]":ens_rel}

df = pd.DataFrame(data)
print(df, end='\n\n')
print(f"Mediana: \nIdx: {len(df)//2}\n", df.median(), end='\n\n')
print("Media: ", df.mean(), end='\n\n')
df.to_csv(args.output.split('.')[0]+'.csv')