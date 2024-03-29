import numpy as np
from scipy.constants import R


def boltzmann_perc(en, T):
    return np.exp(-(en-np.min(en))/R*T)/np.sum(np.exp(-(en-np.min(en))/R*T))


def get_positions(xyz):
    at = [atom for atom in xyz.splitlines()[2:]]
    return np.array([a.strip().split()[1:] for a in at], dtype=np.float64)



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