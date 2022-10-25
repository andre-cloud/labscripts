import os, sys


import numpy as np
from functions.io import parse_ensemble
np.set_printoptions(threshold=sys.maxsize)

THR = 1

def distance(v, w):
    return np.linalg.norm(w - v)


def create_graph(n_atoms, zeros=False):
    if zeros: return np.zeros((n_atoms, n_atoms))    
    return np.empty((n_atoms, n_atoms))


def check_if_calculated(mat, i, j):
    return mat[j,i]


def get_geometry_matrix(points):
    return np.array([i.split()[1:] for i in points.splitlines()[2:]], dtype=np.float64)


def fill_dist_matx(n_atoms, geom):
    dist = create_graph(n_atoms)
    for i in range(n_atoms):
        for j in range(n_atoms):
            if i == j: 
                dist[i, j] = 0
            elif not check_if_calculated(dist, i, j):
                dist[i, j] = distance(geom[i], geom[j])
            else:
                dist[i, j] = dist[j, i]
    return dist


def fill_adj_matrix(n_atoms, dist):
    adj = create_graph(n_atoms, zeros=True)
    idx = 0
    for x in dist:
        ma = list(np.nonzero(x<THR)[0])
        for m in ma:
            if idx == m:
                continue
            adj[idx, m] = m
        idx += 1
    return adj




def create_topo(fname):
    for idx, fl in enumerate(parse_ensemble(fname)):
        geom = get_geometry_matrix(fl)
        n_atoms = len(geom)
        dist = fill_dist_matx(n_atoms, geom)
        adj = fill_adj_matrix(n_atoms, dist)
        np.save(f'topo_{os.path.split(fname)[-1].split(".")[0]}_{idx}.npy', adj)


        
def check_topo(fnames):
    topos = []
    for fname in fnames:
        topos.append(np.load(fname))
    
    for i in topos: 
        for j in range(topos.index(i), len(topos)):
            print(topos.index(i), j)
            if not np.allclose(i, topos[j]):
                return False
        
    return True







if __name__ == '__main__':
    create_topo('/Users/andrea/Desktop/scratch/scan/xtbscan_bb_gp.xyz')
    create_topo('/Users/andrea/Desktop/scratch/scan/gp.xyz')
    fs = [
        f'topo_xtbscan_bb_gp_{i}.npy' for i in range(0, 36)
    ]+['topo_gp_0.npy']
    print(check_topo(fs))