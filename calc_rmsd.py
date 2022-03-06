from ase.io import read, write
from ase.build import minimize_rotation_and_translation
import numpy as np
import os, argparse, re
from rmsd import *
from rmsd import kabsch
from periodictable import core, covalent_radius, mass

pt = core.PeriodicTable(table="H=1")
covalent_radius.init(pt)
mass.init(pt)

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--reference', help='The referement molecule [xyz format].', required=True)
parser.add_argument('-f', '--files', help='The files to compare.', required=True, nargs='+')

args = parser.parse_args()


def calc_rmsd(ref: np.array, conf: np.array, *args, **kwargs):
    return np.sqrt(1/ref.size * np.linalg.norm(ref-conf))


def align_structures(structures:np.array, indexes=None):

    reference = structures[0]
    targets = structures[1]
    if isinstance(indexes, (list, tuple)):
        indexes = np.array(indexes)

    indexes = slice(0,len(reference)) if indexes is None or len(indexes) == 0 else indexes.ravel()

    reference -= np.mean(reference[indexes], axis=1)
    targets -= np.mean(targets[indexes], axis=1)

    matrix = kabsch(reference, targets)
    return (matrix @ targets.T).T


def rmsd_with_center_mass(ref, file, numbers):

    masses = [pt[i].mass for i in numbers]
    mass = np.array([masses for _ in range(3)]).T
    mass_ref = np.sum(ref*mass, axis=0)/np.sum(mass)
    mass_file = np.sum(file*mass, axis=0)/np.sum(mass)

    dist_rif = np.sqrt(np.sum(np.linalg.norm(ref-mass_ref, axis=1 )**4)/len(ref))
    dist_file = np.sqrt(np.sum(np.linalg.norm(file-mass_file, axis=1 )**4)/len(ref))

    return calc_rmsd(dist_rif, dist_file)


def weight_rmsd(ref, file, numbers):
    masses = np.array([pt[i].mass for i in numbers]).T
    return np.sqrt(np.sum(masses*(np.linalg.norm((file - ref), axis=1))**4)/(len(ref)*np.sum(masses)))
    # mass = np.array([masses for _ in range(3)]).T



ref = read(args.reference)
ref_p = ref.positions[:]

print('{:20s}\t{:10s}\t{:10s}\t{:10s}\t{:10s}'.format('CONF', 'Classic', 'CW', 'Aligned', 'VMD'))
for i in args.files:

    file = read(i)
    file_p = file.positions[:]


    minimize_rotation_and_translation(ref, file)
    classic = calc_rmsd(ref_p, file_p)


    nwp = []
    matrix = kabsch(ref_p, file_p)
    nwp.append([matrix @ vector for vector in file_p])
    file.positions = np.array(nwp)




    rmsd_m = rmsd_with_center_mass(ref_p, file_p, file.numbers)
    out = align_structures((ref.positions, file.positions), (62,63, 129))
    out = align_structures((ref.positions, out), (62,63, 129))
    file.positions = np.array(out)
    align = calc_rmsd(ref_p, file.positions)

    vmd_rmsd = weight_rmsd(ref.positions, file.positions, file.numbers)


    # print(re.split('\\\\|/', i))
    print('{:20s}\t{:.5f} \t{:.5f} \t{:.5f} \t{:.5f}'.format(re.split('\\\\|/', i)[-1], classic , rmsd_m, align, vmd_rmsd))
