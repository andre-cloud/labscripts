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

parser = argparse.ArgumentParser(description='''
This script allows to calculate the RMSD on multiple structures. The RMSD is calculated with three different methods:\n
- Classic: two structures are aligned on avarege and the RMSD is calculated for each atoms on couple (ix, jx)
- CW: after the calculation of the center of weight, the RMSD is calculated considering the difference on the atom's distance to the center of weight
- VMD: this algorithm is the same of the software VMD. It's a "classic" algorithm bet weighted on the atom mass
''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('reference', help='The referement molecule [xyz format].')
parser.add_argument('files', help='The files to compare [xyz format].', nargs='+')
parser.add_argument('-i', '--indexes', help='Index for alignment of the molecules. Count starts a 0. Default %(default)s', nargs=3, type=int, default=[])

args = parser.parse_args()


def calc_rmsd(ref: np.array, conf: np.array, *args, **kwargs):
    return np.sqrt(1/ref.size) * np.linalg.norm(ref-conf)


def align_structures(structures:np.array, indexes=None):

    reference = structures[0]
    targets = structures[1]
    if isinstance(indexes, (list, tuple)):
        indexes = np.array(indexes)

    if not(indexes is None or len(indexes) == 0):
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

    dist_rif = np.sqrt(np.sum(np.linalg.norm(ref-mass_ref, axis=1 )**2)/len(ref))
    dist_file = np.sqrt(np.sum(np.linalg.norm(file-mass_file, axis=1 )**2)/len(ref))

    return calc_rmsd(dist_rif, dist_file)


def weight_rmsd(ref, file, numbers):
    masses = np.array([pt[i].mass for i in numbers]).T
    return np.sqrt(np.sum(masses*(np.linalg.norm((file - ref), axis=1))**2)/(len(ref)*np.sum(masses)))
    # mass = np.array([masses for _ in range(3)]).T



ref = read(args.reference)
ref_p = ref.positions[:]


print('Number of structurs evaluated: '+str(len(args.files)))
if args.indexes:
    print('Molecule for RMSD-Aligned calculations aligned on ' + str(tuple(args.indexes)))
else:
    print('Molecule align with no index')

print('='*30)
print('{:20s}\t{:10s}\t{:10s}\t{:10s}'.format('CONF', 'Classic', 'CW', 'VMD'))
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
    out = align_structures((ref.positions, file.positions), args.indexes)
    out = align_structures((ref.positions, out), args.indexes)
    file.positions = np.array(out)

    vmd_rmsd = weight_rmsd(ref.positions, file.positions, file.numbers)


    # print(re.split('\\\\|/', i))
    print('{:20s}\t{:.5f} \t{:.5f} \t{:.5f}'.format(re.split('\\\\|/', i)[-1], classic , rmsd_m, vmd_rmsd))
