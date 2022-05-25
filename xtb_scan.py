import argparse
import re
import sys

import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument('file', help='Scan log file with 1D scan', nargs='+')
parser.add_argument('type', choices=['parser','writer'], help='Define which functionality of the program you want to use')
parser.add_argument('-i', '--index', help='FOR WRITER PART: define the index (starting at 1) of the atoms', nargs='+', required=True)

parser.add_argument('-s', '--start', help='The initial point of the scan', type=float, default=0)
parser.add_argument('-e', '--end', help='The final point of the scan', type=float)
parser.add_argument('-st', '--step', help='Insert the number of step you\'ve scanned', type=int, required=True)
parser.add_argument('--output', default=' ', help='Output filename without extention. Default: scan_INPUTFILE.dat')
parser.add_argument('-fc', '--force_constant', default=1.5, help='Define the force constant for the scan')
parser.add_argument('-in', '--increase', help='Define the increase from the auto parameter', type=float)

args = parser.parse_args()

def dihedral(p):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""
    p0 = p[0]
    p1 = p[1]
    p2 = p[2]
    p3 = p[3]

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


def angle(p):
    a = p[0]
    b = p[1]
    c = p[2]

    ab = np.linalg.norm(a-b)
    ac = np.linalg.norm(a-c)
    bc = np.linalg.norm(b-c)

    return np.cosh((ac**2+bc**2-ab**2)/(2*ab*bc))

def distance(p):
    a = p[0]
    b = p[1]
    
    return np.linalg.norm(a-b)


def get_position(file, indexes, read=True):
    if read:
        with open(file) as f:
            fl = f.readlines()[2:]
    else:
        fl = file.split('\n')[1:]
    
    return np.array([fl[int(i)-1].split()[1:] for i in indexes], dtype=np.float64)


def get_energy(filename):
    with open(filename, 'r') as f:
        fl = f.read()
    
    first = fl.split('\n')[0]
    regex = r'energy:'
    structurs = re.split(regex, fl, maxsplit=args.step)
    data = []
    for i in structurs[1:]:
        atoms = get_position(i, args.index, read=False)
        if len(args.index)==2:
            auto = distance(atoms)
        if len(args.index)==3:
            auto = angle(atoms)
        if len(args.index)==4:
            auto = dihedral(atoms)
        energy = i.split('\n')[0]
        data.append((auto, energy.split()[0]))
    return data


def write_file(energies, idx):
    if args.output != ' ':
        out = args.output+str(idx)+'.dat' if len(args.file) > 1 else args.output+'.dat'
    else:
        out = f'scan_{args.file[idx].split("/")[0]}.dat'
    with open(out, 'w') as f:
        for x, y in energies:
            f.write(f'{x}   {y}\n')
    return None

text ='''$constrain 
 force constant={fc} 
 {prm}: {index} {auto} 
$scan 
 1: {auto}, {end}, {step} 
$end '''


if __name__=='__main__':
    if args.type == 'parser':
        for idx, i in enumerate(args.file):
            en = get_energy(i)
            write_file(en, idx)
    
    if args.type == 'writer':
        if args.index:
            for idx, i in enumerate(args.file):
                atoms = get_position(i, args.index)
                if len(args.index)==2:
                    auto = distance(atoms)
                    prm = 'distance'
                if len(args.index)==3:
                    auto = angle(atoms)
                    prm = 'angle'
                if len(args.index)==4:
                    auto = dihedral(atoms)
                    prm = 'dihedral'
                    
                if args.start != 0 : auto = args.start

                if (args.end or args.increase) and args.step:
                    with open('scan.inp', 'w') as f:
                        f.write(text.format(
                            fc = args.force_constant,
                            prm = prm,
                            index = ', '.join(args.index),
                            auto = auto, 
                            end = str(args.end) if args.end else str(float(auto)+args.increase),
                            step = args.step
                        ))
                
