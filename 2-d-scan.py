#! /usr/bin/python3

import argparse
import os
import sys

import numpy as np

from functions.io import parse_ensemble, mkdir

prm = {2: 'distance', 3: 'angle', 4:'dihedral'}




def required_length(nmin,nmax):
    class RequiredLength(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not nmin<=len(values)<=nmax:
                msg='argument "{f}" requires between {nmin} and {nmax} arguments'.format(
                    f=self.dest,nmin=nmin,nmax=nmax)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return RequiredLength


def file_validator():
    class RequiredGeometryExstention(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not ('xyz' in values and os.path.exists(values)):
                msg='argument "{f}" requires an exsiting xyz geometry file'.format(
                    f=self.dest)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return RequiredGeometryExstention


def scan(scan_prm, const, fix, first=False, add=None) -> None:

    *idx_1, start, end, step = scan_prm
    *idx_2, val = const
    print(idx_1, idx_2, start, end, step, val, fix)
    p1, p2 = prm[len(idx_1)], prm[len(idx_2)]
    fix_txt = f'\n$fix\natoms: {", ".join([str(i) for i in fix])}' if fix else ''
    add_txt = f'\n{add}' if add else ''
    txt = f'''$constrain
force constant=1.5
{p1}: {', '.join(list(np.array(idx_1, dtype=str)))}, {start}
{p2}: {', '.join(list(np.array(idx_2, dtype=str)))}, {val}'''+add_txt+fix_txt+f'''
$scan
1: {start}, {end}, {step}
end'''
    
    with open('input', 'w') as f:
        f.write(txt)

    return None


def parse():

    parser = argparse.ArgumentParser()

    parser.add_argument('file', help='Initial geometry. XYZ file required', action=file_validator())
    parser.add_argument('-idx_1', '--index_1', nargs='+', action=required_length(2, 4), help='Indexes for the scan. Starts the count at 1', required=True)
    parser.add_argument('-v_1', '--values_1', nargs='+', action=required_length(2, 2), help='Starting and Ending value for the scan of parameter 1', required=True)
    parser.add_argument('-idx_2', '--index_2', nargs='+', action=required_length(2, 4), help='Indexes for the scan. Starts the count at 1', required=True)
    parser.add_argument('-v_2', '--values_2', nargs='+', action=required_length(2, 2), help='Starting and Ending value for the scan of parameter 2', required=True)
    parser.add_argument('-s', '--step', nargs='+', action=required_length(2,2), help='Number of step for each scan', required=True)
    parser.add_argument('-f', '--fix', nargs='+', help='Index of atoms to fix')
    parser.add_argument('-add', '--add', help='Add strings in the constrain block. No sanity check is run')

    parser.add_argument('-l', '--launch', action='store_true')


    return parser.parse_args()

def get_topo():
    with open('gfnff_adjacency') as f:
        return f.read()

def check_topo(scan, idx, io, launch, TOPO):
    mkdir(f'topo_{idx}')
    os.chdir(f'topo_{idx}')
    with open('geom.xyz', 'w') as f:
        f.write(io)
    if launch: os.system('xtb topo geom.xyz > xtb.out')
    topo = get_topo()
    if not TOPO == topo:
        print(scan, idx, 'Not the same topo', sep=' - ')
        sys.exit()
    os.chdir('..')
    

        

def main():

    args = parse()

    
    idx_x = np.array(args.index_1, dtype=np.int32)
    idx_y = np.array(args.index_2, dtype=np.int32)
    v_x = np.array(args.values_1, dtype=float)
    v_y = np.array(args.values_2, dtype=float)
    grid = np.meshgrid(np.linspace(*v_x, int(args.step[0])), np.linspace(*v_y, int(args.step[1])))
    z = np.zeros(tuple(np.array(args.step, dtype=np.int32)))


    np.save('grid.npy', grid)

    PATH = os.getcwd()
    
    mkdir('topo')
    os.chdir('topo')
    if args.launch: os.system(f'cp ../{args.file} geom.xyz; xtb topo geom.xyz > xtb.out')
    TOPO = get_topo()

    os.chdir(PATH)
    
    mkdir('first_scan')
    os.chdir('first_scan')
    os.system(f'cp ../{args.file} geom.xyz')
    scan(list(idx_x)+list(v_x)+[args.step[0]], list(idx_y)+[v_y[0]], first=True, fix=args.fix, add=args.add)
    if args.launch: os.system(f'xtb geom.xyz --opt --alpb toluene --input input > xtb.out')
    for idx, j in enumerate(parse_ensemble('xtbscan.log')):
        check_topo('fisrt', idx, j, args.launch, TOPO)
    os.chdir(PATH)

    
    if args.launch:

        os.chdir(PATH)
        geoms = parse_ensemble('first_scan/xtbscan.log')
        for idx, i in enumerate(list(np.linspace(*v_x, int(args.step[0])))):

            d = 'scan_{:.3f}'.format(i)
            mkdir(f'{d}')
            os.chdir(d)

            with open('geom.xyz', 'w') as f:
                f.write(geoms[idx])
            scan(list(idx_y)+list(v_y)+[args.step[1]], list(idx_x)+[i], fix=args.fix, add=args.add)
            if args.launch: os.system(f'xtb geom.xyz --opt --alpb toluene --input input > xtb.out')
            for idx, j in enumerate(parse_ensemble('xtbscan.log')):
                check_topo(i, idx, j, args.launch, TOPO)

            en = np.array([float(i.splitlines()[1].split()[1]) for i in parse_ensemble('xtbscan.log')])
            z[:, idx] = en

            np.save(os.path.join(PATH,'en.npy'), z)
            os.chdir(PATH)

    return





if __name__=='__main__':
    
    main()
