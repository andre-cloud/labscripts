import argparse
from collections import Counter
import re, os
import sys
import shutil


import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()

parser.add_argument('file', help='Scan log file with 1D scan', nargs='+')
parser.add_argument('type', choices=['parser','writer', 'p', 'w'], help='Define which functionality of the program you want to use')
parser.add_argument('-i', '--index', help='FOR WRITER PART: define the index (starting at 1) of the atoms', nargs='+', required=True)


parser.add_argument('-gp', action='store_true', help='Obtain the two ground potential geometries')
parser.add_argument('-s', '--start', help='The initial point of the scan. Default auto: read from the input file the parameter', type=float, default=0)
parser.add_argument('-e', '--end', help='The final point of the scan', type=float)
parser.add_argument('-st', '--step', help='Insert the number of step you\'ve scanned', type=int, required=True)
parser.add_argument('--output', default=' ', help='Output filename without extention. Default: scan_INPUTFILE.dat')
parser.add_argument('-fc', '--force_constant', default=1.5, help='Define the force constant for the scan')
parser.add_argument('-in', '--increase', help='Define the increase from the auto parameter', type=float)
parser.add_argument('-ex', '--estention', help="Set the estantion of the dihedral. Default %(default)s", default=40, type=int)
parser.add_argument('-m', '--maxs', help="Set the number of maxies required. Default %(default)s", default=2, type=int)
parser.add_argument('-thr', '--threshold', help='Energy value below which no maximum can be found (in kcal/mol). Default: %(default)s', default=0.005, type=float)

args = parser.parse_args()
sys.setrecursionlimit(2000)

def find_max(x, y, th=10, add=False):
    if th!=10:
        th -= 0.4 if not add else -0.4
    m = min(y)
    maxs = []
    for idx, i in enumerate(x):
        back, fowr = idx-1, idx+1
        back_2, fowr_2 = idx-2, idx+2
        if back<0: back_2=-2; back=-1
        if back_2<0 and back>=0: back_2=-1
        if fowr>=len(x): fowr=0; fowr_2=1
        if fowr_2==len(x): fowr_2=0 
        if y[idx] > y[back] and y[idx] > y[fowr] and y[idx] > m+args.threshold and y[idx] > y[back_2] and y[idx] > y[fowr_2] and (np.abs(y[idx]-y[fowr])>th/627.51 or np.abs(y[idx]-y[back])>th/627.51): maxs.append(i)
    if len(maxs) != args.maxs and 0.4 < th < 25:
        flag = len(maxs) < args.maxs
        if th == 10: th -= 0.4 if flag else -0.4
        maxs = find_max(x, y, th, add=not flag)
    return maxs

def check_x(x, y):
    x_ = list(x)
    for idx, i in enumerate(x_):
        if i<-180 or i>180: 
            new_x = (360-abs(i)) * (1 if i<-180 else -1)
            x_[idx] = new_x
    dict_ = {k:v for k,v in zip(x_, y)}
    dict_ = {k: v for k, v in sorted(list(dict_.items()))}
    round_x = np.round(np.array(list(dict_.keys())), 3)
    if len(round_x) != len(set(round_x)):
        keys=[]
        dup = [item for item, count in Counter(round_x).items() if count > 1]
        for i in dup:
            for j in dict_:
                if i-j < 0.002:
                    keys.append(j)
            min_value=min([dict_[key] for key in keys])
            dict_[keys[0]] = min_value
            dict_.pop(keys[1])
            keys=[]

    return np.array(list(dict_.keys())), np.array(list(dict_.values()))



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


def get_position(file, indexes, read=True, geom=False):
    if read:
        with open(file) as f:
            fl = f.readlines()[2:]
    elif not geom:
        fl = file.split('\n')[1:]
    else:
        fl = file.split('\n')[2:]

    return np.array([fl[int(i)-1].split()[1:] for i in indexes], dtype=np.float64)


def write_gsm_input(ndoes=15):
    with open('ograd', 'w') as f:
        str_ = '''#!/bin/bash

if [ -z $2 ]
then
  echo " need two arguments! "
  exit
fi

ofile=orcain$1.in
ofileout=orcain$1.out
molfile=structure$1
ncpu=$2
basename="${ofile%.*}"

########## XTB settings: #################
cd scratch
wc -l < $molfile > $ofile.xyz
echo "Dummy for XTB/TM calculation" >> $ofile.xyz
cat $molfile >> $ofile.xyz

xtb $ofile.xyz --grad --chrg 0 > $ofile.xtbout

tm2orca.py $basename
rm xtbrestart
cd ..
'''     
        f.write(str_)
    
    with open('inpfileq', 'w') as f:
        str_ = f'''------------- QCHEM Scratch Info ------------------------
$QCSCRATCH/    # path for scratch dir. end with "/"
GSM_go1q       # name of run
---------------------------------------------------------

------------ String Info --------------------------------
SM_TYPE                 GSM    # SSM, FSM or GSM
RESTART                 0      # read restart.xyz
MAX_OPT_ITERS           160    # maximum iterations
STEP_OPT_ITERS          30     # for FSM/SSM
CONV_TOL                0.0005 # perp grad
ADD_NODE_TOL            0.1    # for GSM
SCALING                 1.0    # for opt steps
SSM_DQMAX               0.8    # add step size
GROWTH_DIRECTION        0      # normal/react/prod: 0/1/2
INT_THRESH              2.0    # intermediate detection
MIN_SPACING             5.0    # node spacing SSM
BOND_FRAGMENTS          1      # make IC's for fragments
INITIAL_OPT             0      # opt steps first node
FINAL_OPT               150    # opt steps last SSM node
PRODUCT_LIMIT           100.0  # kcal/mol
TS_FINAL_TYPE           0      # 0=no bond breaking, 1=breaking of bond
NNODES                  {ndoes}     # including endpoints
---------------------------------------------------------

'''
        f.write(str_)

    os.system('chmod +x ograd')


def mkdir(directory):
    if os.path.exists(directory):
        shutil.rmtree(directory)
    os.mkdir(directory)

def get_energy(filename):
    with open(filename, 'r') as f:
        fl = f.read()
    
    first = fl.split('\n')[0]
    regex = r'energy:'
    structurs = re.split(regex, fl, maxsplit=args.step)
    data = []
    ens = []
    geoms = []
    old_x = []
    for i in structurs[1:]:
        atoms = get_position(i, args.index, read=False)
        if len(args.index)==2:
            auto = distance(atoms)
        if len(args.index)==3:
            auto = angle(atoms)
        if len(args.index)==4:
            auto = dihedral(atoms)
        energy = i.split('\n')[0]
        n = str(len([j for j in i.split('\n')[1:] if j and len(j.split())==4]))
        geom =  f' {n}'
        # tmp = '\n'.join(i.split('\n')[1:-1])
        tmp = '\n'.join([j for j in i.split('\n')[1:] if j and len(j.split())==4])
        geom += '\n\n' + tmp
        ens.append(float(energy.split()[0]))
        geoms.append(geom)
        old_x.append(auto)
        data.append((auto, energy.split()[0]))
    
    if args.gp:
        x, y = check_x(np.array(data, dtype=np.float64)[:, 0], np.array(data, dtype=np.float64)[:, 1])
        st, en = x[0], x[-1]
        step = (en-st)/args.step
        th = int(args.estention/step)
        max_y = find_max(x, y)
        max_x = [list(old_x).index(i) for i in max_y] # index of x for maxes
        print(f'{filename} has {str(len(list(max_y)))} max(s) at: {", ".join([str(i) for i in max_y])}')
        return data
        DIR = os.getcwd()
        try:
            os.chdir(os.path.split(filename)[0])
        except Exception:
            pass

        if len(max_x) > 2: print(f'CAUTION: {filename} has more than two TSs. Remove the TSs in which you\'re not interested into.')
        mkdir('gsm')
        os.chdir('gsm')
        for index, i in enumerate(max_x): 
            start_geom_idx = i-th
            end_geom_idx = i+th if (i+th)<len(geoms) else th-(len(geoms)-1-idx)
            mkdir(f'ts_{index}')
            os.chdir(f'ts_{index}')
            mkdir(f'scratch')
            try:
                with open(os.path.join('scratch','initial0000.xyz'), 'w') as f:
                    str_ = geoms[start_geom_idx].strip() + '\n' + geoms[end_geom_idx].strip()
                    f.write(str_)
            except IndexError:
                print((start_geom_idx, end_geom_idx, len(geoms)-1))
            write_gsm_input()
            os.chdir('..')
            
        
        os.chdir('..')
        os.chdir(DIR)

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


def write_input(i):
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



text ='''$constrain 
 force constant={fc} 
 {prm}: {index}, {auto} 
$scan 
 1: {auto}, {end}, {step} 
$end '''


if __name__=='__main__':
    if args.type in ['parser', 'p']:
        for idx, i in enumerate(args.file):
            en = get_energy(i)
            write_file(en, idx)
    
    if args.type in ['writer', 'w']:
        if args.index:
            for idx, i in enumerate(args.file):
                write_input(i)
                
