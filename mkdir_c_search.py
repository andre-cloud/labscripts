import argparse
import os, re
import shutil

dirs = ['xtb','crest','censo']
# print(os.getlogin())

PATH_CENSORC = '/Users/andrea/Desktop/labscripts/tests/.censorc' if os.getenv('DEBUG') else f'/data/{os.getlogin()}/.censorc'
with open(PATH_CENSORC) as f:
    _censo = f.readlines()

parser = argparse.ArgumentParser()

parser.add_argument('-solvs', help='Define the solvents for the calculations', nargs='+')
parser.add_argument('-diastero', help='Define the diasteroisomers on wich you want to conduct the search', nargs='+')

parser.add_argument('-func', help='Define the final DFT level, as ORCA required', required=True, nargs='+')
parser.add_argument('-basis', help='Define the final basis set, as ORCA required', required=True, nargs='+')
parser.add_argument('-charge', help='Define charge', type=int, default=0)

parser.add_argument('--no_solv', action='store_true')
parser.add_argument('--no_diastero', action='store_true')

args = parser.parse_args()

if len(args.basis) == 1:
    args.basis *= len(args.func)
else:
    assert len(args.basis) == len(args.func), f"Number of basis exspected the same of functinals ({len(args.func)}), got {len(args.basis)}"


sub_dir = {
    'charge:' : args.charge,
    'solvent:' : None,
    'func:' : None,
    'basis:' : args.basis,
}

subs = list(sub_dir.keys())

def replace_string(string):
    txt = string if string not in ['-','/','_','(',')'] else ''
    return txt

def filter_function(line):
    for i in subs:
        if i in line:
            return line


tmp = [''.join(list(map(replace_string, i))) for i in args.func]

dirs += tmp
# print(dirs)

def mkdir(directory):
    if os.path.exists(directory):
        shutil.rmtree(directory)
    os.mkdir(directory)
    return None


def modify_censorc(dir_, solvent='gas', func='r2scan-3c', basis='automatic', not_chg_func=True):
    _cens = _censo[:]
    sub_dir_tmp = sub_dir.copy()
    sub_dir_tmp['solvent:'] = solvent
    sub_dir_tmp['basis:'] = basis
    sub_dir_tmp['func:'] = func
    if not not_chg_func: sub_dir_tmp['part0:'] = 'off'

    for idx, i in enumerate(_cens):
        if '$' in i or i.strip() == '':
            continue
        if i.split()[0] not in sub_dir_tmp:
            continue
        _cens[idx] = re.sub(f'{" ".join(i.split()[:2])}', f'{i.split()[0].split()[0]} {sub_dir_tmp[i.split()[0].split()[0]]}', i)

    with open(os.path.join(dir_, '.censorc'), 'w') as f:
        f.writelines(_cens)

    return None

def make_dirs(solvent='gas'):
    for idx, x in enumerate(dirs):
        mkdir(x)
        if idx >= 2:
            modify_censorc(x, 
                solvent,  
                func='r2scan-3c' if idx == 2 else args.func[idx-3], 
                basis='automatic' if idx == 2 else args.basis[idx-3],
                not_chg_func= idx==2,
            )


if args.no_solv:
    if args.no_diastero:
        make_dirs()
    else:
        for i in args.diastero:
            mkdir(i)
            PATH = os.getcwd()
            os.chdir(i)
            make_dirs()
            os.chdir(PATH)
else:
    for j in args.solvs:
        PATH = os.getcwd()
        mkdir(j)
        os.chdir(j)
        if args.no_diastero:
            make_dirs(j)
        else:
            for i in args.diastero:
                mkdir(i)
                PATH_ = os.getcwd()
                os.chdir(i)
                make_dirs(j)
                os.chdir(PATH_)
        os.chdir(PATH)

