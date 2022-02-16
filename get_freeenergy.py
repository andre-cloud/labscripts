try:
    import cclib, os, argparse
    from alive_progress import alive_bar
    import numpy as np
except ImportError:
    os.system('pip install -r requirements.txt')
    import cclib, os, argparse
    from alive_progress import alive_bar
    import numpy as np


parser = argparse.ArgumentParser(
    description='''This script is meant to read defined log files and get ∆G from them or do a path walk gaining from all log files encountered the ∆G'''
)

parser.add_argument('-f', '--files', nargs='+', default='.')
parser.add_argument('-o', '--output', default='/Users/andreapellegrini/Public/gibbs_energy.txt')
parser.add_argument('-T', '--test', action='store_true')
parser.add_argument('-t', '--title', nargs='+', default=['Entry', 'Filename', 'Filepath', 'G (Hertree)', '∆G (kcal/mol)'])

args = parser.parse_args()

H = 627.5095

def walk():
    files = {}
    for directory, dirnames, filename in os.walk(os.getcwd()):
        with alive_bar(len(filename)) as bar:
            for file_ in filename:
                if (file_.endswith('.log') and ('scan' not in file_.lower() or 'ts' in file_.lower())):
                    fp = os.path.join(directory, file_)
                    parser = cclib.io.ccopen(fp)
                    data = parser.parse()
                    try:
                        files[fp] = data.freeenergy
                        print(f'{file_.ljust(30)} ∆G={data.freeenergy} Hartree')
                    except AttributeError:
                        pass
                bar()

    return files


def define_files():
    files = {}
    with alive_bar(len(args.files)) as bar:
        for file_ in args.files:
            fp = os.path.join(os.path.join(os.getcwd(), file_))
            parser = cclib.io.ccopen(fp)
            data = parser.parse()
            try:
                files[fp] = data.freeenergy
                print(f'{file_.ljust(30)} ∆G={data.freeenergy} Hartree')
            except AttributeError:
                pass
            bar()

    return files

def get_energys():
    if args.test:
        files = np.array([['/Volumes/shares/gau/AP/nn_boc2_aly2/clockwise/trans_ts_b3lyp_gas.log',
        '-1036.622672'],
       ['/Volumes/shares/gau/AP/nn_boc2_aly2/clockwise/cis_ts_b3lyp_gas.log',
        '-1036.617465'],
       ['/Volumes/shares/gau/AP/nn_boc2_aly2/opt/opt_b3lyp.log',
        '-1036.66465']], dtype='<U68')
    elif args.files == '.':
        files = walk()
    else:
        files = define_files()

    files = np.array(list(files.items())) if type(files) is dict else files

    e_ = files[:, -1].astype(np.float64)
    de = (e_-min(e_)) * H
    filepaths = list(files[:,0])
    for i in range(len(filepaths)):
        filepaths[i] = os.path.join(os.getcwd(), filepaths[i])
    tmp = np.vstack((e_, de)).T
    for idx, i in enumerate(list(files[:, 0])):
        _, tail = os.path.split(i)
        files[idx, 0] = tail
        files[idx, 1] = filepaths[idx]
    files = np.concatenate((files, tmp), axis=1)
    return files


def write_file(files):
    with open(args.output, 'w') as f:
        f.write(','.join(args.title))
        f.write('\n')
        for idx, i in enumerate(files):
            f.write(','.join([str(idx)] + list(i)))
            f.write('\n')

def main():
    files = get_energys()
    write_file(files)


if __name__ == '__main__':
    main()
