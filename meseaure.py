import argparse
from functions import *



func = {
    2: geom_prm.distance,
    3: geom_prm.angle,
    4: geom_prm.dihedral,
}
parser = argparse.ArgumentParser()

parser.add_argument('file', nargs='+')
parser.add_argument('-idx', '--indexes', help='Define indexes for the measure. One measure at a time', nargs='+', action=argparse_class.required_length(2, 4))

args = parser.parse_args()


if __name__ == '__main__':
    ruler = func[len(args.indexes)]
    for file in args.file:
        confs = io.parse_ensemble(file)
        for idx, i in enumerate(confs):
            p = get_position(i, args.indexes, read=False, geom=True)
            print(f'{file} - {idx} - {ruler(p):.2f}')




