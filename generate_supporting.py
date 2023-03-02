import argparse, cclib, getpass, os


DEBUG = getpass.getuser() == 'andrea'

def parse():
    parser = argparse.ArgumentParser()

    parser.add_argument('file', nargs='+')
    parser.add_argument('-n', '--name', nargs='+', required=True)
    
    parser.add_argument('--output', help='Define the output file. Default "%(default)s"', default=f'/data/{getpass.getuser()}/supporting.txt' if not DEBUG else 'supporting.txt')
    parser.add_argument('--separater', help='Define the separator between geometries. Defualt 2 new line', default='\n\n')

    return parser.parse_args()


def parse_output(file, name, sep):
    data = cclib.io.ccread(file)
    geom = "\n".join(data.writexyz().splitlines()[2:])
    txt = f'''{name:<}
Inner Energy (Hartree)      : {data.scfenergies[-1]:<5.10f}
Zero Point Energy (Hartree) : {data.zpve:<5.10f}
Enthalpy (Hartree)          : {data.enthalpy:<5.10f}
Entropy (Hartree)           : {data.entropy:<5.10f}
Gibbs Energy (Hartree)      : {data.freeenergy:<5.10f}
Immaginary Frequency (cm-1) : {','.join([f'{i:.2f}' for i in data.vibfreqs if i<0]) if data.vibfreqs[data.vibfreqs<0] else "None"}

{geom}
{sep}
    '''
    return txt


def write_in_file(geom, file):
    with open(file, 'a') as f:
        f.write(geom)
    return


if __name__=='__main__':
    args = parse()
    assert len(args.name) == len(args.file), f"Lenght of list FILES (input lenght {len(args.file)}) and NAMES (given {len(args.name)}) must match."

    if not os.path.exists(args.output): 
        with open(args.output, 'w') as f:
            pass 
        print(f'File {args.output} have been generated.')
    else:
        print(f'File {args.output} have been found, so going to append geometries on that file.')

    for file, name in zip(args.file, args.name):
        geom = parse_output(file, name, args.separater)
        write_in_file(geom, args.output)
    
    print('All geometries have been parsed.')
