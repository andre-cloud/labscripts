import matplotlib.pyplot as plt 
import numpy as np
import argparse, sys

parser = argparse.ArgumentParser()
parser.add_argument('filename', help='The interpolation file generated by ORCA', default='neb.interp', nargs='+')
parser.add_argument('-s', '--silent', help='Flag if you want to only see the graph generated without saving it', action='store_true')
parser.add_argument('-o', '--output', help='Output filename for the graph')


args = parser.parse_args()

H = 627.5096080305927


def parse_file(filename):
    with open(filename) as f:
        file = f.read()

    iters = file.split('Iteration')[1:]
    imags, interpols = [] , []
    for i in iters:
        if not i: continue
        imag, interp, *_ = i.split('''

''')
        img = np.array( [i.split() for i in imag.split('\n')[2:]])
        x, _, y = img[:, 0], img[:, 1], img[:, 2]
        imags.append(np.array([(float(i), float(j)*H) for i, j in zip(x,y)]))
        intr = np.array( [i.split() for i in interp.split('\n')[2:] if i])
        x, _, y = intr[:, 0], intr[:, 1], intr[:, 2]
        interpols.append(np.array([(float(i), float(j)*H) for i, j in zip(x,y)]))
    
    return imags, interpols



for file in args.filename:

    imags, interpols = parse_file(file)

    for idx, i in enumerate(interpols):
        x, y = i[:, 0], i[:, 1]
        c = 'gray' if (idx != 0 and idx != len(interpols)-1) else 'black' if idx == 0 else 'red'
        alpha = 0.1 if (idx != 0 and idx != len(interpols)-1) else 0.8 if idx == 0 else 1
        w = 0.2 if (idx != 0 and idx != len(interpols)-1) else 0.8 if idx == 0 else 1
        plt.plot(x, y, color=c, alpha=alpha, linewidth=w)
        
    for idx, i in enumerate(imags):
        x, y = i[:, 0], i[:, 1]
        c = 'gray' if (idx != 0 and idx != len(imags)-1) else 'black' if idx == 0 else 'red'
        alpha = 0.15 if (idx != 0 and idx != len(imags)-1) else 0.8 if idx == 0 else 1
        s = 20 if (idx != 0 and idx != len(interpols)-1) else 30 if idx == 0 else 35
        plt.scatter(x, y, color=c, alpha=alpha, s=s)


    plt.title(f'NEB Analysis - {len(imags)} iterations')
    plt.xlabel('Reaction coordinate')
    plt.ylabel('∆E [kcal/mol]')

    # plt.xlim((0, 1))

    plt.tight_layout()
    if args.silent:
        plt.show()
    else:
        plt.savefig(f'{file.split(".")[0]}.png', dpi=300)


