import numpy as np


def get_geometry(points):
    geom = np.array([i.split()[1:] for i in points], dtype=np.float64)
    return geom*[-1, 1, 1]


def read_xyz(file): 
    with open(file) as f:
        fl = f.read()
    fl = fl.splitlines()
    points = []
    prev_i = 0
    for i in range(int(fl[0].strip())+2, len(fl)+1, int(fl[0].strip())+2):
        if prev_i != 0:
            if fl[prev_i:i]: points.append('\n'.join(fl[prev_i:i])) 
        else:
            points.append('\n'.join(fl[:i])) 
        prev_i=i
    return points

def write_xyz(inp, out, atoms, en_geom):
    with open(out, 'a') as f:
        f.write(str(len(en_geom))+'\n')
        f.write('Enantiomer of {}\n'.format(inp))
        for a, i in zip(atoms, en_geom):
            f.write('%s \t %.5f \t %.5f \t %.5f \n' % (a, *i))
        
            
if __name__=='__main__':
    file = 'i.xyz'
    gs = read_xyz(file)
    for g in gs:
        g = g.splitlines()
        atoms = [i.split()[0] for i in g[2:] if i]
        en = get_geometry(g[2:])
        write_xyz(file, 'en_i.xyz', atoms, en)
    
