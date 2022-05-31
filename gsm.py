import argparse
from collections import Counter
import os
import sys
import shutil

import matplotlib.pyplot as plt
import numpy as np


parser = argparse.ArgumentParser()

parser.add_argument('file', nargs='+')
parser.add_argument('-i', '--index', help='Define the index (starting at 1) of the atoms', nargs='+', required=True)


args = parser.parse_args()


def read_input(filename):
    with open(filename) as f:
        fl = f.read()
    fl = fl.splitlines()
    points = []
    prev_i = 0
    for i in range(0, len(fl), int(fl[0])+2):
        if fl[prev_i:i]: points.append('\n'.join(fl[prev_i:i])) 
        prev_i=i
    data = []
    ens = []
    for point in points:
        at = get_position(point, args.index)
        dh = dihedral(at)
        en = float(point.splitlines()[1])
        ens.append(en)
        data.append((dh, en))

    geom_max = points[ens.index(max(ens))]
    return np.array(data, dtype=np.float64), geom_max, points


def get_position(file, indexes):
    fl = file.split('\n')[2:]
    return np.array([fl[int(i)-1].split()[1:] for i in indexes], dtype=np.float64)


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

def check_for_maxs(y, geom):
    idx = list(y).index(max(y))
    flag = idx != 0 and idx != len(y)
    at = get_position(geom[idx], args.index)
    dih = np.abs(dihedral(at))
    th = 35
    return (flag and (dih<=th or dih>=180-th))

def write_output(geom, del_=False):
    if os.path.exists("ensemble.xyz") and del_:
        os.remove("ensemble.xyz")
    with open('ensemble.xyz', 'a+') as f:
        f.write(geom+'\n')

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





xs, ys, labels = [], [] , []
if __name__ == '__main__':
    plots = []
    labels = []
    pops = []
    for idx, i in enumerate(args.file):
        data, geom, points = read_input(i)
        x, y = data[:, 0], data[:, 1]
        if check_for_maxs(y, points):
            x, y = check_x(x, y)
            xs.append(x)
            ys.append(y)
            labels.append(i)
            write_output(geom, del_=True if idx == 0 else False)
        else: pops.append(idx)

    for i in pops[::-1]: args.file.pop(i)
    print(len(xs))

    for idx, i in enumerate(zip(xs, ys)):
        x, y = i
        l = os.path.join(args.file[idx].split('/')[0],args.file[idx].split('/')[2])
        plot = plt.plot(x, y, '-o', label=args.file[idx])
        plots.append(plot)
        labels.append(args.file[idx])
    
    fig = plt.gcf()
    legend = plt.legend(
        loc='upper center', bbox_to_anchor=(0.5, -.125), fancybox=True, shadow=True, ncol=3
    )

    for obj in legend.__dict__['legendHandles']:
        obj.set_picker(5)

    def on_pick(event):
        flag = False
        groups = {key: val for key, val in zip(labels, plots)}
        label = event.artist.get_label()
        # print(label)
        # print(label in groups)
        for key in groups:            
            if key == label:
                if not groups[key][0].__dict__['_alpha'] or groups[key][0].__dict__['_alpha'] != 1:
                    groups[key][0].set_alpha(1.0)
                else:
                    flag = True
            else:
                groups[key][0].set_alpha(0.1)
        if flag:
            for key in groups:
                groups[key][0].set_alpha(None)
            flag = False

        fig.canvas.draw()

    plt.connect('pick_event', on_pick)


    plt.tight_layout()
    plt.show()