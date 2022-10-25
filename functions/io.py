import os, re
import shutil

def parse_ensemble(file:str) -> list:
    """
    Parse an ensemble xyzfile

    param: file: filename to be parse

    return: list of poses
    """

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



def mkdir(directory) -> None:
    """
    Create a directory. If it exists, it will overwrite the directory

    param: directory: name of the directory

    return: None
    """
    if os.path.exists(directory):
        shutil.rmtree(directory)
    os.mkdir(directory)
    return None