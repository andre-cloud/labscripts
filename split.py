import argparse, os, shutil
from alive_progress import alive_bar
import sys


parser = argparse.ArgumentParser()

parser.add_argument('file', help='Add an xyz file to parse', nargs='+')
parser.add_argument('-d', '--directory', help='Change the name of the directory. If more file selected it creates more directory, one for each file', default='structures')


args = parser.parse_args(['/Users/andreapellegrini/Desktop/scratch/crest_conformers.xyz'])

def split(filename):
    flag = False 

    if len(args.file) == 1:
        directory = args.directory
    else:
        directory = args.directory + str(args.file.index(filename))
        
    if os.path.exists(directory):
        if 'y' in input(f'A directory named {directory} already exists. Existing directory will be deleted, wanna procede? [y/n]').lower():
            shutil.rmtree(directory)   
        else:
            sys.exit()
    os.mkdir(directory)

    with open(filename) as f:
        file = f.read()
    
    parser = file.split('\n')[0]
    files = file.split(parser)
    
    if 'CONF' in file.split('\n')[1]:
        flag = True

    with alive_bar(len(files)-1, title=filename) as bar:
        for idx, file in enumerate(files, start=0):
            if not file or file == '':
                continue
            fn = file.split('\n')[1].split()[-1].strip('!')+'.xyz' if flag else 'CONF'+str(idx)+'.xyz'
            with open(os.path.join(directory, fn), 'w') as f:
                f.write(parser)
                f.write(file)
            bar()


    
    

if __name__=='__main__':
    for i in args.file:
        split(i)
    

