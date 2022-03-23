import argparse
from copyreg import pickle
import os
import re
import shutil
from statistics import mean
import subprocess
import sys


from alive_progress import alive_bar
import cclib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle


#TODO: Comparison between functionals


parser = argparse.ArgumentParser()

parser.add_argument('file', help='Log file(s) of the TD-SCF calculation', nargs='+')
parser.add_argument('-d', '--directory', help='Name for the new directory to store the files', default='logs')
parser.add_argument('-solv', '--solvation', help='Solvent used in the calculation', default='1')
parser.add_argument('-pv', '--python_version', help='Set the python command. Default: "%(default)s"', default='python')
parser.add_argument('-il', '--initial_lambda', help='Inital wavelenght. Default: %(default)s nm', default=200, type=int)
parser.add_argument('-fl', '--final_lambda', help='Final wavelenght. Default: %(default)s nm', default=600, type=int)
parser.add_argument('-def', '--definition', help='Definition of the spectra. Add more points on which calculate the absortion. (MAX-MIN)*100^d. Default: %(default)s', default=2, type=float)
parser.add_argument('-si', '--sigma', help='Peak width. Defaul %(default)s', default=10, type=float)
parser.add_argument('-p', '--pop', help='Define the population of the conformers indicated. Use with caution and be sure of the order. MAX population = 1', nargs='+', type=float)
parser.add_argument('-r', '--reference', help='File xy of the ECD plot sperimental')
parser.add_argument('-t', '--title', help='Set the title of the final plot. Default %(default)s', default='ECD graph')
parser.add_argument('--no_goodvibes', help='Don\'t execute goodvibes. If not -p, equal population will be considered', action='store_true')
parser.add_argument('--skip_data', help='Skip data to find the pick wavelenght', default=15, type=int)
parser.add_argument('-n', '--normalisation', help='Set the normalisation range. Default: [-%(default)s, %(default)s', default=1)
parser.add_argument('-sh', '--shift', help='Manually set the wavelenght shift in order to match the resulting spectra with the reference', type=float)
parser.add_argument('--save', help='Save pickle of the graph.', action='store_true')


args = parser.parse_args()

FILETOANALYSE = []
X = np.linspace(args.initial_lambda-100, args.final_lambda+100, (args.final_lambda-args.initial_lambda)*100**args.definition)

columns = ['fln', 'pop', 'R', 'l', 'conv']
DF = pd.DataFrame(columns=columns)


class InputError(Exception):
    def __init__(self, message):
        super().__init__(message)

if args.pop:
    if not 0.99 < sum(args.pop) < 1.01:
        print(f'Sum population = {sum(args.pop)*100}%. {"Not complete population give" if sum(args.pop) < 0.99 else "Too much population given"}. Prociding with unputted populations')


def get_rvel(filename:str):
    
    with open(filename) as f:
        fl = f.read()
    if re.findall(' <0|del|b> * <b|rxdel|0> + <0|del|b> * <b|delr+rdel|0>', fl) == ['del']:
        raise InputError(f'File {filename} it\'s not a TD-DFT calculation. Please check the input and re-run the code.')
    fl = fl.split(' <0|del|b> * <b|rxdel|0> + <0|del|b> * <b|delr+rdel|0>')[-1]
    fl = fl.split(' 1/2[<0|r|b>*<b|rxdel|0> + (<0|rxdel|b>*<b|r|0>)*]')[0]
    fl = fl.split('\n')[3:]
    return np.array([float(i.split()[4]) for i in fl if i])

def get_wavelenght(filename):
    with open(filename) as f:
        fl = f.read()
    if re.findall('Excitation energies and oscillator strengths', fl) == ['del']:
        raise InputError(f'File {filename} it\'s not a TD-DFT calculation. Please check the input and re-run the code.')
    fl = fl.split('Excitation energies and oscillator strengths')[-1]
    fl = fl.split('Population analysis using the SCF density.')[0]
    l = []
    for i in fl.split('\n'):
        if 'Excited State' in i:
            l.append(float(i.split()[6].strip()))
    return l

   

def split_file(filename: str):
    head_directory = os.path.split(filename)[0]
    directory = os.path.join(os.getcwd(), head_directory, args.directory)
    if len(args.file)>1:
        directory += '_'+os.path.split(filename)[1].split('.')[0]

    
    with open(filename) as f:
        file = f.read()
    if len(re.findall('Copyright', file))>1:
        files = file.split('Copyright')
        if os.path.exists(directory):
            if 'y' in input(f'A directory named {directory} already exists. Existing directory will be deleted, wanna procede? [y/n]').lower():
                shutil.rmtree(directory)
            else:
                sys.exit()
        os.mkdir(directory)
        for idx, i in enumerate(files):
            fln = os.path.join(directory, os.path.split(filename)[1].split('.')[0]+'_'+str(idx)+'.log')
            with open(fln, 'w') as f:
                f.write(i)
            if idx != 0: FILETOANALYSE.append(fln)
    else:
        FILETOANALYSE.append(filename)
    

def get_absolute_path(filename):
    for i in FILETOANALYSE:
        if os.path.split(i)[1].strip('.log') == filename:
            return i 

def get_filename(abs_path):
    return os.path.split(abs_path)[1]


def equalpop():
    with alive_bar(len(FILETOANALYSE), title='Getting λ R(vel) and population for files') as bar:
        for i in FILETOANALYSE:
            DF.loc[len(DF)] = [get_filename(i), 1/len(FILETOANALYSE), get_rvel(i), get_wavelenght(i), None]
            bar()

def get_population():

    if args.pop:
        with alive_bar(len(FILETOANALYSE), title='Getting λ R(vel) and population for files') as bar:
            for idx, i in enumerate(FILETOANALYSE):
                DF.loc[len(DF)] = [get_filename(i), args.pop[idx], get_rvel(i), get_wavelenght(i), None]
                bar()

    elif args.no_goodvibes or len(FILETOANALYSE)==1:
        equalpop()

    else:
        p = subprocess.Popen(get_cmd_line(args.python_version).split())
        p.wait()
        with open('Goodvibes_output.dat') as f:
            file = f.read()
        if 'Warning! Couldn\'t find frequency information' in file:
            equalpop()
        else:
            file = file.split('***************************************************************************************************************************************')
            for i in file[1].split('\n'):
                if i.strip():
                    j = i.split()
                    DF.loc[len(DF)] = [j[1], j[9], get_rvel(get_absolute_path(j[1])), get_wavelenght(get_absolute_path(j[1])), None]


def get_cmd_line(python):
    cmd = python + ' -m goodvibes ' + ' '.join(FILETOANALYSE) + ' --boltz'
    if args.solvation != '1': cmd += ' --freespace ' +  args.solvation
    return cmd


def gaussian(x, sigma, r, l):
    f = r/(sigma*np.sqrt(2*np.pi)) * np.exp(-0.5*((x-l)/sigma)**2)
    return f


def normalize(array):
    return array/np.max(np.abs(array)) * args.normalisation


def convolution():
    with alive_bar(len(list(DF.iterrows())), title='Greating plots') as bar:
        for index, row in DF.iterrows():
            conv = 0
            for l, r in zip(row['l'], row['R']):
                conv += gaussian(X, args.sigma, r, l)
            DF.at[index, 'conv'] = np.hstack((np.array(list([i] for i in X)), np.array(list([i] for i in normalize(conv)))))
            bar()
        
    if not (args.reference is None):
        shift(ref)
        

def weight_plot():
    conv = 0
    with alive_bar(len(list(DF.iterrows())), title='Weighting plots') as bar:
        for index, row in DF.iterrows():
            g = row['conv'][:, 1] * float(row['pop'])
            conv += g
            plt.plot(row['conv'][:, 0], g, color='gray', alpha=.3, label=(row['fln'].strip('.log').title()[:10]+'...') if len(args.file) > 1 else None)
            bar()
    plt.plot(row['conv'][:, 0], normalize(conv), color='salmon', label='Weigthed computational graph')


def get_reference(filename):
    if args.reference:
        data = np.loadtxt(filename)
        try:
            l_cut = np.argwhere(data[:, 0] == args.initial_lambda)[0, 0]
            u_cut = np.argwhere(data[:, 0] == args.final_lambda)[0, 0]
            data = data[l_cut:u_cut, :]
        except IndexError:
            print('No splice of the set for the file {0} is made'.format(filename))
        x, y = data[:, 0], data[:, 1]
        plt.plot(x, normalize(y), color='brown', label='Experimental Graph')
        return data 
    return None

def show_plot():
    plt.xlabel('Wavelenght [nm]')
    plt.ylabel('Δε')
    plt.title(args.title)
    plt.xlim([args.initial_lambda, args.final_lambda])
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -.125), fancybox=True, shadow=True, ncol=3)
    plt.tight_layout()
    fig = plt.gcf()
    if args.save:      
        with alive_bar(1, title='Saving plot') as bar:      
            with open('ecd.pickle', 'wb') as f:
                pickle.dump(fig, f)
            bar()
        # fig = pickle.load(open('ecd.pickle','rb'))    
    plt.show()


def x_max(ref, value=False):
    idx = args.skip_data
    if not value:
        try:
            max_ = np.max(ref[idx:, 1])
            min_ = np.min(ref[idx:, 1])
        except ValueError:
            raise Exception(f'{idx} is out of bound: there is no so many data in array.')
        set = max_ if max_ > -min_ else min_
        return (ref[np.argwhere(ref[:, 1] == set )[0], 0], 'pos' if set == max_ else 'neg')

    if value == 'pos':
        set = np.max(ref[idx:, 1])
    else:
        set = np.min(ref[idx:, 1])

    return ref[np.argwhere(ref[:, 1] == set )[0], 0]

        
def shift(ref):
    if not args.shift:
        shs = []
        x_ref, sign = x_max(ref)
        with alive_bar(len(list(DF.iterrows())), title='Shifting plots') as bar:
            for index, row in DF.iterrows():
                x_row = x_max(row['conv'], sign)
                shift = x_ref - x_row
                shs.append(shift)
                x_shifted = row['conv'][:, 0] + shift
                DF._set_value(index, 'conv', np.hstack((np.array(list([i] for i in x_shifted)), np.array(list([i] for i in row['conv'][:, 1])))))
                bar()
        print(f'Plots shifted on average {np.mean(shs)} nm')

    else:
        with alive_bar(len(list(DF.iterrows())), title='Shifting plots') as bar:
            for index, row in DF.iterrows():
                x_shifted = row['conv'][:, 0] + args.shift
                DF._set_value(index, 'conv', np.hstack((np.array(list([i] for i in x_shifted)), np.array(list([i] for i in row['conv'][:, 1])))))
                bar()
        print(f'Every plot is shifted by {args.shift} nm, as asked by the user.')






if __name__ == '__main__':
    with alive_bar(len(args.file), title='Splitting files') as bar:
        for i in args.file:
            split_file(i)
            bar()

    get_population()

    ref = get_reference(args.reference)
    convolution()
    weight_plot()
    show_plot()
