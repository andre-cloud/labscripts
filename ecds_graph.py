import argparse
from copyreg import pickle
import os
import re
import shutil
import subprocess
import sys


from alive_progress import alive_bar
import cclib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle


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
parser.add_argument('-gd','--graph_directory', help='Define the directory in which you want to save the files of the graph. Default: %(default)s', default='ecd_graphs')

parser.add_argument('--save', help='Save pickle and csvs of the graph', action='store_true')
parser.add_argument('--compare', help='Get the graph with the comparison of more functional over the reference. Suggested to give as input file generated with this program. Name for the graph will take filname-dft.txt', action='store_true')


args = parser.parse_args()

FILETOANALYSE = []
X = np.linspace(args.initial_lambda-100, args.final_lambda+100, (args.final_lambda-args.initial_lambda)*100**args.definition)

columns = ['fln', 'pop', 'R', 'l', 'conv', 't']
DF = pd.DataFrame(columns=columns)


class InputError(Exception):
    def __init__(self, message):
        super().__init__(message)

if args.pop:
    if not 0.99 < sum(args.pop) < 1.01:
        print(f'Sum population = {sum(args.pop)*100}%. {"Not complete population give" if sum(args.pop) < 0.99 else "Too much population given"}. Prociding with unputted populations')


def get_thoery(file):
    repeated_theory = 0
    with open(file) as f:
        data = f.readlines()
    level, bs = 'none', 'none'

    for line in data:
        if line.strip().find('External calculation') > -1:
            level, bs = 'ext', 'ext'
            break
        if '\\Freq\\' in line.strip() and repeated_theory == 0:
            try:
                level, bs = (line.strip().split("\\")[4:6])
                repeated_theory = 1
            except IndexError:
                pass
        elif '|Freq|' in line.strip() and repeated_theory == 0:
            try:
                level, bs = (line.strip().split("|")[4:6])
                repeated_theory = 1
            except IndexError:
                pass
        if '\\SP\\' in line.strip() and repeated_theory == 0:
            try:
                level, bs = (line.strip().split("\\")[4:6])
                repeated_theory = 1
            except IndexError:
                pass
        elif '|SP|' in line.strip() and repeated_theory == 0:
            try:
                level, bs = (line.strip().split("|")[4:6])
                repeated_theory = 1
            except IndexError:
                pass
        if 'DLPNO BASED TRIPLES CORRECTION' in line.strip():
            level = 'DLPNO-CCSD(T)'
        if 'Estimated CBS total energy' in line.strip():
            try:
                bs = ("Extrapol." + line.strip().split()[4])
            except IndexError:
                pass
        # Remove the restricted R or unrestricted U label
        if level[0] in ('R', 'U'):
            level = level[1:]
        # Remuve the TD FC from level
        level = level.split('TD')[0]
        # level = level.strip('')
    level_of_theory = '_'.join([level, bs])
    return level_of_theory

    # with open(filename) as f:
    #     fl = f.read()
    # fl = fl.split('----------------------------------------------------------------------')[1]
    # th = fl.strip().split('/')[0].split(')')[1].strip()
    # if not th:
    #     th = 'NA'
    # return th
    


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


def df_row(filename, pop):
    fl = os.path.split(filename)[1].strip('.log')
    return [fl, pop , get_rvel(get_absolute_path(fl)), get_wavelenght(get_absolute_path(fl)), None, get_thoery(get_absolute_path(fl)) ]


def equalpop():
    with alive_bar(len(FILETOANALYSE), title='Getting λ R(vel) and population for files') as bar:
        for i in FILETOANALYSE:
            DF.loc[len(DF)] = df_row(i.strip('.log'), 1/len(FILETOANALYSE))
            bar()

def get_population():

    if args.pop:
        with alive_bar(len(FILETOANALYSE), title='Getting λ R(vel) and population for files') as bar:
            for idx, i in enumerate(FILETOANALYSE):
                DF.loc[len(DF)] = df_row(i, args.pop[idx])
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
                    DF.loc[len(DF)] = df_row(j[1], j[9])


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
            plt.plot(row['conv'][:, 0], g, alpha=.3, label=(row['fln'].strip('.log').title()[:5]+'...-'+row['t']) if len(args.file) > 1 else None)
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
            pass
        x, y = data[:, 0], data[:, 1]
        plt.plot(x, normalize(y), color='brown', label='Experimental Graph')
        return data 
    return None

def show_plot(compare=False):
    plt.xlabel('Wavelenght [nm]')
    plt.ylabel('Δε')
    title = args.title
    if title == 'ECD graph' and len(set(DF['t'])) == 1:
        title += ' '+DF['t'][0]
    plt.title(title)
    plt.xlim([args.initial_lambda, args.final_lambda])
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -.125), fancybox=True, shadow=True, ncol=1)
    plt.tight_layout()
    fig = plt.gcf()

    if args.save and not compare:

        directory = args.graph_directory   

        if os.path.exists(directory):
            if 'y' in input(f'A directory named {directory} already exists. Existing directory  will be deleted, wanna procede? [y/n]').lower():
                shutil.rmtree(directory)   
            else:
                sys.exit()
        os.mkdir(directory)

        with alive_bar(1+len(list(DF.iterrows())), title='Saving plot') as bar:      
            with open(os.path.join(directory, 'ecd.pickle'), 'wb') as f:
                pickle.dump(fig, f)
            bar()
            for index, row in DF.iterrows():

                np.savetxt(os.path.join(directory, f"{row['fln'].strip('.log').title()}-graph-{row['t']}.txt"), row['conv'], newline='\n')
                bar()
            
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
        print(f'Plots shifted on average by {np.mean(shs)} nm')

    else:
        with alive_bar(len(list(DF.iterrows())), title='Shifting plots') as bar:
            for index, row in DF.iterrows():
                x_shifted = row['conv'][:, 0] + args.shift
                DF._set_value(index, 'conv', np.hstack((np.array(list([i] for i in x_shifted)), np.array(list([i] for i in row['conv'][:, 1])))))
                bar()
        print(f'Every plot is shifted by {args.shift} nm, as asked by the user.')



def compare_graphs():
    gf = {}
    with alive_bar(len(args.file), title='Loading files') as bar:
        for i in args.file:
            dft = i.split('-')[-1].split('.')[0] if  i.split('-')[-2] != 'cam' else 'cam-'+ i.split('-')[-1].split('.')[0]
            if dft in gf:
                print(f"\33[1m\33[33mWarnig\33[0m: {dft} is already present in the comparison. Check file {i}")
            gf[dft] = np.loadtxt(i)
            bar()
    with alive_bar(len(args.file), title='Plotting files') as bar:
        for i in gf:
            plt.plot(gf[i][:, 0], normalize(gf[i][:, 1]), label=i, alpha=.5)
            bar()
    
    if args.reference:
        with alive_bar(1, title='Plotting files') as bar:
            ref = get_reference(args.reference)
            bar()
    
    show_plot(compare=True)



if __name__ == '__main__':


    if args.compare:
        compare_graphs()
        sys.exit()
    
    with alive_bar(len(args.file), title='Splitting files') as bar:
        for i in args.file:
            split_file(i)
            bar()

    get_population()

    ref = get_reference(args.reference)
    convolution()
    weight_plot()
    show_plot()
