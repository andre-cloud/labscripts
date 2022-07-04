import argparse
import datetime
import sys

import numpy as np
from scipy.constants import Planck , R, k


KB = k
H = Planck
CONV = 4186 #kcal -> J

parser = argparse.ArgumentParser()

parser.add_argument('order', help='Define the reaction order', choices=[1, 2], type=int)
parser.add_argument('-T', help='Absolute temperature. Default %(default)s', type=float, default=298.15)
parser.add_argument('-G', '-g', '--gibbs_energy', help='Define Gibbs energy (kcal/mol) of the reaction TS to get the half-time, in kcal/mol', type=float, default=0)
parser.add_argument('-t', '--time', help='Define hal-life time to get the corrisponding ΔG', type=float, default=0)
parser.add_argument('-c', '--conc', help='Concentration of the substrate', default=0, type=float)
parser.add_argument('-j', '--kj', help='Get the energies in KJ/mol', action='store_true')


args = parser.parse_args()

text = '''{frame}
Results:
- Reaction order: {order}
- ΔG: {g} {unit}
- t1/2: {t}
- Temperature: {Temp} K
- k : {const}
{frame}'''

def eyring(g, T):
    return KB*T/H * np.exp(-(g*CONV)/(R*T))

def inverse_eyring(k, T):
    f = 1000 if args.kj else CONV
    return -R*T*np.log((k*H)/(KB*T)) / f

def get_time(g, T, order, c=None):
    if order == 1:
        return trasform_time(np.log(2)/eyring(g, T))
    if order == 2:
        return trasform_time(1/(eyring(g, T) * c))

def get_g_from_t(t, T, order, c=None):
    if order == 1:
        k = np.log(2)/t
        return inverse_eyring(k, T)
    if order == 2:
        k = (t*c)**-1
        return inverse_eyring(k, T)

def trasform_time(timestamp):
    return datetime.timedelta(seconds=timestamp)

if __name__=='__main__':
    if args.time != 0:
        output = get_g_from_t(args.time, args.T, args.order, args.conc)
        print(text.format(frame='='*20, order=args.order,g = output, unit= 'kJ/mol' if args.kj else 'kcal/mol', t = args.time, Temp=args.T, const=inverse_eyring(output, args.T)))
    if args.gibbs_energy != 0:
        output = get_time(args.gibbs_energy, args.T, args.order, args.conc)
        print(text.format(frame='='*20, order=args.order,g = args.gibbs_energy, unit= 'kJ/mol' if args.kj else 'kcal/mol', t = output, Temp=args.T, const=eyring(args.gibbs_energy, args.T)))