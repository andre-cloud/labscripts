import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 

parser = argparse.ArgumentParser()

parser.add_argument('energies', nargs='+')
parser.add_argument('-t','--title', default='DIM Analysis')

args = parser.parse_args()

data_ = list(np.array(args.energies, dtype=float))
key = ['NN', 'BOC1', 'BOC2', 'R3', 'R4', 'EI']
dict = {v:k for k,v in zip(data_, key)}

df = pd.DataFrame.from_dict(dict, orient='index', columns=['Fraction'])
ei = data_.pop(-1)

df['TS'] = 0
df['TS']['NN'] = sum(data_)+ei

ax = df[['Fraction', 'TS']].T.plot.bar(stacked=True, legend=True, width=0.4, align='center', color=['#b3e0ff', '#4db8ff', '#0099ff', '#006bb3', '#003d66', '#000f1a'])
ax.properties()['children'][1].set_color('#228b22')
ax.set_title(args.title)
ax.set_ylabel(r'$\Delta$E (kcal/mol)')
plt.axhline(y=0, xmin=0, xmax=2, color='black', alpha=0.7, linewidth=.75)

plt.xticks(rotation=0)
plt.show()
