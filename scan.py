import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def lines(array, x, color, x_pad, y_pad, pad=0):
    
    plt.vlines(x=x[np.argwhere(array[:len(array)//2] == np.max(array[:len(array)//2]))[0]], ymin=0, ymax=1.5, colors=color)
    plt.vlines(x=x[np.argwhere(array[:len(array)//2] == np.max(array[:len(array)//2]))[0]], ymin=np.max(array[:len(array)//2])-1.5+pad, ymax=np.max(array[:len(array)//2])+1.5+pad, colors=color)

    plt.vlines(x=x[np.argwhere(array[len(array)//2:] == np.max(array[len(array)//2:]))[0]+len(array)//2], ymin=0, ymax=1.5, colors=color)
    plt.vlines(x=x[np.argwhere(array[len(array)//2:] == np.max(array[len(array)//2:]))[0]+len(array)//2], ymin=np.max(array[len(array)//2:])-1.5+pad, ymax=np.max(array[len(array)//2:])+1.5+pad, colors=color)

    plt.text(x=x[np.argwhere(array[:len(array)//2] == np.max(array[:len(array)//2]))[0]]+x_pad, y=np.max(array[:len(array)//2])+y_pad+pad, s=f'∆E = {round(float(np.max(array[:len(array)//2])), 2)}\n{round(float(x[np.argwhere(array[:len(array)//2] == np.max(array[:len(array)//2]))[0]]), 2)}°', color=color)

    plt.text(x=x[np.argwhere(array[len(array)//2:] == np.max(array[len(array)//2:]))[0]+len(array)//2]+x_pad, y=np.max(array[len(array)//2:])+y_pad+pad, s=f'∆E = {round(float(np.max(array[len(array)//2:])), 2)}\n{round(float(x[np.argwhere(array[len(array)//2:] == np.max(array[len(array)//2:]))[0]+len(array)//2]), 2)}°', color=color)




data = np.loadtxt('nn_boc2_bn2/scan.txt')
x1 = data[:, 0]
b3lyp = data[:, 1]
x2 = data[:, 2] 
m062x = data[:, 3] 
x3 = data[:, 4]
wB97xd = data[:, 5] 
# x = data[:, 0]
# b3lyp = data[:, 1]
# m06 = data[:, 2]
# wb97x = data[:, 3]


plt.scatter(x1, b3lyp, lw=.5, color='salmon', alpha=1, label='B3LYP(BJ)')
plt.scatter(x2, m062x + 10, lw=.5, color='brown', alpha=1, label='M06-2x')
plt.scatter(x3, wB97xd +20 , lw=.5, color='darkgrey', alpha=1, label='ωB97x-D')

lines(b3lyp, x1, 'salmon', 20, -4, 0)
lines(m062x, x2, 'brown', 15, 0, 10)
lines(wB97xd, x3, 'darkgrey', 20, 1, 20)


plt.legend(loc='upper center', bbox_to_anchor=(0.5, -.125), fancybox=True, shadow=True, ncol=3)

axes = plt.gca()
axes.set_xlim([-125, 300])
axes.set_ylim([0, 70])

plt.title('Scan of BnN(Boc)N(Boc)Bn')
plt.ylabel('∆E [kcal/mol]')

plt.tight_layout()
# plt.show()
plt.savefig('nn_boc2_bn2/scan.png', dpi=700, bbox_inches = 'tight')