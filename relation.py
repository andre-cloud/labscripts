import matplotlib.pyplot as plt
import numpy as np


x = [24.1, 24.8, 25.8, 26.9, 21.5, 21.0]
y = [24.3, 25.1, 26.0, 25.6, 20.6, 23.1]
lx = np.linspace(21, 27)
n = ['A1', 'B2', 'C3', 'D4', 'E5', 'F6']

xp, yp = -0.2, 0.2 

def lin_reg(x):
    m = 0.739519350823562
    q = 6.365935366873560
    return m*x+q 


fig, ax = plt.subplots()


ax.scatter(x, y, color="#1e90ff")
for i, txt in enumerate(n):
    ax.annotate(txt, (x[i]+xp, y[i]+yp))

ax.plot(lx, list(map(lin_reg, lx)), '--', color='forestgreen')

ax.annotate('$R^2$=0.7391', (24, 23.5), )

ax.set_xlim(20, 28)
ax.set_ylim(20, 27)

plt.xlabel(r'$\Delta E_{d, NN}$ [kcal/mol]')
plt.ylabel(r'$\Delta G_{exp}$ [kcal/mol]')
plt.title('Correlation between $\Delta G_{exp}$ and $\Delta E_{d, NN}$')

plt.tight_layout()

plt.savefig('correlazione.png', dpi=300)
plt.show()