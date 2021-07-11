import matplotlib
matplotlib.use('PDF')
params = {'font.size': 10, 'font.family': 'Times New Roman', 'mathtext.fontset': 'stix'}
matplotlib.rcParams.update(params)
import matplotlib.pyplot as plt
import numpy as np

pl = 'pot.xvg'

fontsize=14
fig, ax = plt.subplots(figsize=(3.25,2.5))

c = np.loadtxt('C36.charmm-gui.gromacs/energy/' + pl, comments=['#', '@'])
s = np.loadtxt('C36.str2gmx.gromacs/energy/' + pl,    comments=['#', '@'])

ax.plot(c[:,0], c[:,1] / 1000, marker='o', ls='None', color='black', label='CHARMM-GUI')
ax.plot(s[:,0], s[:,1] / 1000, color='red', label='STR2GMX')

ax.set_xlabel('t (ps)', fontsize=fontsize)
ax.set_ylabel('U (MJ/mol)', fontsize=fontsize)
ax.legend(frameon=False)

fig.tight_layout()
fig.savefig('plot.pdf')

