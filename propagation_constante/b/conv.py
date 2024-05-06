import numpy as np
import subprocess
import matplotlib.pyplot as plt
from matplotlib.animation import ArtistAnimation

repertoire = ''  # Path to the compiled code
executable = 'ex5.exe'  # Name of the compiled code
# executable = 'Exercice5.exe'  # Name of the compiled code
input_filename = 'propagation_constante\\b\\droite.in' # Input file

nsteps=5000
CFL=np.geomspace(0.999,1.001,50)

paramstr = 'CFL'  # Parameter name to scan
param=CFL

# Simulations
for j in range(len(param)):
    output_file = f"propagation_constante/b/{paramstr}={param[j]}_{j}.out"
    cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param[j]:.15g} output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')

# output_file = f"propagation_constante/b/{paramstr}={param}.out"
# cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param:.15g} output={output_file}"
# print(cmd)
# subprocess.run(cmd, shell=True)
# print('Done.')

colors = [ 'black','blue','coral', 'purple', 'orange', 'brown', 'pink', 'gray', 'olive', 'cyan', 'navy', 'teal', 'violet', 'magenta', 'turquoise', 'tan', 'salmon', 'gold', 'indigo', 'maroon', 'ivory', 'chartreuse', 'wheat', 'azure', 'lavender', 'rose', 'lemon', 'scarlet', 'sand', 'cream', 'khaki', 'plum', 'orchid', 'lime', 'pear', 'cerulean', 'slate', 'steel', 'apricot', 'brass', 'bronze']
fs=15

dt=[0.0259219,0.0287733]
k=[3,2,1]
y_f=np.empty(np.shape(CFL))
lengths=np.empty(np.shape(CFL))

# Analyse
for j in range(len(param)):
    output_file = f"propagation_constante/b/{paramstr}={param[j]}_{j}.out"
    lbl_f=output_file+'_f'
    lbl_v=output_file+'_v'
    lbl_x=output_file+'_x'
    data_f=np.loadtxt(lbl_f)
    data_v=np.loadtxt(lbl_v)
    data_x=np.loadtxt(lbl_x)
    lengths[j]=len(data_f[-1,:])
    y_f[j]=data_f[-1,80]
print(lengths)
plt.plot(CFL,y_f,'+-')
# plt.xscale('log')
plt.legend(fontsize=fs-3)
plt.xlabel(r'$\beta$',fontsize=fs)
plt.ylabel(r'$f(t_{\text{fin}},x=4)$',fontsize=fs)
plt.tight_layout()
plt.savefig(f"propagation_constante/b/{paramstr}={param[j]}_{j}_conv.png",dpi=200)
plt.show()