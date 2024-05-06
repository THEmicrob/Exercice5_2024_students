import numpy as np
import subprocess
import matplotlib.pyplot as plt
from matplotlib.animation import ArtistAnimation

repertoire = ''  # Path to the compiled code
executable = 'ex5.exe'  # Name of the compiled code
# executable = 'Exercice5.exe'  # Name of the compiled code
input_filename = 'propagation_constante\\b\\droite.in' # Input file

nsteps=5000
CFL=np.array([1.0,1.0007,1.0008])

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

# Analyse
for j in range(len(param)):
    output_file = f"propagation_constante/b/{paramstr}={param[j]}_{j}.out"
    lbl_f=output_file+'_f'
    lbl_v=output_file+'_v'
    lbl_x=output_file+'_x'
    data_f=np.loadtxt(lbl_f)
    data_v=np.loadtxt(lbl_v)
    data_x=np.loadtxt(lbl_x)
    
    plt.plot(data_x, np.delete(data_f[-1,:],0), label=r"$\beta=$"+str(param[j]), color=colors[j],zorder=k[j])
    print(f"nx={param[j]}")
plt.legend(fontsize=fs-3)
# artists = [ plt.plot(data_x, np.delete(data_f[i,:],0), label=f"t={i}") for i in range(0, nsteps, 100)]
plt.xlabel('x',fontsize=fs)
plt.ylabel('y',fontsize=fs)
# plot
# fig=plt.figure()
# image=ArtistAnimation(fig,artists, interval=100, blit=True)
# image.save(f"propagation_constante/b/{paramstr}={param[j]}_{j}.gif", writer='imagemagick')
plt.tight_layout()
plt.savefig(f"propagation_constante/b/{paramstr}={param[j]}_{j}.png",dpi=200)
plt.show()