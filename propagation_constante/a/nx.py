import numpy as np
import subprocess
import matplotlib.pyplot as plt
from matplotlib.animation import ArtistAnimation

repertoire = ''  # Path to the compiled code
executable = 'ex5.exe'  # Name of the compiled code
# executable = 'Exercice5.exe'  # Name of the compiled code
input_filename = 'propagation_constante\\a\\droite.in' # Input file

nsteps=5000
nx=np.linspace(0,1000,10)

paramstr = 'nx'  # Parameter name to scan
param=nx

# Simulations
# for j in range(len(param)):
#     output_file = f"propagation_constante/a/{paramstr}={param[j]}_{j}.out"
#     cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param[j]:.15g} output={output_file}"
#     print(cmd)
#     subprocess.run(cmd, shell=True)
#     print('Done.')

# output_file = f"propagation_constante/a/{paramstr}={param}.out"
# cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param:.15g} output={output_file}"
# print(cmd)
# subprocess.run(cmd, shell=True)
# print('Done.')

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'coral', 'lime', 'aqua', 'fuchsia', 'yellow', 'black', 'red', 'green', 'blue', 'purple', 'orange', 'brown', 'pink', 'gray', 'olive', 'cyan', 'navy', 'teal', 'violet', 'magenta', 'turquoise', 'tan', 'salmon', 'gold', 'indigo', 'maroon', 'ivory', 'chartreuse', 'wheat', 'azure', 'lavender', 'rose', 'lemon', 'scarlet', 'sand', 'cream', 'khaki', 'plum', 'orchid', 'lime', 'pear', 'cerulean', 'slate', 'steel', 'apricot', 'brass', 'bronze']

# Analyse
for j in range(len(param)):
    output_file = f"propagation_constante/a/{paramstr}={param[j]}_{j}.out"
    lbl_f=output_file+'_f'
    lbl_v=output_file+'_v'
    lbl_x=output_file+'_x'
    data_f=np.loadtxt(lbl_f)
    data_v=np.loadtxt(lbl_v)
    data_x=np.loadtxt(lbl_x)
    plt.plot(data_x, np.delete(data_f[2,:],0), label=f"nx={param[j]}", color=colors[j])
    print(f"nx={param[j]}")
plt.legend()
# artists = [ plt.plot(data_x, np.delete(data_f[i,:],0), label=f"t={i}") for i in range(0, nsteps, 100)]

# plot
fig=plt.figure()
# image=ArtistAnimation(fig, artists, interval=100, blit=True)
# image.save(f"propagation_constante/a/{paramstr}={param[j]}_{j}.gif", writer='imagemagick')
plt.show()