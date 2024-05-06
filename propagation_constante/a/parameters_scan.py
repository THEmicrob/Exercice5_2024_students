import numpy as np
import subprocess
import matplotlib.pyplot as plt
from matplotlib.animation import ArtistAnimation

repertoire = ''  # Path to the compiled code
executable = 'ex5.exe'  # Name of the compiled code
# executable = 'Exercice5.exe'  # Name of the compiled code
input_filename = ['propagation_constante\\a\\droite.in', 'propagation_constante\\a\\gauche.in','propagation_constante\\a\\repos.in']  # Input file

nsteps=5000
# nx=78

paramstr = 'nsteps'  # Parameter name to scan
param = nsteps  # Parameter values to scan

# Simulations
for i in range(len(input_filename)):
    output_file = f"propagation_constante/a/{paramstr}={param}_{i}.out"
    cmd = f"{repertoire}{executable} {input_filename[i]} {paramstr}={param:.15g} output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')

# output_file = f"propagation_constante/a/{paramstr}={param}.out"
# cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param:.15g} output={output_file}"
# print(cmd)
# subprocess.run(cmd, shell=True)
# print('Done.')

# Analyse
for i in range(len(input_filename)):
    output_file = f"propagation_constante/a/{paramstr}={param}_{i}.out"
    lbl_f=output_file+'_f'
    lbl_v=output_file+'_v'
    lbl_x=output_file+'_x'
    data_f=np.loadtxt(lbl_f)
    data_v=np.loadtxt(lbl_v)
    data_x=np.loadtxt(lbl_x)

    # artists are f at different times
    # artists = [ plt.plot(data_x, np.delete(data_f[i,:],0), label=f"t={i}") for i in range(0, nsteps, 100)]
    artists = [ plt.plot(data_x, np.delete(data_f[i,:],0), label=f"t={i}") for i in range(0, nsteps, 100)]

    # plot
    fig=plt.figure()
    image=ArtistAnimation(fig, artists, interval=100, blit=True)
    image.save(f"propagation_constante/a/{paramstr}={param}_{i}.gif", writer='imagemagick')
    plt.show()
# lbl_f=output_file+'_f'
# lbl_v=output_file+'_v'
# lbl_x=output_file+'_x'
# data_f=np.loadtxt(lbl_f)
# data_v=np.loadtxt(lbl_v)
# data_x=np.loadtxt(lbl_x)

# # artists are f at different times
# # artists = [ plt.plot(data_x, np.delete(data_f[i,:],0), label=f"t={i}") for i in range(0, nsteps, 100)]
# artists = [ plt.plot(data_x, np.delete(data_f[i,:],0), label=f"t={i}") for i in range(0, nsteps, 100)]

# # plot
# fig=plt.figure()
# image=ArtistAnimation(fig, artists, interval=50, blit=True)
# image.save(f"propagation_constante/a/{paramstr}={param}.gif", writer='imagemagick')
# plt.show()