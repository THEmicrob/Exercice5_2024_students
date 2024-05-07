import numpy as np
import subprocess
import matplotlib.pyplot as plt
from matplotlib.animation import ArtistAnimation

repertoire = ''  # Path to the compiled code
executable = 'ex5.exe'  # Name of the compiled code
# executable = 'Exercice5.exe'  # Name of the compiled code
input_filename = 'vague\\a\\input.in' # Input file

nsteps=5000

paramstr = 'nsteps'  # Parameter name to scan
param=nsteps

# Simulations
output_file = f"vague/a/{paramstr}={param}.out"
cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param:.15g} output={output_file}"
print(cmd)
subprocess.run(cmd, shell=True)
print('Done.')

# output_file = f"vague/a/{paramstr}={param}.out"
# cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param:.15g} output={output_file}"
# print(cmd)
# subprocess.run(cmd, shell=True)
# print('Done.')

hL = 7000.0
hR = 200.0
hC = 35.0
h00= 3.0
xa= 3e5
xb= 7e5
xc= 7.2e5
xd= 8.5e5

colors = [ 'black','blue','coral', 'purple', 'orange', 'brown', 'pink', 'gray', 'olive', 'cyan', 'navy', 'teal', 'violet', 'magenta', 'turquoise', 'tan', 'salmon', 'gold', 'indigo', 'maroon', 'ivory', 'chartreuse', 'wheat', 'azure', 'lavender', 'rose', 'lemon', 'scarlet', 'sand', 'cream', 'khaki', 'plum', 'orchid', 'lime', 'pear', 'cerulean', 'slate', 'steel', 'apricot', 'brass', 'bronze']
fs=15

dt=[0.0259219,0.0287733]
k=[3,2,1]

# Analyse
output_file = f"vague/a/{paramstr}={param}.out"
lbl_f=output_file+'_f'
lbl_v=output_file+'_v'
lbl_x=output_file+'_x'
lbl_h=output_file+'_h'
data_f=np.loadtxt(lbl_f)
data_v=np.loadtxt(lbl_v)
data_x=np.loadtxt(lbl_x)
h=np.loadtxt(lbl_h)

# plt.plot(data_x, np.delete(data_f[-1,:],0), label=r"$\beta=$"+str(param))
# print(f"nx={param}")
# plt.legend(fontsize=fs-3)
# xb=np.linspace(xd)
plt.plot(data_x,-h)
artists = [ plt.plot(data_x, np.delete(data_f[i,:],0), label=f"t={i}") for i in range(0,int(nsteps/3), 100)]
plt.xlabel('x',fontsize=fs)
plt.ylabel('y',fontsize=fs)
# plot
fig=plt.figure()
image=ArtistAnimation(fig,artists, interval=100, blit=True)
image.save(f"vague/a/{paramstr}={param}.gif", writer='imagemagick')
# plt.tight_layout()
# plt.savefig(f"vague/a/{paramstr}={param}.png",dpi=200)
plt.show()