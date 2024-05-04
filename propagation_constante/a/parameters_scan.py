import numpy as np
import subprocess
import matplotlib.pyplot as plt

repertoire = ''  # Path to the compiled code
executable = 'Exercice5_students.exe'  # Name of the compiled code
input_filename = 'propagation_constante\\a\\input.in'  # Name of the input file

nsteps = 5000

paramstr = 'nsteps'  # Parameter name to scan
param = nsteps  # Parameter values to scan

# Simulations
output_file = f"propagation_constante/a/{paramstr}={param}.out"
cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param:.15g} output={output_file}"
print(cmd)
subprocess.run(cmd, shell=True)
print('Done.')

# Analyse
lbl_f=output_file+'_f'
lbl_v=output_file+'_v'
lbl_x=output_file+'_x'
data_f=np.loadtxt(lbl_f)
data_v=np.loadtxt(lbl_v)
data_x=np.loadtxt(lbl_x)

# vague part vers la droite
plt.figure()
plt.plot(data_x, data_v)
plt.xlabel('x')
plt.ylabel('v')
plt.title('v(x)')
plt.grid()
plt.savefig(f'propagation_constante/a/v_x_{paramstr}={param}.png')
plt.show()
