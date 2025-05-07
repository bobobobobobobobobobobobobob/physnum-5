import numpy as np
import matplotlib
#matplotlib.use("MacOSX")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import subprocess
from multiprocessing import Pool

plt.rcParams.update({
    #"text.usetex": True,
    "font.family": "sans-serif",
    "font.size": 16,
})
fs=16
ls=14

repertoire = "./"
executable = "physnum5"
input_filename = "configuration.in"
# === Paramètres ===
output_base = "test.out"



g  = 9.81
f_hat = 1.0
L=15
tfin=10
nx=20
n_init=4
CFL = 1.0
h00= 4.0
initialization="mode"


u = np.sqrt(h00*g)
omega_n = (0.5+n_init) * np.pi * u/L
Tn = 2*np.pi/omega_n
tfin = Tn
nsteps=int(u*tfin/(L/nx)+1)

paramstr="nsteps"
nsimul=150
param = np.arange(1, nsimul+1) * nsteps
param2 = np.arange(1, nsimul+1) * nx


def _simulate(idx):
    output = f"{output_base}_{idx}"
    cmd = f"{repertoire}{executable} {input_filename} impose_nsteps=true tfin={tfin} f_hat={f_hat} h00={h00} CFL={CFL} initialization={initialization} L={L} tfin={tfin} nsteps={param[idx]} nx={param2[idx]} n_init={n_init} output={output}"
    subprocess.run(cmd, shell=True)
    data_f = np.loadtxt(output + "_f")
    times = data_f[:, 0]
    fx = data_f[:, 1:]
    x = np.loadtxt(output + "_x")

    p_t = np.tile(times, (fx.shape[1], 1)).T
    p_x = np.tile(np.linspace(0, L, fx.shape[1]), (times.shape[0], 1))
    p_fx = f_hat * np.cos(omega_n * p_t) * np.cos(p_x * omega_n / u)
    
    dx = L/(fx.shape[1]-1)
    #print(f"{output}: CFL={u*(tfin/param[idx])/(L/param2[idx])} norm:{np.linalg.norm(fx[-1,:]-p_fx[-1,:])} dx:{dx}")
    error = np.sum(np.abs(fx[-1,:] - p_fx[-1,:]))*dx



    print(f'Done. ({paramstr}={param[idx]})')
    return error

errors = np.zeros(nsimul)

def plot():
    plt.figure()
    plt.scatter(np.arange(1, nsimul+1), errors)
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig("conv.svg", bbox_inches="tight")
    # plt.show()

if __name__ == '__main__':
    batch_size = 10
    for nbatch in range(0, nsimul, batch_size):
        with Pool(processes=batch_size) as pool:
            results = pool.map(_simulate, range(nbatch, min(nsimul, nbatch + batch_size)))
            for i, res in enumerate(results):
                errors[nbatch + i] = res
    plot()

