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

g  = 9.81;
f_hat = 1.0
L=15
tfin=10
nx=350
n_init=4
CFL = 1.0
h00= 4.0
initialization="mode"

u = np.sqrt(h00*g)
omega_n = (0.5+n_init) * np.pi * u/L
Tn = 2*np.pi/omega_n

tfin = Tn


cmd = f"{repertoire}{executable} {input_filename} tfin={tfin} f_hat={f_hat} h00={h00} CFL={CFL} initialization={initialization} L={L} tfin={tfin} nx={nx} n_init={n_init} output={output_base}"
subprocess.run(cmd, shell=True)


# === Charger les données ===
data_f = np.loadtxt(output_base + "_f")
times = data_f[:, 0]
fx = data_f[:, 1:]

x = np.loadtxt(output_base + "_x")

# === Préparer l'animation ===
fig, ax = plt.subplots(figsize=(10, 5))
line, = ax.plot(x, fx[0], lw=2, label="Simul")
# mode propre
p_t = np.tile(times, (fx.shape[1], 1)).T
p_x = np.tile(np.linspace(0, L, nx+1), (times.shape[0], 1))
#print(p_t.shape, times.shape, fx.shape[1], p_x.shape)
p_fx = f_hat * np.cos(omega_n * p_t) * np.cos(p_x * omega_n / u)
p_line, = ax.plot(x, p_fx[0], lw=2, label="Analytique")

time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

ax.set_title("Propagation de l'onde : f(x,t)")
ax.set_xlabel("x")
ax.set_ylabel("f(x)")
ax.grid(True)
ax.set_ylim(np.min(fx)*1.1, np.max(fx)*1.1)

def update(frame):
    line.set_ydata(fx[frame])
    p_line.set_ydata(p_fx[frame])
    time_text.set_text(f"t = {times[frame]:.3f} s")
    return line, p_line, time_text

ani = animation.FuncAnimation(fig, update, frames=len(times), interval=15, blit=True)
plt.tight_layout()


plt.figure()
mp = plt.imshow(fx, cmap='hot', extent=(0,L,tfin,0))
plt.colorbar(mp, location='right', label="Hauteur [m]")
plt.xlabel("Temps [s]")
plt.ylabel("Position [m]")
plt.savefig("heatmap.svg")

# Pour CFL>1
if CFL>1:
    x_idx = 100
    plt.figure()
    plt.plot(times, fx[:,x_idx])
    plt.xlabel("Temps [s]")
    plt.ylabel("Hauteur [m]")
    plt.yscale('log')
    plt.savefig(f"CFL={CFL}_x={x_idx*L/nx}.svg")

# Question 5.3c)
if initialization=="mode":
    pass


plt.show()

