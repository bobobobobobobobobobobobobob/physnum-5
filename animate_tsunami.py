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
L=1000 * 1000
tfin=12000
nx=2000
n_init=4
CFL = 0.99
h00= 4.0
initialization="autre"
equation_type="C"

hL=8000.0
hR=20.0
xa=450e3
xb=950e3
def h0(x_):
    if(type(x_) is float):
        x_ = np.array([x_])
        print(x_)
    out = np.zeros_like(x_)
    out[x_ <= xa] = hL
    out[x_ >= xb] = hR
    for i in range(out.size):
        if not xa < x_[i] < xb:
            continue
        out[i] = 0.5*(hL+hR)+0.5*(hL-hR)*np.cos(np.pi*(x_[i]-xa)/(xb-xa))
    return out

#u = np.sqrt(h00*g)
#omega_n = (0.5+n_init) * np.pi * u/L
#Tn = 2*np.pi/omega_n

#tfin = 


cmd = f"{repertoire}{executable} {input_filename} equation_type={equation_type} v_uniform=false tfin={tfin} f_hat={f_hat} h00={h00} CFL={CFL} initialization={initialization} L={L} tfin={tfin} x1=50000 x2=350000 nx={nx} n_init={n_init} hL={hL} hR={hR} xa={xa} xb={xb} output={output_base}"
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
#p_t = np.tile(times, (fx.shape[1], 1)).T
#p_x = np.tile(np.linspace(0, L, nx+1), (times.shape[0], 1))
#print(p_t.shape, times.shape, fx.shape[1], p_x.shape)
#p_fx = f_hat * np.cos(omega_n * p_t) * np.cos(p_x * omega_n / u)
#p_line, = ax.plot(x, p_fx[0], lw=2, label="Analytique")

time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

ax.set_title("Propagation de l'onde : f(x,t)")
ax.set_xlabel("x")
ax.set_ylabel("f(x)")
ax.grid(True)
ax.set_ylim(np.min(fx)*1.1, np.max(fx)*1.1)

def update(frame):
    line.set_ydata(fx[frame])
    #p_line.set_ydata(p_fx[frame])
    time_text.set_text(f"t = {times[frame]:.3f} s")
    #return line, p_line, time_text
    return line, time_text

ani = animation.FuncAnimation(fig, update, frames=len(times), interval=15, blit=True)
plt.tight_layout()
#ani.save("tsunami.gif", fps=30)


wave_height = np.max(fx, axis=0)
wave_x = x
wave_t = np.zeros_like(x)
#mask = (0.2 <= wave_x) & (wave_x < (1-0.03)*1e6)
k=2*2
fit_range=50
for i in np.arange(0,x.size):
    t_idx = np.argmax(fx[:,i])
    mask = np.arange(max(0, t_idx-fit_range),min(fx.shape[0],t_idx+fit_range+1))
    f_fit = fx[mask,i]
    t_fit = times[mask]
    fit = np.polyfit(t_fit, f_fit, 2)
    wave_t[i] = -fit[1]/(2*fit[0])
    #wave_t[i] = np.mean(t_fit)
    if i == 1000:
        plt.figure()
        plt.plot(t_fit, f_fit)
        plt.plot(t_fit, fit[0]*(t_fit**2) + fit[1]*(t_fit) + fit[2], label="fit")
        plt.legend()



#print(wave_t[x.size//2: x.size//2 + 10])
plt.figure()
plt.scatter(wave_x, wave_height, s=4,label="Simulation")
# WKB
x_init_idx = np.argmax(wave_x>0.2*1e6)
x_init = float(wave_x[x_init_idx])
u = np.sqrt(g*h0(wave_x))
if(equation_type=="B"):
    wkb_A = (np.sqrt(u)[x_init_idx] * fx[0,x_init_idx] / np.sqrt(u))
elif(equation_type=="A"):
    wkb_A =  (fx[0,x_init_idx]/np.sqrt(u)[x_init_idx])*np.sqrt(u)
elif(equation_type=="C"):
    wkb_A = ((u[x_init_idx]**1.5) * fx[0,x_init_idx]) / (u**1.5)
plt.scatter(wave_x, wkb_A,s=4,label="Solution WKB")
plt.xlabel("Position [m]")
plt.ylabel("Hauteur [m]")
plt.legend()
plt.savefig("wave_height.svg", bbox_inches="tight")

diff_x = wave_x[k:] - wave_x[:-k]
diff_t = wave_t[k:] - wave_t[:-k]
v = diff_x / diff_t

#print((wave_x[k:-k]).size, v.size)

plt.figure()
pltx = wave_x[int(k/2):-int(k/2)]
mask = (0.2*1e6 < pltx)
plt.scatter(pltx[mask], v[mask], s=4,label="Simulation")
plt.scatter(pltx[mask], np.sqrt(g*h0(pltx[mask])),s=4,label="Solution WKB") 
plt.xlabel("Position [m]")
plt.ylabel("Vitesse [m/s]")
plt.legend()
plt.savefig("wave_spd.svg", bbox_inches="tight")


plt.figure()
#mp = plt.imshow(fx, cmap='hot', extent=(0,L,tfin,0),norm=matplotlib.colors.LogNorm())
mp = plt.imshow(fx[::-1,], cmap='hot', extent=(0,L,0,tfin), aspect="auto")
plt.colorbar(mp, location='right', label="Hauteur [m]")
plt.xlabel("Position [m]")
plt.ylabel("Temps [s]")
plt.savefig("tsunami_hm.svg", bbox_inches="tight")



# Pour CFL>1
if CFL>1:
    x_idx = 128
    plt.figure()
    plt.plot(times, fx[:,x_idx])
    plt.xlabel("Temps [s]")
    plt.ylabel("Hauteur [m]")
    plt.yscale('log')
    plt.savefig(f"CFL={CFL}_x={x_idx*L/nx}_fuck.svg", bbox_inches="tight")

# Question 5.3c)
if initialization=="mode":
    pass


plt.show()

