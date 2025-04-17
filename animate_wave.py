import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# === Paramètres ===
output_base = "test.out"

# === Charger les données ===
data_f = np.loadtxt(output_base + "_f")
times = data_f[:, 0]
fx = data_f[:, 1:]

x = np.loadtxt(output_base + "_x")

# === Préparer l'animation ===
fig, ax = plt.subplots(figsize=(10, 5))
line, = ax.plot(x, fx[0], lw=2)
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

ax.set_title("Propagation de l'onde : f(x,t)")
ax.set_xlabel("x")
ax.set_ylabel("f(x)")
ax.grid(True)
ax.set_ylim(np.min(fx)*1.1, np.max(fx)*1.1)

def update(frame):
    line.set_ydata(fx[frame])
    time_text.set_text(f"t = {times[frame]:.3f} s")
    return line, time_text

ani = animation.FuncAnimation(fig, update, frames=len(times), interval=50, blit=True)
plt.tight_layout()
plt.show()

