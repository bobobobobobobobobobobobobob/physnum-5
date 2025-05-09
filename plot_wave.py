import numpy as np
import matplotlib.pyplot as plt

# === Chargement des données ===
output_base = "test.out"

data_f = np.loadtxt(output_base + "_f")
times = data_f[:, 0]
fx = data_f[:, 1:]

data_en = np.loadtxt(output_base + "_en")
t_en = data_en[:, 0]
E = data_en[:, 1]

x = np.loadtxt(output_base + "_x")

# === Heatmap (f(x, t)) ===
plt.figure(figsize=(10, 6))
extent = [x[0], x[-1], times[-1], times[0]]  # y inversé pour que t=0 soit en haut
plt.imshow(fx, extent=extent, aspect='auto', cmap='inferno', origin='upper')
plt.colorbar(label="Amplitude f(x,t)")
plt.title("Amplitude de l'onde f(x,t) en fonction du temps et de l'espace")
plt.xlabel("x")
plt.ylabel("temps t")
plt.tight_layout()
plt.savefig("heatmap_b.svg", bbox_inches="tight")

# === Tracé de l'énergie E(t) ===
plt.figure(figsize=(8, 4))
plt.plot(t_en, E, label="E(t)")
plt.title("Énergie totale du système")
plt.xlabel("temps t")
plt.ylabel("E(t)")
plt.grid()
plt.tight_layout()

