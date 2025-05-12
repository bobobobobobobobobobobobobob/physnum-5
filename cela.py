import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.animation import FuncAnimation

# === PARAMÈTRES DE LA GRILLE ===
Nx, Ny = 100, 100
X = np.linspace(0, 1, Nx)
Y = np.linspace(0, 1, Ny)
X, Y = np.meshgrid(X, Y)

# === LIRE UN FICHIER wave_N.dat ===
def load_frame(filename):
    data = np.loadtxt(filename)
    Z = data[:,2].reshape((Nx, Ny))
    return Z

# === LISTE DES FICHIERS À CHARGER ===
frames = [load_frame(f"wave_{n}.dat") for n in range(0, 301, 10)]

# === FIGURE 3D ===
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')

surf = [ax.plot_surface(X, Y, frames[0], cmap=cm.inferno)]

def update(frame):
    ax.collections.clear()
    surf[0] = ax.plot_surface(X, Y, frame, cmap=cm.inferno)
    ax.set_zlim(-1, 1)  # adapte selon amplitude
    return surf[0],

ani = FuncAnimation(fig, update, frames=frames, interval=1000, blit=False)
plt.title("Propagation de l'onde en 2D")
# === GRAPHIQUE 3D SPATIO-TEMPOREL AVEC COULEUR ===
from mpl_toolkits.mplot3d import Axes3D  # déjà inclus par projection='3d'

# On construit les coordonnées x, y, t pour chaque frame
num_frames = len(frames)
T_vals = np.linspace(0, 1, num_frames)  # échelle du temps normalisée
X_, Y_, T_, F_ = [], [], [], []

for k, Z in enumerate(frames):
    t = T_vals[k]
    for i in range(Nx):
        for j in range(Ny):
            X_.append(X[i,j])
            Y_.append(Y[i,j])
            T_.append(t)
            F_.append(Z[i,j])

X_, Y_, T_, F_ = map(np.array, (X_, Y_, T_, F_))

# Affichage 3D avec couleur = amplitude
fig2 = plt.figure(figsize=(10,7))
ax2 = fig2.add_subplot(111, projection='3d')

p = ax2.scatter(X_, Y_, T_, c=F_, cmap='inferno', s=1, alpha=0.7)
fig2.colorbar(p, label="Amplitude f(x,y,t)")
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_zlabel("Temps t")
ax2.set_title("Évolution spatio-temporelle de l'amplitude")
plt.tight_layout()


time_indices = [0, len(frames)//4, len(frames)//2, 3*len(frames)//4, -1]
x_line = X[:,0]
colors = ['k', 'b', 'g', 'orange', 'red']

plt.figure()
for idx, color in zip(time_indices, colors):
    plt.plot(x_line, frames[idx][:,Ny//2], color=color, label=f"t = {idx}")
plt.xlabel("Position x")
plt.ylabel("Amplitude f")
plt.title("Profils d'onde à différents instants (y = milieu)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

