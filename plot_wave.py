import numpy as np
import matplotlib.pyplot as plt


output_base = "test.out"


data_f = np.loadtxt(output_base + "_f")
times = data_f[:, 0]
fx = data_f[:, 1:]


data_en = np.loadtxt(output_base + "_en")
t_en = data_en[:, 0]
E = data_en[:, 1]


x = np.loadtxt(output_base + "_x")


plt.figure(figsize=(10, 5))
for i in range(0, len(times), max(1, len(times)//10)):
    plt.plot(x, fx[i], label=f't = {times[i]:.2f}')
plt.title("f(x,t) à différents instants")
plt.xlabel("x")
plt.ylabel("f(x)")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(t_en, E)
plt.title("Énergie E(t)")
plt.xlabel("temps t")
plt.ylabel("E(t)")
plt.grid()
plt.tight_layout()
plt.show()
