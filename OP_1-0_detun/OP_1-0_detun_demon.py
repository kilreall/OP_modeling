import numpy as np
from scipy. integrate import odeint
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 22
plt.rcParams["font.family"] = "Century Gothic"
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['lines.linewidth'] = 1.5


TK0 = np.load(r"OP_1-0_detun\T_d_0.npy")
P0 = np.load(r"OP_1-0_detun\P_d_0.npy")
TKY = np.load(r"OP_1-0_detun\T_d_Y.npy")
PY = np.load(r"OP_1-0_detun\P_d_Y.npy")
TK2Y = np.load(r"OP_1-0_detun\T_d_2Y.npy")
P2Y = np.load(r"OP_1-0_detun\P_d_2Y.npy")
TK3Y = np.load(r"OP_1-0_detun\T_d_3Y.npy")
P3Y = np.load(r"OP_1-0_detun\P_d_3Y.npy")


fig, ax = plt.subplots(figsize=(8, 8))


ax.plot(TK0, P0, color="b")
ax.plot(TKY, PY, color="r")
ax.plot(TK2Y, P2Y, color="black")
ax.plot(TK3Y, P3Y, color="g")

ax.set_xlabel("Heat, [μK]")
ax.set_ylabel("Population")
ax.set_xlim(-0.03, 0.55)
plt.show()