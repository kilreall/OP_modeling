import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

plt.rcParams['font.size'] = 22
plt.rcParams["font.family"] = "Century Gothic"
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['lines.linewidth'] = 1.5

t22 = np.load(r"pictures\OP22t.npy")
P22 = np.load(r"pictures\OP22P.npy")
TK22 = np.load(r"pictures\OP22TK.npy")

t10 = np.load(r"pictures\OP10t.npy")
P10 = np.load(r"pictures\OP10P.npy")
TK10 = np.load(r"pictures\OP10TK.npy")

tII = np.load(r"pictures\OPIIt.npy")
PII = np.load(r"pictures\OPIIP.npy")
TKII = np.load(r"pictures\OPIITK.npy")

tIII = np.load(r"pictures\OPIIIt.npy")
PIII = np.load(r"pictures\OPIIIP.npy")
TKIII = np.load(r"pictures\OPIIITK.npy")


fig, axex = plt.subplots(1, 4, figsize=(28,7))
plt.subplots_adjust(wspace=0.7, hspace=0.0)  # Увеличиваем промежутки
blackconst = np.zeros(1000)+1
blackty = np.linspace(0, 1, 1000)

dx = axex[0]
dx.set_box_aspect(1)
dx1 = dx.twinx()
dx1.set_box_aspect(1)

dx.plot(t22, P22, color="blue")





dx.text(350, 0.92, "P = 1")

dx.set_xlabel("time, [μs]")
dx.set_ylabel("Population |F=1›", color="blue")
dx.set_title("(a)                                                          ")
dx.text(320, 0.25, "Params:\nIp = 1/5\nδω = 0Г")

dx1.plot(t22, TK22, color="red")

dx1.text(270, 0.09, "ΔT = 0.25 μK")

dx1.set_ylabel("Heat, [μK]", color="red")
dx1.set_ylim(-0.08,1.5)


ax = axex[1]
ax.set_box_aspect(1)
ax1 = ax.twinx()
ax1.set_box_aspect(1)

ax.plot(t10, P10, color="blue")

t10x = 17.61
ax.plot(t10x*blackconst, blackty, color="black")

P10y = 0.9527
blackPx = np.linspace(0, t10x, 1000)
blackPy = P10y*blackconst
ax.plot(blackPx, blackPy, color="black")
ax.text(0.42, 0.964, "P = 0.95")

ax.set_xlabel("time, [μs]")
ax.set_ylabel("Population |F=1, mF=0›", color="blue")
ax.set_title("(b)                                                          ")
ax.text(5, 0.25, "Params:\nIp = 1/5\nδω = 0")

ax1.plot(t10, TK10, color="red")

TK10y = 0.3624
blackTKx = np.linspace(t10x, t10[-1], 1000)
blackTKy = blackconst*TK10y
ax1.plot(blackTKx, blackTKy, color="black")
ax1.text(17.95, 0.256, "ΔT = 0.36 μK")

ax1.set_ylabel("Heat, [μK]", color="red")
ax1.set_ylim(-0.08,1.5)

bx = axex[2]
bx.set_box_aspect(1)
bx1 = bx.twinx()
bx1.set_box_aspect(1)


bx.plot(tII, PII, color="blue")

tIIx = 148.1
bx.plot(tIIx*blackconst, blackty, color="black")

PIIy = 0.8758
blackPx = np.linspace(0, tIIx, 1000)
blackPy = PIIy*blackconst
bx.plot(blackPx, blackPy, color="black")
bx.text(0.42, 0.885, "P = 0.88")

bx.set_xlabel("time, [μs]")
bx.set_ylabel("Population |F=2, mF=0›", color="blue")
bx.set_title("(c)                                                          ")
bx.text(57, 0.25, "Params:\nIp = 1/10\nδω = 0")

bx1.plot(tII, TKII, color="red")

TKIIy = 1.1944
blackTKx = np.linspace(tIIx, tII[-1], 1000)
blackTKy = blackconst*TKIIy
bx1.plot(blackTKx, blackTKy, color="black")
bx1.text(150.4, 1.088, "ΔT = 1.2 μK")

bx1.set_ylabel("Heat, [μK]", color="red")
bx1.set_ylim(-0.08,1.5)

cx = axex[3]
cx.set_box_aspect(1)
cx1 = cx.twinx()
cx1.set_box_aspect(1)


cx.plot(tIII, PIII, color="blue")

tIIIx = 105.2
cx.plot(tIIIx*blackconst, blackty, color="black")

PIIIy = 0.8702
blackPx = np.linspace(0, tIIx, 1000)
blackPy = PIIIy*blackconst
cx.plot(blackPx, blackPy, color="black")
cx.text(0.42, 0.881, "P = 0.87")

cx.set_xlabel("time, [μs]")
cx.set_ylabel("Population |F=1, mF=0›", color="blue")
cx.set_title("(d)                                                          ")
cx.text(43, 0.15, "Params:\nIp = 1/5\nδω = 0")

cx1.plot(tIII, TKIII, color="red")

TKIIIy = 1.1026
blackTKx = np.linspace(tIIIx, tIII[-1], 1000)
blackTKy = blackconst*TKIIIy
cx1.plot(blackTKx, blackTKy, color="black")
cx1.text(115.4, 1, "ΔT = 1.1 μK")

cx1.set_ylabel("Heat, [μK]", color="red")
cx1.set_ylim(-0.08,1.5)

plt.tight_layout()

plt.savefig(r"pictures\OP_model.eps", format="eps")
plt.show()
