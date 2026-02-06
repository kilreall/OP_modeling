import numpy as np
from scipy. integrate import odeint
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 22
plt.rcParams["font.family"] = "Century Gothic"
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['lines.linewidth'] = 1.5

def R(x): # F = 1
    gFx = gJS*((1*(1+1)-3/2*(3/2+1)+1/2*(1/2+1))/(2*1*(1+1)))
    Ex = gFx*x*uB*B
    det0_MB = det0_M + Ex/h
    return Y/2*i0_M/(1+i0_M+(2*det0_MB/Y)**2)

def Rl():
    return Y/2*i0_M/(1+i0_M+(2*det0_M/Y)**2)

def f(G, t):
    return M@G

def Num(P, mf, t):
    y = 1/3*P*(s*R(mf)*(abs(mf) != 0) + l*R(mf)*int(mf==0))
    N = np.zeros(len(t))
    for i in range(1, len(t)):
        N[i] = N[i-1] + (t[i] - t[i-1]) * (y[i] + y[i-1]) / 2
    return N


h = 1e-34
uB = 927.4*1e-26 # SI-26 SGS-23
gJS = 2.00233113 # СО неизвестно, скорее СИ
gJP = 1.3362 # СО не известна, скорее СИ
Y = 2*np.pi*6.0666*1e6
B = 20*1e-6*0
i0_M = 4/5*1
det0_M = 0
T = 0.00001
s, l = 0.98, 0.02

# sigma polar contribute
M0 = np.zeros((3,3))

M0[0,0] = -2*R(1)
M0[1,0] = 1*R(1)
M0[2,0] = 1*R(1)

M0[0,2] = 1*R(-1)
M0[1,2] = 1*R(-1)
M0[2,2] = -2*R(-1)

Ms = M0/9

# linear polar contribute
b = 1/9 # b в квадрате

L = np.zeros((3,3))
L[0,1] = b*Rl()
L[1,1] = -2*b*Rl()
L[2, 1] = b*Rl()

M = s*Ms + l*L

GP0 = np.array([1/3, 1/3, 1/3])
N = 1000
t = np.linspace(0, T, N)

Sol = odeint(f, GP0, t)
Sol = np.array(Sol)
 
G0 = Sol.transpose()[1]
G1 = Sol.transpose()[0]
G_1 = Sol.transpose()[2]


fig, ax = plt.subplots(figsize=(8, 8))
ax1 = ax.twinx()
ax.plot(t*1e6, G0, color="blue")

plt.xlabel("t, [мкc]")
ax.set_ylabel("Population |F=1, mF=0›", color="blue")


print(np.sum(Sol[-1]))
print(Sol[-1][1])
        

N1 = Num(G1, 1, t)
N_1 = Num(G_1, -1, t)
N0 = Num(G0, 0, t)
N = N0+N1+N_1


kTN = 1.2*1e-7
Tr = 362e-9

#ax1.plot(t*1e6, N*kTN*1e6, color="red")
ax1.plot(t*1e6, N*Tr/3*1e6, color="red")

ax1.set_ylabel("Heat, μK", color="red")
ax1.set_ylim(-0.08,1.5)



plt.show()