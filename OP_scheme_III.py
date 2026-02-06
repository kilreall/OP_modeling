import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 20
plt.rcParams["font.family"] = "Century Gothic"
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['lines.linewidth'] = 1.5
#rcParams["legend.frameon"] = False

def R2(p): # F = 2
    I2 = i2*p
    return Y/2*I2/(1+I2+(2*det/Y)**2)

def R1(p): # F = 1
    I1 = i1*p
    return Y/2*I1/(1+I1+(2*det/Y)**2)

def f(G, t):
    return M@G

h = 1.054e-34
uB = 927.4*1e-26 # SI-26 SGS-23
gJS = 2.00233113 # СО неизвестно, думаю СИ
gJP = 1.3362 # СО не известна, думаю СИ
Y = 2*np.pi*6.0666*1e6 # это значение взято из методички, нужно взять более точное
det = 0
i2 = 1/10
i1 = 1/10
B = 0.5*1e-4*0 # Gauss 0.5

#Djkl = 0
#R = Y/2*i/(1+i+(2*Djkl/Y)**2)

sp1, l1, sn1 = 0, 1, 0
sp2, l2, sn2 = 1/3, 1/3, 1/3 

S0 = np.array([1/5, 1/5, 1/5, 1/5, 1/5, 0, 0, 0])
n = 100
T = 0.00015

# b - F = 2, a - F = 1

b11 = 1/20
b00 = 1/15
b01 = 1/60
b10 = 1/20
b21 = 1/10
a11 = 5/12
a10 = 5/12
a01 = 5/12


# F = 2

# linear

L2 = np.zeros((8,8))

L2[0,1] = b11*b21
L2[1,1] = -b11*(b21+b01+a11+a01)
L2[2,1] = b11*b01
L2[5,1] = b11*a11
L2[6,1] = b11*a01
L2[:,1] = L2[:,1]*R2(l2)

L2[1,2] = b00*b10
L2[2,2] = -b00*(2*b10 + 2*a10)
L2[3,2] = b00*b10
L2[5,2] = b00*a10
L2[7,2] = b00*a10
L2[:,2] = L2[:,2]*R2(l2)

L2[2,3] = b11*b01
L2[3,3] = -b11*(b21+b01+a11+a01)
L2[4,3] = b11*b21
L2[6,3] = b11*a01
L2[7,3] = b11*a11
L2[:,3] = L2[:,1]*R2(l2)

# sigma+

P2 = np.zeros((8,8))

P2[0,2] = b01*b21
P2[1,2] = b01*b11
P2[2,2] = -b01*(b11 + b21 + a11 + a01)
P2[5,2] = b01*a11
P2[6,2] = b01*a01
P2[:,2] = P2[:,1]*R2(sp2)

P2[1,3] = b01*b10
P2[2,3] = b01*b00
P2[3,3] = -b10*(b00 + b10 + 2*a10)
P2[5,3] = b01*a10
P2[7,3] = b01*a10
P2[:,2] = P2[:,1]*R2(sp2)

P2[2,4] = b21*b01
P2[3,4] = b21*b11
P2[4,4] = -b21*(b11 + b01 + a11 + a01)
P2[6,4] = b21*a01
P2[7,4] = b21*a11
P2[:,2] = P2[:,1]*R2(sp2)

# sigma -

N2 = np.zeros((8,8))

N2[0] = P2[4]
N2[1] = P2[3]
N2[2] = P2[2]
N2[3] = P2[1]
N2[4] = P2[0]
N2 = N2.transpose() # это для удобства
N2_ = np.zeros((8,8))
N2_[0] = N2[4]
N2_[1] = N2[3]
N2_[2] = N2[2]
N2_[3] = N2[1]
N2_[4] = N2[0]

N2 = N2_.transpose()

# F = 1

# linear

L1 = np.zeros((8,8))

L1[0,5] = a11*b21
L1[1,5] = a11*b11
L1[2,5] = a11*b01
L1[5,5] = -a11*(b21+b11+b01+a01)
L1[6,5] = a11*a01
L1[:,5] = L1[:,5]*R1(l2)

L1[2,7] = a11*b01
L1[3,7] = a11*b11
L1[4,7] = a11*b21
L1[6,7] = a11*a01
L1[7,7] = a11*(b01+b11+b21+a01)
L1[:,7] = L1[:,7]*R1(l2)

#sigma +

P1 = np.zeros((8,8))

P1[0,6] = a01*b21
P1[1,6] = a01*b11
P1[2,6] = a01*b01
P1[5,6] = a01*a11
P1[2,6] = -a01*(b21+b11+b01+a11)
P1[:,2] = P1[:,2]*R1(sp1)

P1[1,7] = a10*b10
P1[2,7] = a10*b00
P1[3,7] = a10*b10
P1[5,7] = a10*a10
P1[1,7] = -a10*(2*b10+b00+a10)
P1[:,2] = P1[:,2]*R1(sp1)

# sigma -

N1 = np.zeros((8,8))


N1[5] = P1[7]
N1[6] = P1[6]
N1[7] = P1[5]
N1 = N2.transpose() # это для удобства
N1_ = np.zeros((8,8))
N1_[5] = N1[7]
N1_[6] = N1[6]
N1_[7] = N1[5]

N1 = N1_.transpose()

M = L1 + P1 + N1 + L2 + P2 + N2

t = np.linspace(0, T, n)

Sol = odeint(f, S0, t)
Sol = np.array(Sol)
#print(Sol)
Sol = Sol.transpose()
G0 = Sol[2]
G2 = Sol[0]
G1 = Sol[1]
G_1 = Sol[3]
G_2 = Sol[4]
H1 = Sol[5]
H0 = Sol[6]
H_1 = Sol[7]
#print(G_0)

tmk = t*1e6

fig, ax = plt.subplots(figsize=(8, 8))
ax1 = ax.twinx()
ax.set_xlabel("time, [μs]")
ax.set_ylabel("Population", color='blue')
# plt.plot(tmk, G0)
# plt.plot(tmk, G2)
# plt.plot(tmk, G1)
# plt.plot(tmk, G_1)
# plt.plot(tmk, G_2)
# plt.plot(tmk, H1)
# plt.plot(tmk, H0)
# plt.plot(tmk, H_1)
ax.plot(tmk, H0, color="blue")

print(np.sum(Sol.transpose()[-1]))
print(Sol.transpose()[-1][2])


plt.show()

