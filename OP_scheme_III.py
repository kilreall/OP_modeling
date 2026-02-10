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
    return Y/2*I2/(1+I2+(2*det2/Y)**2)

def R1(p): # F = 1
    I1 = i1*p
    return Y/2*I1/(1+I1+(2*det1/Y)**2)

def f(G, t):
    return M@G

bmfl = [1/15, 1/20, 0]
bmfp = [1/60, 0, 0, 1/20, 1/10]
bmfn = [1/60, 1/20, 1/10, 0, 0]
amfl = [0, 5/12]    
amfp = [5/12, 0, 5/12]
amfn = [5/12, 5/12, 0]

def Num(P, F, mf, t):
    if F == 2:
        y = P*(R2(l2)*bmfl[abs(mf)]+R2(sp2)*bmfp[mf]+R2(sn2)*bmfn[mf])
        N = np.zeros(len(t))
        for j in range(1, len(t)):
            N[j] = N[j-1] + (t[j] - t[j-1]) * (y[j] + y[j-1]) / 2
        return N
    elif F == 1:
        y = P*(R1(l1)*amfl[abs(mf)]+R1(sp1)*amfp[mf]+R1(sn1)*amfn[mf])
        N = np.zeros(len(t))
        for j in range(1, len(t)):
            N[j] = N[j-1] + (t[j] - t[j-1]) * (y[j] + y[j-1]) / 2
        return N


h = 1.054e-34
uB = 927.4*1e-26 # SI-26 SGS-23
gJS = 2.00233113 # СО неизвестно, думаю СИ
gJP = 1.3362 # СО не известна, думаю СИ
Y = 2*np.pi*6.0666*1e6
det1 = 0
det2 = 0
i2 = 1/10
i1 = 1/10
B = 0.5*1e-4*0 # Gauss 0.5

#Djkl = 0
#R = Y/2*i/(1+i+(2*Djkl/Y)**2)

cj = 0.007
sp1, l1, sn1 = cj, 1 - cj*2, cj
sp2, l2, sn2 = 1/2, 0, 1/2 

S0 = np.array([1/5, 1/5, 1/5, 1/5, 1/5, 0, 0, 0])
n = 1000
T = 0.0003

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
L2[:,3] = L2[:,3]*R2(l2)

# sigma+

P2 = np.zeros((8,8))

P2[0,2] = b01*b21
P2[1,2] = b01*b11
P2[2,2] = -b01*(b11 + b21 + a11 + a01)
P2[5,2] = b01*a11
P2[6,2] = b01*a01
P2[:,2] = P2[:,2]*R2(sp2)

P2[1,3] = b10*b10
P2[2,3] = b10*b00
P2[3,3] = -b10*(b00 + b10 + 2*a10)
P2[5,3] = b10*a10
P2[7,3] = b10*a10
P2[:,3] = P2[:,3]*R2(sp2)

P2[2,4] = b21*b01
P2[3,4] = b21*b11
P2[4,4] = -b21*(b11 + b01 + a11 + a01)
P2[6,4] = b21*a01
P2[7,4] = b21*a11
P2[:,4] = P2[:,4]*R2(sp2)

# sigma -

N2 = np.zeros((8,8))

N2[0] = P2[4]
N2[1] = P2[3]
N2[2] = P2[2]
N2[3] = P2[1]
N2[4] = P2[0]
N2[5] = P2[7]
N2[6] = P2[6]
N2[7] = P2[5]
N2 = N2.transpose() # это для удобства
N2_ = np.zeros((8,8))
N2_[0] = N2[4]
N2_[1] = N2[3]
N2_[2] = N2[2]
N2_[3] = N2[1]
N2_[4] = N2[0]
N2_[5] = N2[7]
N2_[6] = N2[6]
N2_[7] = N2[5]


N2 = N2_.transpose()

# F = 1

# linear

L1 = np.zeros((8,8))

L1[0,5] = a11*b21
L1[1,5] = a11*b11
L1[2,5] = a11*b01
L1[5,5] = -a11*(b21+b11+b01+a01)
L1[6,5] = a11*a01
L1[:,5] = L1[:,5]*R1(l1)

L1[2,7] = a11*b01
L1[3,7] = a11*b11
L1[4,7] = a11*b21
L1[6,7] = a11*a01
L1[7,7] = -a11*(b01+b11+b21+a01)
L1[:,7] = L1[:,7]*R1(l1)

#sigma +

P1 = np.zeros((8,8))

P1[0,6] = a01*b21
P1[1,6] = a01*b11
P1[2,6] = a01*b01
P1[5,6] = a01*a11
P1[6,6] = -a01*(b21+b11+b01+a11)
P1[:,6] = P1[:,6]*R1(sp1)

P1[1,7] = a10*b10
P1[2,7] = a10*b00
P1[3,7] = a10*b10
P1[5,7] = a10*a10
P1[7,7] = -a10*(2*b10+b00+a10)
P1[:,7] = P1[:,7]*R1(sp1)

# sigma -

N1 = np.zeros((8,8))

N1[0] = P1[4]
N1[1] = P1[3]
N1[2] = P1[2]
N1[3] = P1[1]
N1[4] = P1[0]
N1[5] = P1[7]
N1[6] = P1[6]
N1[7] = P1[5]
N1 = N1.transpose() # это для удобства
N1_ = np.zeros((8,8))
N1_[0] = N1[4]
N1_[1] = N1[3]
N1_[2] = N1[2]
N1_[3] = N1[1]
N1_[4] = N1[0]
N1_[5] = N1[7]
N1_[6] = N1[6]
N1_[7] = N1[5]

N1 = N1_.transpose()

M = L1 + P1 + L2 + P2 + N1 + N2

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
print(Sol.transpose()[-1][6])

N2 = Num(G2, 2, 2, t)
N1 = Num(G1, 2, 1, t)
N0 = Num(G0, 2, 0, t)
N_1 = Num(G_1, 2, -1, t)
N_2 = Num(G_2, 2, -2, t)
K1 = Num(H1, 1, 1, t)
K0 = Num(H0, 1, 0, t)
K_1 = Num(H_1, 1, -1, t)

N = N2+N1+N_1+N_2+K1+K0+K_1+N0


Tr = 362e-9
ax1.plot(tmk, N*Tr/3*1e6, color='red')
ax1.set_ylim(-0.03, 1.5)
#ax.plot(tmk, N2)
#ax.plot(tmk, N1)
#ax.plot(tmk, N_1)
#ax.plot(tmk, N_2)
#ax.plot(tmk, K1)
#ax.plot(tmk, K0)
#ax.plot(tmk, K_1)


ax1.set_ylabel("Heat, [μK]", color="red")
#ax.legend(loc='center right', fontsize='medium')
#plt.grid(True)
#plt.savefig('Ph_2_2.png', dpi=300, bbox_inches='tight')

plt.show()

