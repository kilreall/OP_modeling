import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 20
plt.rcParams["font.family"] = "Century Gothic"
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['lines.linewidth'] = 1.5
#rcParams["legend.frameon"] = False

def R2(x, y, p): # F = 2
    I2 = i2*p
    gFx = gJS*(2*(2+1)-3/2*(3/2+1)+1/2*(1/2+1))/(2*2*(2+1))
    gFy = gJP*(2*(2+1)-3/2*(3/2+1)+3/2*(3/2+1))/(2*2*(2+1)) 
    Ex = gFx*x*uB*B
    Ey = gFy*y*uB*B
    det = (Ex-Ey)/h + d02
    return Y/2*I2/(1+I2+(2*det/Y)**2)

def R1(x, y, p): # F = 1
    I1 = i1*p
    gFx = gJS*(1*(1+1)-3/2*(3/2+1)+1/2*(1/2+1))/(2*1*(1+1))
    gFy = gJP*(2*(2+1)-3/2*(3/2+1)+3/2*(3/2+1))/(2*2*(2+1)) 
    Ex = gFx*x*uB*B
    Ey = gFy*y*uB*B
    det = (Ex-Ey)/h + d01
    return Y/2*I1/(1+I1+(2*det/Y)**2)

def f(G, t):
    return M@G


    
bmfl = [0, 1/12, 1/3]
bmfp = [1/4, 1/6, 0, 1/6, 1/4]
bmfn = [1/4, 1/4, 1/6, 0, 1/6]
amfl = [1/3, 1/4]    
amfp = [1/4, 1/2, 1/12]
amfn = [1/4, 1/12, 1/2]

def Num(P, F, mf, t):
    if F == 2:
        y = P*(R2(mf, mf, sl2)*bmfl[abs(mf)]+R2(mf, mf+1, sp2)*bmfp[mf]+R2(mf, mf-1, sn2)*bmfn[mf])
        N = np.zeros(len(t))
        for j in range(1, len(t)):
            N[j] = N[j-1] + (t[j] - t[j-1]) * (y[j] + y[j-1]) / 2
        return N
    elif F == 1:
        y = P*(R1(mf, mf, sl1)*amfl[abs(mf)]+R1(mf, mf+1, sp1)*amfp[mf]+R1(mf, mf-1, sn1)*amfn[mf])
        N = np.zeros(len(t))
        for j in range(1, len(t)):
            N[j] = N[j-1] + (t[j] - t[j-1]) * (y[j] + y[j-1]) / 2
        return N
    
h = 1.054e-34
uB = 927.4*1e-26 # SI-26 SGS-23
gJS = 2.00233113 # СО неизвестно, думаю СИ
gJP = 1.3362 # СО не известна, думаю СИ
Y = 2*np.pi*6.0666*1e6 # это значение взято из методички, нужно взять более точное
d01 = 0
d02 = 0
i1 = 0.5
i2 = 0.07
B = 0.2*1e-4 # Gauss 0.2

#Djkl = 0
#R = Y/2*i/(1+i+(2*Djkl/Y)**2)

cj = 0.01
sp2, sl2, sn2 = cj, 1 - cj*2, cj
sp1, sl1, sn1 = cj, 1 - cj*2, cj

S0 = np.array([1/5, 1/5, 1/5, 1/5, 1/5, 0, 0, 0])
n = 10000
T = 0.00025



# pi 
# F = 2
L2 = np.zeros((8,8))

L2[0,0] = -32*R2(2,2, sl2)
L2[1,0] = 8*R2(2,2, sl2)
L2[5,0] = 24*R2(2,2, sl2)

L2[0,1] = 2*R2(1,1, sl2)
L2[1,1] = -11*R2(1,1, sl2)
L2[2,1] = 3*R2(1,1, sl2)
L2[5,1] = 3*R2(1,1, sl2)
L2[6,1] = 3*R2(1,1, sl2)

L2[2,3] = 3*R2(-1,-1, sl2)
L2[3,3] = -11*R2(-1,-1, sl2)
L2[4,3] = 2*R2(-1,-1, sl2)
L2[6,3] = 3*R2(-1,-1, sl2)
L2[7,3] = 3*R2(-1,-1, sl2)

L2[3,4] = 8*R2(-2,-2, sl2)
L2[4,4] = -32*R2(-2,-2, sl2)
L2[7,4] = 24*R2(-2,-2, sl2)

L2 = L2/144
#F = 1
L1 = np.zeros((8,8))
L1[0,5] = 6*R1(1,1, sl1)
L1[1,5] = 3*R1(1,1, sl1)
L1[2,5] = 9*R1(1,1, sl1)
L1[5,5] = -27*R1(1,1, sl1)
L1[6,5] = 9*R1(1,1, sl1)

L1[1,6] = 12*R1(0,0, sl1)
L1[3,6] = 12*R1(0,0, sl1)
L1[5,6] = 4*R1(0,0, sl1)
L1[6,6] = -32*R1(0,0, sl1)
L1[7,6] = 4*R1(0,0, sl1)

L1[2,7] = 9*R1(-1,1, sl1)
L1[3,7] = 3*R1(-1,1, sl1)
L1[4,7] = 6*R1(-1,1, sl1)
L1[6,7] = 9*R1(-1,1, sl1)
L1[7,7] = -27*R1(-1,1, sl1)

L1 = L1/144

# coef a and b
b01 = 1/4
b11 = 1/12
b12 = 1/6
b22 = 1/3
a00 = 1/3
a01 = 1/4
a10 = 1/12
a11 = 1/4
a12 = 1/2

# sigma +
# F = 2
P2 = np.zeros((8,8))

P2[0, 1] = b12*b22*R2(1,2, sp2)
P2[1, 1] = -b12*(b22+a12)*R2(1, 2, sp2)  
P2[5, 1] =  b12*a12*R2(1, 2, sp2)            

P2[0, 2] =  b01*b12*R2(0, 1, sp2)
P2[1, 2] =  b01*b11*R2(0, 1, sp2)
P2[2, 2] =  -b01*(b11+b12+a01+a11)*R2(0, 1, sp2)
P2[5, 2] =  b01*a11*R2(0, 1, sp2)
P2[6, 2] =  b01*a01*R2(0, 1, sp2)

P2[2, 3] =  b01*b01*R2(-1, 0, sp2)
P2[3, 3] =  -b01*(b01+a10+a00+a10)*R2(-1, 0, sp2)
P2[5, 3] =  b01*a10*R2(-1, 0, sp2)
P2[6, 3] =  b01*a00*R2(-1, 0, sp2)
P2[7, 3] =  b01*a10*R2(-1, 0, sp2)

P2[1, 4] =  b12*b01*R2(-2, -1, sp2)
P2[3, 4] =  b12*b11*R2(-2, -1, sp2)
P2[4, 4] =  -b12*(b11+b01+a11+a01)*R2(-2, -1, sp2)
P2[6, 4] =  b12*a01*R2(-2, -1, sp2)
P2[7, 4] =  b12*a11*R2(-2, -1, sp2)

# F = 1
P1 = np.zeros((8,8))
P1[0, 5] =  a12*b22*R1(1, 2, sp1)
P1[1, 5] =  a12*b12*R1(1, 2, sp1)
P1[5, 5] =  -a12*(b22+b12)*R1(1, 2, sp1)

P1[0, 6] =  a01*b12*R1(0, 1, sp1)
P1[1, 6] =  a01*b11*R1(0, 1, sp1)
P1[2, 6] =  a01*b01*R1(0, 1, sp1)
P1[5, 6] =  a01*a11*R1(0, 1, sp1)
P1[6, 6] =  -a01*(b01+b11+b12+a11)*R1(0, 1, sp1)

P1[1, 7] =  a10*b01*R1(-1, 0, sp1)
P1[3, 7] =  a10*b01*R1(-1, 0, sp1)
P1[5, 7] =  a10*a10*R1(-1, 0, sp1)
P1[6, 7] =  a10*a00*R1(-1, 0, sp1)
P1[7, 7] =  -a10*(b01+b01+a10+a00)*R1(-1, 0, sp1)

#print(P)

# sigma -

N = np.zeros((8,8))

N[0, 1] = b12*b22*R2(-1,-2, sn2)
N[1, 1] = -b12*(b22+a12)*R2(-1, -2, sn2)  
N[5, 1] =  b12*a12*R2(-1, -2, sn2)             

N[0, 2] =  b01*b12*R2(0, -1, sn2)
N[1, 2] =  b01*b11*R2(0, -1, sn2)
N[2, 2] =  -b01*(b11+b12+a01+a11)*R2(0, -1, sn2)
N[5, 2] =  b01*a11*R2(0, -1, sn2)
N[6, 2] =  b01*a01*R2(0, -1, sn2)

N[2, 3] =  b01*b01*R2(1, 0, sn2)
N[3, 3] =  -b01*(b01+a10+a00+a10)*R2(1, 0, sn2)
N[5, 3] =  b01*a10*R2(1, 0, sn2)
N[6, 3] =  b01*a00*R2(1, 0, sn2)
N[7, 3] =  b01*a10*R2(1, 0, sn2)

N[1, 4] =  b12*b01*R2(2, 1, sn2)
N[3, 4] =  b12*b11*R2(2, 1, sn2)
N[4, 4] =  -b12*(b11+b01+a11+a01)*R2(2, 1, sn2)
N[6, 4] =  b12*a01*R2(2, 1, sn2)
N[7, 4] =  b12*a11*R2(2, 1, sn2)

N[0, 5] =  a12*b22*R1(-1, -2, sn1)
N[1, 5] =  a12*b12*R1(-1, -2, sn1)
N[5, 5] =  -a12*(b22+b12)*R1(-1, -2, sn1)

N[0, 6] =  a01*b12*R1(0, -1, sn1)
N[1, 6] =  a01*b11*R1(0, -1, sn1)
N[2, 6] =  a01*b01*R1(0, -1, sn1)
N[5, 6] =  a01*a11*R1(0, -1, sn1)
N[6, 6] =  -a01*(b01+b11+b12+a11)*R1(0, -1, sn1)

N[1, 7] =  a10*b01*R1(1, 0, sn1)
N[3, 7] =  a10*b01*R1(1, 0, sn1)
N[5, 7] =  a10*a10*R1(1, 0, sn1)
N[6, 7] =  a10*a00*R1(1, 0, sn1)
N[7, 7] =  -a10*(b01+b01+a10+a00)*R1(1, 0, sn1)

# тут как бы порядок нарушен, его нужно возвратить
N1 = np.zeros((8,8))
N1[0] = N[4]
N1[1] = N[3]
N1[2] = N[2]
N1[3] = N[1]
N1[4] = N[0]
N1[5] = N[7]
N1[6] = N[6]
N1[7] = N[5]
N1 = N1.transpose() # это для удобства
N2 = np.zeros((8,8))
N2[0] = N1[4]
N2[1] = N1[3]
N2[2] = N1[2]
N2[3] = N1[1]
N2[4] = N1[0]
N2[5] = N1[7]
N2[6] = N1[6]
N2[7] = N1[5]

N = N2.transpose()
#print(N)

N1 = np.zeros((8,8))
N2 = np.zeros((8,8))
N = N.transpose() # столбцы в строки
N1[0] = N[0]
N1[1] = N[1]
N1[2] = N[2]
N1[3] = N[3]
N1[4] = N[4]
N1 = N1.transpose()
N2[5] = N[5]
N2[6] = N[6]
N2[7]= N[7]
N2 = N2.transpose()


M = P1+P2+L1+L2+N1 + N2


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
ax.plot(tmk, G0, color="blue")
# np.save(r"pictures\OPIIt.npy", tmk)
# np.save(r"pictures\OPIIP.npy", G0)

print(np.sum(Sol.transpose()[-1]))
print(Sol.transpose()[-1][2])


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
# np.save(r"pictures\OPIITK.npy", N*Tr/3*1e6)
ax1.set_ylim(-0.15, 3)
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
plt.savefig('OP_II.svg', dpi=300, bbox_inches='tight')

plt.show()
