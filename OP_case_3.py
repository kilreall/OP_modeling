import numpy as np
from scipy. integrate import odeint
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 20
plt.rcParams["font.family"] = "Century Gothic"
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['lines.linewidth'] = 1.5
#rcParams["legend.frameon"] = False

def R2(x, y): # F = 2
    gFx = gJS*(2*(2+1)-3/2*(3/2+1)+1/2*(1/2+1))/(2*2*(2+1))
    gFy = gJP*(2*(2+1)-3/2*(3/2+1)+3/2*(3/2+1))/(2*2*(2+1)) 
    Ex = gFx*x*uB*B
    Ey = gFy*y*uB*B
    det = (Ex-Ey)/h
    return Y/2*i/(1+i+(2*det/Y)**2)

def R1(x, y): # F = 1
    gFx = gJS*(1*(1+1)-3/2*(3/2+1)+1/2*(1/2+1))/(2*1*(1+1))
    gFy = gJP*(2*(2+1)-3/2*(3/2+1)+3/2*(3/2+1))/(2*2*(2+1)) 
    Ex = gFx*x*uB*B
    Ey = gFy*y*uB*B
    det = (Ex-Ey)/h
    return Y/2*i/(1+i+(2*det/Y)**2)*0

def f(G, t):
    return M@G
    
h = 1e-34
uB = 927.4*1e-26 # SI-26 SGS-23
gJS = 2.00233113 # СО неизвестно, думаю СИ
gJP = 1.3362 # СО не известна, думаю СИ
Y = 2*np.pi*6.0666*1e6 # это значение взято из методички, нужно взять более точное
i = 1/5
B = 0.5*1e-4 # Gauss 0.5

#Djkl = 0
#R = Y/2*i/(1+i+(2*Djkl/Y)**2)

sp, sl, sn = 0, 1, 0

S0 = np.array([1/5, 1/5, 1/5, 1/5, 1/5, 0, 0, 0])
n = 100
T = 0.0001



# pi 
L = np.zeros((8,8))

L[0,0] = -32*R2(2,2)
L[1,0] = 8*R2(2,2)
L[5,0] = 24*R2(2,2)

L[0,1] = 2*R2(1,1)
L[1,1] = -11*R2(1,1)
L[2,1] = 3*R2(1,1)
L[5,1] = 3*R2(1,1)
L[6,1] = 3*R2(1,1)

L[2,3] = 3*R2(-1,-1)
L[3,3] = -11*R2(-1,-1)
L[4,3] = 2*R2(-1,-1)
L[6,3] = 3*R2(-1,-1)
L[7,3] = 3*R2(-1,-1)

L[3,4] = 8*R2(-2,-2)
L[4,4] = -32*R2(-2,-2)
L[7,4] = 24*R2(-2,-2)

L[0,5] = 6*R1(1,1)
L[1,5] = 3*R1(1,1)
L[2,5] = 9*R1(1,1)
L[5,5] = -27*R1(1,1)
L[6,5] = 9*R1(1,1)

L[1,6] = 12*R1(0,0)
L[3,6] = 12*R1(0,0)
L[5,6] = 4*R1(0,0)
L[6,6] = -32*R1(0,0)
L[7,6] = 4*R1(0,0)

L[2,7] = 9*R1(-1,1)
L[3,7] = 3*R1(-1,1)
L[4,7] = 6*R1(-1,1)
L[6,7] = 9*R1(-1,1)
L[7,7] = -27*R1(-1,1)

#print(L)
L = L/144

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
P = np.zeros((8,8))

P[0, 1] = b12*b22*R2(1,2)
P[1, 1] = -b12*(b22+a01+a11)*R2(1, 2)  
P[5, 1] =  b12*a12*R2(1, 2)            

P[0, 2] =  b01*b12*R2(0, 1)
P[1, 2] =  b01*b11*R2(0, 1)
P[2, 2] =  -b01*(b11+b12+a01+a11)*R2(0, 1)
P[5, 2] =  b01*a11*R2(0, 1)
P[6, 2] =  b01*a01*R2(0, 1)

P[2, 3] =  b01*b01*R2(-1, 0)
P[3, 3] =  -b01*(b01+a10+a00+a10)*R2(-1, 0)
P[5, 3] =  b01*a10*R2(-1, 0)
P[6, 3] =  b01*a00*R2(-1, 0)
P[7, 3] =  b01*a10*R2(-1, 0)

P[1, 4] =  b12*b01*R2(-2, -1)
P[3, 4] =  b12*b11*R2(-2, -1)
P[4, 4] =  -b12*(b11+b01+a11+a01)*R2(-2, -1)
P[6, 4] =  b12*a01*R2(-2, -1)
P[7, 4] =  b12*a11*R2(-2, -1)

P[0, 5] =  a12*b22*R1(1, 2)
P[1, 5] =  a12*b12*R1(1, 2)
P[5, 5] =  -a12*(b22+b12)*R1(1, 2)

P[0, 6] =  a01*b12*R1(0, 1)
P[1, 6] =  a01*b11*R1(0, 1)
P[2, 6] =  a01*b01*R1(0, 1)
P[5, 6] =  a01*a11*R1(0, 1)
P[6, 6] =  -a01*(b01+b11+b12+a11)*R1(0, 1)

P[1, 7] =  a10*b01*R1(-1, 0)
P[3, 7] =  a10*b01*R1(-1, 0)
P[5, 7] =  a10*a10*R1(-1, 0)
P[6, 7] =  a10*a00*R1(-1, 0)
P[7, 7] =  -a10*(b01+b01+a10+a00)*R1(-1, 0)

#print(P)

# sigma -

N = np.zeros((8,8))

N[0, 1] = b12*b22*R2(-1,-2)
N[1, 1] = -b12*(b22+a01+a11)*R2(-1, -2)  
N[5, 1] =  b12*a12*R2(-1, -2)             

N[0, 2] =  b01*b12*R2(0, -1)
N[1, 2] =  b01*b11*R2(0, -1)
N[2, 2] =  -b01*(b11+b12+a01+a11)*R2(0, -1)
N[5, 2] =  b01*a11*R2(0, -1)
N[6, 2] =  b01*a01*R2(0, -1)

N[2, 3] =  b01*b01*R2(1, 0)
N[3, 3] =  -b01*(b01+a10+a00+a10)*R2(1, 0)
N[5, 3] =  b01*a10*R2(1, 0)
N[6, 3] =  b01*a00*R2(1, 0)
N[7, 3] =  b01*a10*R2(1, 0)

N[1, 4] =  b12*b01*R2(2, 1)
N[3, 4] =  b12*b11*R2(2, 1)
N[4, 4] =  -b12*(b11+b01+a11+a01)*R2(2, 1)
N[6, 4] =  b12*a01*R2(2, 1)
N[7, 4] =  b12*a11*R2(2, 1)

N[0, 5] =  a12*b22*R1(-1, -2)
N[1, 5] =  a12*b12*R1(-1, -2)
N[5, 5] =  -a12*(b22+b12)*R1(-1, -2)

N[0, 6] =  a01*b12*R1(0, -1)
N[1, 6] =  a01*b11*R1(0, -1)
N[2, 6] =  a01*b01*R1(0, -1)
N[5, 6] =  a01*a11*R1(0, -1)
N[6, 6] =  -a01*(b01+b11+b12+a11)*R1(0, -1)

N[1, 7] =  a10*b01*R1(1, 0)
N[3, 7] =  a10*b01*R1(1, 0)
N[5, 7] =  a10*a10*R1(1, 0)
N[6, 7] =  a10*a00*R1(1, 0)
N[7, 7] =  -a10*(b01+b01+a10+a00)*R1(1, 0)

N1 = np.zeros((8,8))
N1[0] = N[4]
N1[1] = N[3]
N1[2] = N[2]
N1[3] = N[1]
N1[4] = N[0]
N1[5] = N[7]
N1[6] = N[6]
N1[7] = N[5]
N1 = N1.transpose()
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

M = sp*P+sl*L+sn*N


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
plt.title("Population components time evolution")
plt.xlabel("time")
plt.ylabel("Population")
plt.plot(tmk, G0)
plt.plot(tmk, G2)
plt.plot(tmk, G1)
plt.plot(tmk, G_1)
plt.plot(tmk, G_2)
plt.plot(tmk, H1)
plt.plot(tmk, H0)
plt.plot(tmk, H_1)
#plt.plot(t, G_1)
#plt.plot(t, G_2)
#plt.show()
print(np.sum(Sol.transpose()[-1]))
print(Sol.transpose()[-1][2])

def Num(P, F, mf, t):
    if F == 2:
        y = P*(sl*R2(mf, mf)+sp*R2(mf, mf+1)*(int(abs(mf)<2))+sn*R2(mf, mf-1)*(int(abs(mf)<2)))
        N = np.zeros(len(t))
        for i in range(1, len(t)):
            N[i] = N[i-1] + (t[i] - t[i-1]) * (y[i] + y[i-1]) / 2
        return N
    elif F == 1:
        y = P*(sl*R1(mf, mf)+sp*R1(mf, mf+1)+sn*R1(mf, mf-1))
        N = np.zeros(len(t))
        for i in range(1, len(t)):
            N[i] = N[i-1] + (t[i] - t[i-1]) * (y[i] + y[i-1]) / 2
        return N



N2 = Num(G2, 2, 2, t)
N1 = Num(G1, 2, 1, t)
N_1 = Num(G_1, 2, -1, t)
N_2 = Num(G_2, 2, -2, t)
K1 = Num(H1, 1, 1, t)
K0 = Num(H1, 1, 0, t)
K_1 = Num(H1, 1, -1, t)

N = N2+N1+N_1+N_2+K1+K0+K_1

fig, ax = plt.subplots(figsize=(11, 8))

ax.plot(tmk, N)
#ax.plot(tmk, N2)
#ax.plot(tmk, N1)
#ax.plot(tmk, N_1)
#ax.plot(tmk, N_2)
#ax.plot(tmk, K1)
#ax.plot(tmk, K0)
#ax.plot(tmk, K_1)

plt.xlabel("t, [мкс]")
plt.ylabel("число фотонов на атом")
#ax.legend(loc='center right', fontsize='medium')
plt.grid(True)
plt.savefig('Ph_2_2.png', dpi=300, bbox_inches='tight')
