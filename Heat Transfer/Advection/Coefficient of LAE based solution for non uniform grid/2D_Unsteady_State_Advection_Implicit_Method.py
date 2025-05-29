import numpy as np
import matplotlib.pyplot as plt

#step 1 User Input
L1=1
L2=1
rho=1000
cp=4180
imax=22
jmax=22
T0=0.5
T_eb=0
T_nb=1
u=-1
v=-1
beta=1.2
epsilon_st=0.0001
Dt=0.5
print("SELECT THE ADVECTION SCHEME(1/2/3)")
print("1.FOU\n")
print("2.SOU\n")
print("3.QUICK\n")
scheme = int(input("Enter the Advection Scheme:"))

#step 2 Grid generation
xi=np.linspace(0,1,imax-1)

#clustering
beta_p1=beta+1
beta_m1=beta-1
beta_p1_div_m1=(beta_p1/beta_m1)**(2*xi-1)
num=beta_p1*(beta_p1_div_m1)-beta_m1
den=1+beta_p1_div_m1
x=L1*(num/den)
y=L2*(num/den)

#define parameters
xc=np.zeros(imax)
yc=np.zeros(jmax)
Dx=np.zeros(imax)
Dy=np.zeros(jmax)

xc[1:imax-1]=0.5*(x[1:imax-1]+x[0:imax-2])
xc[0]=x[0]
xc[imax-1]=x[imax-2]
yc[1:jmax-1]=0.5*(y[1:jmax-1]+y[0:jmax-2])
yc[0]=y[0]
yc[jmax-1]=y[jmax-2]


XC,YC=np.meshgrid(xc,yc)
X,Y=np.meshgrid(x,y)

Dx[1:imax -1]=x[1:imax -1]-x[0:imax -2] 
Dx[0]=0
Dx[imax-1]=0
Dy=Dx

#step 3 weight calcualtion

def weights(scheme, Ds_D, Ds_U, Ds_UU):
    w = np.zeros(3)
    if scheme == 1:
        w[1] = 1
    elif scheme == 2:
        w[1] = (2*Ds_U + Ds_UU) / (Ds_U + Ds_UU)
        w[2] = -Ds_U / (Ds_U + Ds_UU)
    elif scheme == 3:
        w[0] = Ds_U * (2*Ds_U + Ds_UU) / ((Ds_D + Ds_U)*(Ds_D + 2*Ds_U + Ds_UU))
        w[1] = Ds_D * (2*Ds_U + Ds_UU) / ((Ds_D + Ds_U)*(Ds_U + Ds_UU))
        w[2] = -Ds_D * Ds_U / ((Ds_U + Ds_UU)*(Ds_D + 2*Ds_U + Ds_UU))
    return w

# Weights arrays
we_p = np.zeros((imax, 3))
we_m = np.zeros((imax, 3))
wn_p = np.zeros((jmax, 3))
wn_m = np.zeros((jmax, 3))

for i in range(1, imax-2):
    we_p[i] = weights(scheme, Dx[i+1], Dx[i], Dx[i-1])
    we_m[i] = weights(scheme, Dx[i], Dx[i+1], Dx[i+2])

for j in range(1, jmax-2):
    wn_p[j] = weights(scheme, Dy[j+1], Dy[j], Dy[j-1])
    wn_m[j] = weights(scheme, Dy[j], Dy[j+1], Dy[j+2])

#step 4 -IC and BCs
T = np.zeros((imax, jmax))
T[1:imax-1, 1:jmax-1] = T0
T[-1, :] = T_eb
T[:, -1] = T_nb

#step 5 - Mass flux
mx = rho * u
my = rho * v
mx_p = max(mx, 0)
mx_m = min(mx, 0)
my_p = max(my, 0)
my_m = min(my, 0)

#step 6 - coefficient of implicit LAE's
aE = np.zeros((imax, jmax))
aN = np.zeros((imax, jmax))
aP0 = np.zeros((imax, jmax))
aP = np.zeros((imax, jmax))

for j in range(1, jmax-1):
    for i in range(imax-1):
        aE[i, j] = -mx_m * Dy[j]
for j in range(jmax-1):
    for i in range(1, imax-1):
        aN[i, j] = -my_m * Dx[i]
for j in range(1, jmax-1):
    for i in range(1, imax-1):
        aP0[i, j] = rho * Dx[i] * Dy[j] / Dt
        aP[i, j] = aP0[i, j] + aE[i, j] + aE[i-1, j] + mx * Dy[j] + aN[i, j] + aN[i, j-1] + my * Dx[i]

#step 7 - function to compute face centre temperature
def temp_f_d(k, wf_p, wf_m, T1, T2, T3, T4):
    Tfd1 = wf_p[k, 0] * T3 + (wf_p[k, 1] - 1) * T2 + wf_p[k, 2] * T1
    Tfd2 = wf_m[k, 0] * T2 + (wf_m[k, 1] - 1) * T3 + wf_m[k, 2] * T4
    return Tfd1, Tfd2

#step 8 -Time marching for implicit method
unsteadiness_nd = 1
n = 0
Ntot = 0


while unsteadiness_nd >= epsilon_st:
    n += 1
    T[0, :] = T[1, :]
    T[:, 0] = T[:, 1]
    T_old = T.copy()
    epsilon = 1e-4
    N = 0
    Error = 1

    while Error >= epsilon:
        T[0, :] = T[1, :]
        T[:, 0] = T[:, 1]
        T_old_iter = T.copy()
        N += 1

        Ted_p = np.zeros((imax, jmax))
        Ted_m = np.zeros((imax, jmax))
        for j in range(1, jmax-1):
            for i in range(1, imax-1):
                if i == 1 or i == imax - 2:
                    Ted_p[i, j] = Ted_m[i, j] = 0
                else:
                    Ted_p[i, j], Ted_m[i, j] = temp_f_d(i, we_p, we_m,
                                                       T_old_iter[i-1, j], T_old_iter[i, j],
                                                       T_old_iter[i+1, j], T_old_iter[i+2, j])

        Tnd_p = np.zeros((imax, jmax))
        Tnd_m = np.zeros((imax, jmax))
        for j in range(jmax-1):
            for i in range(1, imax-1):
                if j == 0 or j == jmax - 2:
                    Tnd_p[i, j] = Tnd_m[i, j] = 0
                else:
                    Tnd_p[i, j], Tnd_m[i, j] = temp_f_d(j, wn_p, wn_m,
                                                       T_old_iter[i, j-1], T_old_iter[i, j],
                                                       T_old_iter[i, j+1], T_old_iter[i, j+2])

        for j in range(jmax-2, 0, -1):
            for i in range(imax-2, 0, -1):
                Qd_adv_old_iter = (mx_p * Ted_p[i, j] + mx_m * Ted_m[i, j]) * Dy[j] - \
                                  (mx_p * Ted_p[i-1, j] + mx_m * Ted_m[i-1, j]) * Dy[j]
                Qd_adv_old_iter += (my_p * Tnd_p[i, j] + my_m * Tnd_m[i, j]) * Dx[i] - \
                                   (my_p * Tnd_p[i, j-1] + my_m * Tnd_m[i, j-1]) * Dx[i]

                b = aP0[i, j] * T_old[i, j] - Qd_adv_old_iter

                T[i, j] = (aE[i, j] * T[i+1, j] + (aE[i-1, j] + mx * Dy[j]) * T[i-1, j] +
                           aN[i, j] * T[i, j+1] + (aN[i, j-1] + my * Dx[i]) * T[i, j-1] + b) / aP[i, j]

        Error = np.max(np.abs(T - T_old_iter))

    Ntot += N
    unsteadiness_nd = np.max(np.abs(T - T_old)) / Dt
    print(f"Time step no. {n:5d}, unsteadiness_nd = {unsteadiness_nd:8.4e}")

import matplotlib.pyplot as plt

# Create meshgrid for plotting
X, Y = np.meshgrid(xc, yc)

# Plotting the final temperature distribution
plt.figure(figsize=(8, 6))
cp = plt.contourf(X, Y, T.T, levels=50, cmap='jet')  # Transpose T to match X and Y orientation
plt.colorbar(cp, label='Temperature')
plt.title('2D Temperature Distribution')
plt.xlabel('X')
plt.ylabel('Y')
plt.grid(True)
plt.tight_layout()
plt.show()

# Centerline indices
i_center = imax // 2
j_center = jmax // 2

# Plot horizontal centerline temperature (T along x at mid y)
plt.figure(figsize=(7, 5))
plt.plot(xc, T[i_center, :], 'r-', marker='o', label='Horizontal Centerline (y = mid)')
plt.xlabel('x')
plt.ylabel('Temperature')
plt.title('Temperature along Horizontal Centerline')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Plot vertical centerline temperature (T along y at mid x)
plt.figure(figsize=(7, 5))
plt.plot(yc, T[:, j_center], 'b-', marker='s', label='Vertical Centerline (x = mid)')
plt.xlabel('y')
plt.ylabel('Temperature')
plt.title('Temperature along Vertical Centerline')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
