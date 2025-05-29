import numpy as np
import matplotlib.pyplot as plt

#Step 1 -user imput
L1=1
L2=1
rho=1000
cp=4180
imax=12
jmax=12
T0=0.5
T_wb=1
T_sb=0
u=1
v=1
epsilon_st=0.0001
print("SELECT THE ADVECTION SCHEME(1/2/3)")
print("1.FOU\n")
print("2.SOU\n")
print("3.QUICK\n")
scheme = int(input("Enter the Advection Scheme:"))

def weights(scheme, k):
    w = [0, 0, 0]  # initialize list with 3 elements

    if scheme == 1:
        w[0] = 0
        w[1] = 1
        w[2] = 0

    elif scheme == 2:
        if k == 2:
            w[0] = 0
            w[1] = 2
            w[2] = -1
        else:
            w[0] = 0
            w[1] = 3/2
            w[2] = -1/2

    elif scheme == 3:
        if k == 2:
            w[0] = 1/3
            w[1] = 1
            w[2] = -1/3
        else:
            w[0] = 3/8
            w[1] = 6/8
            w[2] = -1/8

    return w


Dx=L1/(imax-2)
Dy=L2/(jmax-2)
Dt=1/((np.abs(u)/Dx)+(np.abs(v)/Dy))

if scheme == 2:
    Dt *= 0.2
elif scheme == 3:
    Dt *= 0.1

T=np.zeros((imax,jmax))
hx_old = np.zeros_like(T)
hy_old=np.zeros_like(T)
Q_adv_old=np.zeros_like(T)

#step 3 IC and BCS
T[1:imax-1,1:jmax-1]=T0
T[0,:]=T_wb
T[:,0]=T_sb

#step 4 calculate mass flux
mx = rho*u
my=rho*v
mx_p=max(mx,0)
mx_m=min(mx,0)
my_p=max(my,0)
my_m=min(my,0)


def temp_f(w, T1, T2, T3, T4):
    Tf_p = w[0]*T3 + w[1]*T2 + w[2]*T1
    Tf_m = w[0]*T2 + w[1]*T3 + w[2]*T4
    return Tf_p, Tf_m

#step 5 time marching equation
unsteadiness_nd=1
n=0

while unsteadiness_nd>=epsilon_st:
    n+=1
    T_old=T.copy()
    T[imax-1,:]=T[imax-2,:]
    T[:,jmax-1]=T[:,jmax-2]

    #step 6 computation of enthalp-flux and temp

    for j in range (1,jmax-1):
        for i in range (0,imax-1):
            if i==0 or i==imax-2:
                Te_p=T_old[i,j]
                Te_m=T_old[i+1,j]
            else:
                 w = weights(scheme, i)
                 Te_p, Te_m = temp_f(w, T_old[i-1, j], T_old[i, j], T_old[i+1, j], T_old[i+2, j])
            hx_old[i, j] = cp * (mx_p*Te_p + mx_m*Te_m)
    for j in range (0,jmax-1):
        for i in range (1,imax-1):
            if j==0 or j==jmax-2:
                Tn_p=T_old[i,j]
                Tn_m=T_old[i,j+1]
            else:
                w = weights(scheme, i)
                Tn_p, Tn_m = temp_f(w, T_old[i, j-1], T_old[i, j], T_old[i, j+1], T_old[i, j+2])
            hy_old[i, j] = cp * (my_p*Tn_p + my_m*Tn_m)
    for j in range (1,jmax-1):
        for i in range (1,imax-1):
            Q_adv_old[i,j]=((hx_old[i,j]-hx_old[i-1,j])*Dy)+((hy_old[i,j]-hy_old[i,j-1])*Dx)
            T[i,j]=T_old[i,j]-(Dt/(rho*cp*Dx*Dy))*Q_adv_old[i,j]
    unsteadiness_nd=np.max(np.abs(T-T_old))/Dt

print(f"Time step {n:5d}, Unsteadiness = {unsteadiness_nd:.4e}") 
    
#plotting
i_center = (imax-2) // 2
theta = (T[i_center, :] - T_sb) / (T_wb - T_sb)
y_vals = np.linspace(0, L2, jmax)

plt.figure(figsize=(6, 5))
plt.plot(theta, y_vals, marker='o', color='blue')
plt.title("Non-Dimensional Temperature along Vertical Centerline")
plt.xlabel("θ = (T - T_sb)/(T_wb - T_sb)")
plt.ylabel("y (Vertical position)")
plt.grid(True)
plt.ylim(0, L2)
plt.tight_layout()
plt.show()