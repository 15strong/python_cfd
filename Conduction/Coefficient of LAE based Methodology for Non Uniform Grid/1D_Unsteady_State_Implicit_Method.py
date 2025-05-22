import numpy as np
import matplotlib.pyplot as plt

#grid generation
L=1.0
imax=12
beta=1.2
xi = np.linspace(0, 1, imax)  # 12 points between 0 and 1

#Clustering
beta_p1 = beta + 1
beta_m1 = beta - 1
beta_p1_div_m1 = (beta_p1 / beta_m1) ** (2*xi-1)
num = ( beta_p1 * beta_p1_div_m1 ) - beta_m1
den = 1 + beta_p1_div_m1
xc=np.zeros(imax)
Dx=np.zeros(imax)
dx=np.zeros(imax-1)
x=L*(num/den)
xc[1:imax-1]=0.5*(x[1:imax-1]+x[:imax-2])
xc[0]=x[0]
xc[imax-1]=x[imax-2]
Dx[1:imax-1]=x[1:imax-1]-x[:imax-2]
dx=xc[1:]-xc[:imax-1]


#Coeff. of LAE based methodology non uniform cartesian grid of 1D

#step 1
rho=7750
cp=500
k=16.2
T0=30
T_wb=0
h=1000
T_inf=100
epsilon_st=0.0001
Q_vol_gen=0
Q_gen=Q_vol_gen*Dx
# Dt=1
alpha=k/(rho*cp)
Dt = 0.25 * np.min(Dx[1:imax-1]**2) / alpha  # Stability condition


#step 2
aE=np.zeros(imax)
aP=np.zeros(imax)
aW=np.zeros(imax)
aP0=np.zeros(imax)
b=np.zeros(imax)

for i in range(1, imax-1):
    aE[i] = k / (xc[i+1] - xc[i]) if i < imax-1 else 0
    aW[i] = k / (xc[i] - xc[i-1]) if i > 0 else 0
    aP0[i] = rho * cp * Dx[i] / Dt
    aP[i] = aP0[i] + aE[i] + aW[i]

#step 3  IC and Boundary Conditions
T=np.zeros(imax)
T[1:imax-1]=T0
T[0]=T_wb
unsteadiness_nd=1
n=0

DTc=T_inf-T_wb


while unsteadiness_nd >= epsilon_st:
    n += 1
    T_old = T.copy()
    
    # Apply boundary conditions
    T[0] = T_wb  # Fixed temperature at left boundary
    # Convective boundary at right
    T[-1] = (k * T[-2] + h * Dx[-2] * T_inf) / (k + h * Dx[-2])
    
    # Update source terms
    for i in range(1, imax-1):
        b[i] = aP0[i] * T_old[i] + Q_gen[i]
    
    # Iterative solution
    epsilon = 0.0001
    error = 1
    N = 0
    
    while error >= epsilon:
        T_old_iter = T.copy()
        
        # Gauss-Seidel iteration
        for i in range(1, imax-1):
            T[i] = (aE[i] * T[i+1] + aW[i] * T[i-1] + b[i]) / aP[i]
        
        # Update boundary condition in iteration
        T[-1] = (k * T[-2] + h * Dx[-2] * T_inf) / (k + h * Dx[-2])
        
        error = np.max(np.abs(T - T_old_iter))
        N += 1
    
    # Check convergence
    unsteadiness = np.max(np.abs(T - T_old))
    unsteadiness_nd = (unsteadiness * L**2) / (alpha * DTc)
print(T,n,N,unsteadiness_nd)

# Plot the results
plt.figure(figsize=(8, 6))
plt.plot(xc, T, 'bo-', label='Temperature')
plt.xlabel('Position (m)')
plt.ylabel('Temperature (°C)')
plt.title('1D Transient Heat Conduction')
plt.grid(True)
plt.legend()
plt.show()
