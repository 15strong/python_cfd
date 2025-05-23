import numpy as np
import matplotlib.pyplot as plt

#grid generation

L1=1
L2=1
imax=12
jmax=12
beta=1.2
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
dx=np.zeros(imax)
dy=np.zeros(jmax)

xc[1:imax-1]=0.5*(x[1:imax-1]+x[0:imax-2])
xc[0]=x[0]
xc[imax-1]=x[imax-2]
yc[1:jmax-1]=0.5*(y[1:jmax-1]+y[0:jmax-2])
yc[0]=y[0]
yc[jmax-1]=y[jmax-2]


XC,YC=np.meshgrid(xc,yc)
X,Y=np.meshgrid(x,y)

Dx[1:imax -1]=x[1:imax -1]-x[0:imax -2] 
Dy=Dx
dx[:imax -1]=xc[1: imax]-xc[: imax -1]
dy=dx


zc=np.zeros_like(xc)
z=np.zeros_like(x)

#Solution
#step 1
rho = 7750
cp=500
k=16.2
T0=30
T_wb=100
T_inf=30                        #at north boundary
h = 100
q_w=10000
Q_vol_gen=0
Q_gen=Q_vol_gen*Dx*Dy
epsilon_st=0.0001
Dt=1000


#step 2
aE=np.zeros((imax,jmax))
aW=np.zeros((imax,jmax))
aN=np.zeros((imax,jmax))
aS=np.zeros((imax,jmax))
aP0=np.zeros((imax,jmax))
aP=np.zeros((imax,jmax))
b=np.zeros((imax,jmax))


for j in range (1,jmax-1):
    for i in range (0,imax-1):
        aE[i,j]=k*Dy[j]/dx[i]
for j in range (0,jmax-1):
    for i in range (1,imax-1):
        aN[i,j]=k*Dx[i]/dy[j]
for j in range (1,jmax-1):
    for i in range (1,imax-1):
        aP0[i,j]=rho*cp*Dx[i]*Dy[j]/Dt
        aP[i,j]=aP0[i,j]+aE[i,j]+aE[i-1,j]+aN[i,j]+aN[i,j-1]

#step 3
T=np.zeros((imax,jmax))

for i in range (1,imax-1):
    for j in range (1,jmax-1):
        T[i,j]=T0
T[1,1:]=T_wb
unsteadiness_nd=1
n=0
alpha=k/(rho*cp)
DTc=T_wb-T_inf


while unsteadiness_nd>=epsilon_st:
    n+=1
    T_old=T.copy()
    for j in range (0,jmax):
        T[imax-1,j]=T[imax-2,j]+(q_w*dx[imax-2]/k)
    for i in range (0,imax):
        T[i,0]=T[i,1]
        T[i,jmax-1]=((k*T[i,jmax-2])+(h*dy[jmax-1]*T_inf))/(k+h*dy[jmax-1])
    
    for j in range (1,jmax-1):
        for i in range(1,imax-1):
            b[i,j]=aP0[i,j]*T_old[i,j]+Q_vol_gen*Dx[i]*Dy[j]

    epsilon=0.0001
    N=0
    error=1
    while error>=epsilon:
        T_old_iter=T.copy()
        
        for j in range (0,jmax):
            T[i,jmax-1]=((k*T[i,jmax-2])+(h*dy[jmax-2]*T_inf))/(k+h*dy[jmax-2])
        
        for j in range (1,jmax-1):
            for i in range (1,imax-1):
                T[i,j]=aE[i,j]*T[i+1,j]+aE[i-1,j]*T[i-1,j]+aN[i,j]*T[i,j+1]+aN[i,j-1]*T[i,j-1]+b[i,j]
                T[i,j]=T[i,j]/aP[i,j]
        N+=1
        error=np.max(np.abs(T-T_old_iter))
    unsteadiness=np.max(np.abs(T-T_old))/Dt
    unsteadiness_nd=unsteadiness*L1*L1/(alpha*DTc)

print(T)

# Return T and grid for plotting
T_plot = T.T  # Transpose for correct orientation
X_plot = xc.T
Y_plot = yc.T

T_plot, X_plot, Y_plot

plt.figure(figsize=(8, 6))
cp = plt.contourf(X_plot, Y_plot, T_plot, levels=50, cmap='plasma')
plt.colorbar(cp, label='Temperature (°C)')
plt.title('Temperature Distribution in the Plate')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.grid(True)
plt.axis('equal')
plt.show()