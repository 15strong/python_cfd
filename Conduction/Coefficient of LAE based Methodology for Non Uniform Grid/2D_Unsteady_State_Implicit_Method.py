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
zc=np.zeros_like(xc)
z=np.zeros_like(x)

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xc, yc, zc, color='blue', s=20, label='Cell centers')       # Blue dots
plt.show()