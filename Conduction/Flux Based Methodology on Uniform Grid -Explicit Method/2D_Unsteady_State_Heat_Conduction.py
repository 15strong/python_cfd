import numpy as np
import matplotlib.pyplot as plt

#Step 1 - User Input

L1 = 1.0
L2 = 1.0
rho = 7750.0
cp = 500.0
k = 16.2
T0 = 30
T_wb = 100
T_sb = 200
T_nb = 400
T_eb = 300
imax = 12
jmax = 12
epsilon_st=0.0001
Q_gen_vol = 0

#Step 2 - Geometrical Parameters
alpha = k/(rho*cp)
Dx= L1/(imax-1)
Dy= L2/(jmax-1)
Dt= 0.7*(1/(2*alpha*((1/Dx**2)+(1/Dy**2))))
DTc= T_nb - T_wb
Q_gen = Q_gen_vol*Dx*Dy
#Step 3 - IC and Boundary COnditions
T = np.zeros((imax,jmax))

#Internal cells
T[1:imax-1,1:jmax-1]=T0
#Boundary Cells
T[0,0:jmax]=T_wb
T[imax-1,0:jmax]=T_eb
T[0:imax,0]=T_sb
T[0:imax,jmax-1]=T_nb
#Time Marching for Explicit LAE's
unsteadines_nd=1
n=0
frames=[]
while unsteadines_nd>=epsilon_st:
    n+=1
    T_old = T.copy()
    # Reapply boundary conditions each iteration
    T[0,0:jmax] = T_wb
    T[imax-1,0:jmax] = T_eb
    T[0:imax,0] = T_sb
    T[0:imax,jmax-1] = T_nb
    #Step 4 - Computaion of conduction flux and Temp
    qx_old=np.zeros((imax,jmax))
    qy_old=np.zeros((imax,jmax))
    Q_cond_old = np.zeros((imax,jmax))
    for j in range(1,jmax-1):
        for i in range(0,imax-1):
            if i ==0:
                qx_old[i,j]= -k*(T_old[i+1,j]-T_old[i,j])/(Dx/2)
            elif i == imax-2:
                qx_old[i,j]= -k*(T_old[i+1,j]-T_old[i,j])/(Dx/2)
            else:
                qx_old[i,j]= -k*(T_old[i+1,j]-T_old[i,j])/Dx
    for j in range(0,jmax-1):
        for i in range(1,imax-1):
            if j ==1:
                qy_old[i,j]= -k*(T_old[i,j+1]-T_old[i,j])/(Dy/2)
            elif j==jmax-2:
                qy_old[i,j]= -k*(T_old[i,j+1]-T_old[i,j])/(Dy/2)
            else:
                qy_old[i,j]= -k*(T_old[i,j+1]-T_old[i,j])/Dy

    for i in range(1,imax-1):
        for j in range(1,jmax-1):
            Q_cond_old[i,j] = ((qx_old[i-1,j]-qx_old[i,j])*Dy)+((qy_old[i,j-1]-qy_old[i,j])*Dx)
            T[i,j] = T_old[i,j]+((Dt/(rho*cp*Dx*Dy))*(Q_cond_old[i,j]+Q_gen))
    frames.append(T.copy())
    unsteadines = np.max(np.abs(T-T_old))/Dt
    unsteadines_nd = unsteadines*L1*L2/(alpha*DTc)


print(T)


import matplotlib.animation as animation

# Non-uniform grid again
x = np.zeros(imax)
x[0] = 0
x[1] = Dx / 2
for i in range(2, imax-1):
    x[i] = x[i-1] + Dx
x[imax-1] = L1

y = np.zeros(jmax)
y[0] = 0
y[1] = Dy / 2
for j in range(2, jmax-1):
    y[j] = y[j-1] + Dy
y[jmax-1] = L2

X, Y = np.meshgrid(x, y)

fig, ax = plt.subplots(figsize=(8, 6))
levels = np.linspace(np.min(T), np.max(T), 20)
contour = ax.contourf(X, Y, frames[0].T, levels=levels, cmap='rainbow')
cbar = plt.colorbar(contour)
cbar.set_label('Temperature')

def update(frame):
    ax.clear()
    contour = ax.contourf(X, Y, frame.T, levels=levels, cmap='rainbow')
    ax.set_title('Temperature Iteration')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_aspect('equal', adjustable='box')
    return contour

ani = animation.FuncAnimation(fig, update, frames=frames, interval=100)

#To display in notebook (if using Jupyter)
#from IPython.display import HTML
#HTML(ani.to_jshtml())

#To save as MP4
#ani.save('heat_conduction_animation.mp4', writer='ffmpeg', fps=10)

# Or to save as GIF
ani.save('heat_conduction_animation.gif', writer='pillow', fps=10)

plt.close()
