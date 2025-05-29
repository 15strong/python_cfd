# import numpy as np
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

# # Parameters
# L1 = L2 = 1
# imax = jmax = 12
# Beta = 1.2

# # Step 1-3: Generate stretched grid
# xi = np.linspace(0, 1, imax - 1)
# Beta_p1 = Beta + 1
# Beta_m1 = Beta - 1
# Beta_ratio = (Beta_p1 / Beta_m1) ** (2 * xi - 1)
# num = (Beta_p1 * Beta_ratio) - Beta_m1
# den = 2 * (1 + Beta_ratio)
# x = L1 * num / den
# y = x.copy()

# # Step 4: Compute cell centers
# xc = np.zeros(imax)
# xc[1:imax-1] = 0.5 * (x[1:] + x[:-1])
# xc[0] = x[0]
# xc[-1] = x[-1]
# yc = xc.copy()

# # Meshgrids
# X, Y = np.meshgrid(x, y)
# Xc, Yc = np.meshgrid(xc, yc)
# Z = np.zeros_like(X)
# Zc = np.zeros_like(Xc)

# # Plotting
# fig = plt.figure(figsize=(8, 6))
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(X, Y, Z, color='lightgray', alpha=0.7, edgecolor='k')  # Flat surface
# ax.scatter(Xc, Yc, Zc, color='blue', s=20, label='Cell centers')       # Blue dots

# ax.set_xlabel("X")
# ax.set_ylabel("Y")
# ax.set_zlabel("Z")
# ax.set_title("Stretched Mesh and Cell Centers")
# ax.legend()
# plt.tight_layout()
# plt.show()






# import numpy as np
# import matplotlib.pyplot as plt

# # Parameters
# L1 = L2 = 1
# imax = jmax = 12
# Beta = 1.2

# # Step 1-3: Generate stretched grid
# xi = np.linspace(0, 1, imax - 1)
# Beta_p1 = Beta + 1
# Beta_m1 = Beta - 1
# Beta_ratio = (Beta_p1 / Beta_m1) ** (2 * xi - 1)
# num = (Beta_p1 * Beta_ratio) - Beta_m1
# den = 2 * (1 + Beta_ratio)
# x = L1 * num / den
# y = x.copy()

# # Step 4: Compute cell centers
# xc = np.zeros(imax)
# xc[1:imax-1] = 0.5 * (x[1:] + x[:-1])
# xc[0] = x[0]
# xc[-1] = x[-1]
# yc = xc.copy()

# # Meshgrids
# X, Y = np.meshgrid(x, y)
# Xc, Yc = np.meshgrid(xc, yc)

# # 2D Plot
# plt.figure(figsize=(8, 8))

# # Plot vertical lines
# for i in range(len(x)):
#     plt.plot([x[i]] * len(y), y, color='gray', linewidth=0.8)

# # Plot horizontal lines
# for j in range(len(y)):
#     plt.plot(x, [y[j]] * len(x), color='gray', linewidth=0.8)

# # Plot cell centers
# plt.scatter(Xc, Yc, color='blue', label='Cell Centers')

# plt.xlabel("X")
# plt.ylabel("Y")
# plt.title("2D View of Stretched Mesh and Cell Centers")
# plt.axis("equal")
# plt.legend()
# plt.grid(True)
# plt.show()





# import numpy as np

# T0=30
# imax=12
# T= np.zeros(imax)
# for i in range [1,3]:
#     T[i]=T0
#     print(T)

# import numpy as np
# import matplotlib.pyplot as plt

# # Grid generation and clustering
# L1 = L2 = 1
# imax = jmax = 12
# beta = 1.2
# xi = np.linspace(0, 1, imax-1)

# beta_p1 = beta + 1
# beta_m1 = beta - 1
# beta_ratio = (beta_p1 / beta_m1) ** (2 * xi - 1)
# num = beta_p1 * beta_ratio - beta_m1
# den = 1 + beta_ratio
# x = L1 * (num / den)
# y = L2 * (num / den)

# # Calculate cell centers and grid spacing
# xc = np.zeros(imax)
# yc = np.zeros(jmax)
# xc[1:imax-1] = 0.5 * (x[1:imax-1] + x[0:imax-2])
# xc[0] = x[0]
# xc[-1] = x[-1]
# yc[1:jmax-1] = 0.5 * (y[1:jmax-1] + y[0:jmax-2])
# yc[0] = y[0]
# yc[-1] = y[-1]

# Dx = np.zeros(imax)
# Dy = np.zeros(jmax)
# Dx[1:-1] = x[1:] - x[:-1]
# Dx[0] = Dx[1]
# Dx[-1] = Dx[-2]
# Dy[1:-1] = y[1:] - y[:-1]
# Dy[0] = Dy[1]
# Dy[-1] = Dy[-2]

# # Constants
# rho = 7750
# cp = 500
# k = 16.2
# T0 = 30
# T_wb = 100
# T_inf = 30
# h = 100
# q_w = 10000
# Dt = 1000
# epsilon_st = 1e-4
# alpha = k / (rho * cp)
# DTc = T_wb - T_inf

# # Initialize temperature
# T = np.full((imax, jmax), T0)
# T[:, 1] = T_wb

# # Coefficients
# aP0 = np.zeros((imax, jmax))
# aP = np.zeros((imax, jmax))
# aE = np.zeros((imax, jmax))
# aW = np.zeros((imax, jmax))
# aN = np.zeros((imax, jmax))
# aS = np.zeros((imax, jmax))

# for j in range(1, jmax-1):
#     for i in range(imax-1):
#         aE[i, j] = k * Dy[j] / Dx[i]
# for j in range(jmax-1):
#     for i in range(1, imax-1):
#         aN[i, j] = k * Dx[i] / Dy[j]
# for j in range(1, jmax-1):
#     for i in range(1, imax-1):
#         aP0[i, j] = rho * cp * Dx[i] * Dy[j] / Dt
#         aP[i, j] = aP0[i, j] + aE[i, j] + aE[i-1, j] + aN[i, j] + aN[i, j-1]

# # Time iteration
# unsteadiness_nd = 1
# while unsteadiness_nd >= epsilon_st:
#     T_old = T.copy()
    
#     # Apply boundary conditions
#     T[-1, :] = T[-2, :] + (q_w * Dx[-2] / k)
#     T[:, 0] = T[:, 1]
#     T[:, -1] = (k * T[:, -2] + h * Dy[-1] * T_inf) / (k + h * Dy[-1])
    
#     b = aP0 * T_old
    
#     error = 1
#     while error > 1e-4:
#         T_iter = T.copy()
#         for j in range(1, jmax-1):
#             for i in range(1, imax-1):
#                 T[i, j] = (
#                     aE[i, j] * T[i+1, j] +
#                     aE[i-1, j] * T[i-1, j] +
#                     aN[i, j] * T[i, j+1] +
#                     aN[i, j-1] * T[i, j-1] +
#                     b[i, j]
#                 ) / aP[i, j]
#         error = np.max(np.abs(T - T_iter))

#     unsteadiness = np.max(np.abs(T - T_old)) / Dt
#     unsteadiness_nd = unsteadiness * L1**2 / (alpha * DTc)

# # Prepare meshgrid for visualization
# Xc, Yc = np.meshgrid(xc, yc)

# # Return T and grid for plotting
# T_plot = T.T  # Transpose for correct orientation
# X_plot = Xc.T
# Y_plot = Yc.T

# T_plot, X_plot, Y_plot
# import matplotlib.pyplot as plt

# plt.figure(figsize=(8, 6))
# cp = plt.contourf(X_plot, Y_plot, T_plot, levels=50, cmap='plasma')
# plt.colorbar(cp, label='Temperature (°C)')
# plt.title('Temperature Distribution in the Plate')
# plt.xlabel('x (m)')
# plt.ylabel('y (m)')
# plt.grid(True)
# plt.axis('equal')
# plt.show()


import numpy as np

# --- Step 1: Parameters ---
rho, cp = 1000.0, 4180.0
L1 = L2 = 1.0
imax = jmax = 12
T0, T_wb, T_sb = 0.5, 1.0, 0.0
u = v = 1.0
epsilon_st = 1e-4

# --- Step 2: Grid & Time Step ---
Dx = L1 / (imax - 2)
Dy = L2 / (jmax - 2)
Dt = 1.0 / (abs(u)/Dx + abs(v)/Dy)

# Choose scheme: 1=FOU, 2=SOU, 3=QUICK
scheme = int(input("Enter advection scheme (1=FOU, 2=SOU, 3=QUICK): "))
if scheme == 2:
    Dt *= 2/3
elif scheme == 3:
    Dt *= 4/9

# --- Step 3: Initialize Temp Field ---
T = np.full((imax+1, jmax+1), T0)
T[0, :] = T_wb
T[:, 0] = T_sb

# --- Step 4: Mass fluxes ---
mx, my = rho*u, rho*v
mx_p, mx_m = max(mx, 0), min(mx, 0)
my_p, my_m = max(my, 0), min(my, 0)

# --- Step 5: Weight function ---
def weights(scheme, k):
    if scheme == 1:
        return [0, 1, 0]
    elif scheme == 2:
        return [0, 2, -1] if k == 2 else [0, 1.5, -0.5]
    elif scheme == 3:
        return [1/3, 1, -1/3] if k == 2 else [3/8, 6/8, -1/8]

# --- Step 6: Face temperature interpolation ---
def temp_f(w, T1, T2, T3, T4):
    Tf_p = w[0]*T3 + w[1]*T2 + w[2]*T1
    Tf_m = w[0]*T2 + w[1]*T3 + w[2]*T4
    return Tf_p, Tf_m

# --- Step 7: Time loop ---
unsteadiness = 1.0
n = 0
while unsteadiness >= epsilon_st:
    n += 1
    T_old = T.copy()

    # Apply outlet BCs
    T[imax, :] = T[imax-1, :]
    T[:, jmax] = T[:, jmax-1]

    # X-direction enthalpy flux
    hx = np.zeros_like(T)
    for j in range(1, jmax):
        for i in range(imax):
            if i == 0 or i == imax - 1:
                Te_p = T_old[i, j]
                Te_m = T_old[i+1, j]
            else:
                w = weights(scheme, i)
                Te_p, Te_m = temp_f(w, T_old[i-1, j], T_old[i, j], T_old[i+1, j], T_old[i+2, j])
            hx[i, j] = cp * (mx_p*Te_p + mx_m*Te_m)

    # Y-direction enthalpy flux
    hy = np.zeros_like(T)
    for j in range(jmax):
        for i in range(1, imax):
            if j == 0 or j == jmax - 1:
                Tn_p = T_old[i, j]
                Tn_m = T_old[i, j+1]
            else:
                w = weights(scheme, j)
                Tn_p, Tn_m = temp_f(w, T_old[i, j-1], T_old[i, j], T_old[i, j+1], T_old[i, j+2])
            hy[i, j] = cp * (my_p*Tn_p + my_m*Tn_m)

    # Update Temperature
    for j in range(1, jmax):
        for i in range(1, imax):
            Q_adv = (hx[i, j] - hx[i-1, j])*Dy + (hy[i, j] - hy[i, j-1])*Dx
            T[i, j] = T_old[i, j] - Dt / (rho*cp*Dx*Dy) * Q_adv

    # Convergence check
    unsteadiness = np.max(np.abs(T - T_old)) / Dt
    print(f"Time step {n:5d}, Unsteadiness = {unsteadiness:.4e}")

import matplotlib.pyplot as plt

# --- 1. 2D Contour Plot of Final Temperature Field ---
x = np.linspace(0, L1, imax+1)
y = np.linspace(0, L2, jmax+1)
X, Y = np.meshgrid(x, y, indexing='ij')

plt.figure(figsize=(8, 6))
cp_plot = plt.contourf(X, Y, T.T, 20, cmap='plasma')
plt.colorbar(cp_plot, label='Temperature')
plt.title("Steady-State Temperature Distribution (Contour)")
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
plt.tight_layout()
plt.show()

# --- 2. Vertical Centerline Non-Dimensional Profile ---
i_center = imax // 2
theta = (T[i_center, :] - T_sb) / (T_wb - T_sb)
y_vals = np.linspace(0, L2, jmax+1)

plt.figure(figsize=(6, 5))
plt.plot(theta, y_vals, marker='o', color='blue')
plt.title("Non-Dimensional Temperature along Vertical Centerline")
plt.xlabel("θ = (T - T_sb)/(T_wb - T_sb)")
plt.ylabel("y (Vertical position)")
plt.grid(True)
plt.ylim(0, L2)
plt.tight_layout()
plt.show()
