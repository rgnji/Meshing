import numpy as np
import matplotlib.pyplot as plt

from lib.laplace import solve_grid_laplace_sor
from lib.cylindermesh import initial_o, initial_square
from lib.plot3dout import bin_fort12

#========== parameters ==========
# cylinder radius
ro = 5
# radius of inner square
ri = 10
# axial length of cylinder
z_len = 6

# number of division in arc, radius and axial direction
arc_number = 10
radius_number = 5
axis_number = 10


"""
#========== create 2D meshing ==========
plt.figure(figsize=(6, 6))

# generate o-grid
for qua in range(1,5):
    x,y = initial_o(qua, ro, ri, arc_number, radius_number)
    x_smooth, y_smooth, iterations = solve_grid_laplace_sor(x, y)
    print(f"SOR of quadrant {qua} converged in {iterations} iterations.")

    # plot mesh
    for i in range(x_smooth.shape[0]):
        plt.plot(x_smooth[i, :], y_smooth[i, :], color='blue', linewidth=0.8)
    for j in range(x_smooth.shape[1]):
        plt.plot(x_smooth[:, j], y_smooth[:, j], color='blue', linewidth=0.8)

# generate inner square
x, y = initial_square(ri, arc_number)
x_smooth, y_smooth, iterations = solve_grid_laplace_sor(x, y)
print(f"SOR of inner square converged in {iterations} iterations.")
# plot mesh
for i in range(x.shape[0]):
    plt.plot(x[i, :], y[i, :], color='blue', linewidth=0.8)
for j in range(x.shape[1]):
    plt.plot(x[:, j], y[:, j], color='blue', linewidth=0.8)


plt.gca().set_aspect('equal')
plt.title("Optimized Grid (Elliptic Smoothing)")
plt.xlabel("x")
plt.ylabel("y")
plt.grid(False)
plt.show()
"""


#========== main program ==========
z = np.linspace(0, z_len, axis_number)
XT = []
YT = []
ZT = []

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# four quadrant
for qua in range(1,5): 
    # create 2D mesh
    x,y = initial_o(qua, ro, ri, arc_number, radius_number)
    # solve laplace equation
    x_smooth, y_smooth, iterations = solve_grid_laplace_sor(x, y)
    print(f"SOR of quadrant {qua} converged in {iterations} iterations.")
    # create 3D mesh
    X = np.full( (axis_number, x_smooth.shape[0], x_smooth.shape[1]), x_smooth)
    Y = np.full( (axis_number, y_smooth.shape[0], y_smooth.shape[1]), y_smooth)
    jj, Z, ii = np.meshgrid(np.arange(x_smooth.shape[0]),   # j, k, i
                          z,
                          np.arange(x_smooth.shape[1]))
    
    XT.append(X)
    YT.append(Y)
    ZT.append(Z)
    

# inner square
x, y = initial_square(ri, arc_number)  

x_smooth, y_smooth, iterations = solve_grid_laplace_sor(x, y)
print(f"SOR of inner square converged in {iterations} iterations.")

X = np.full( (axis_number, x_smooth.shape[0], x_smooth.shape[1]), x_smooth)
Y = np.full( (axis_number, y_smooth.shape[0], y_smooth.shape[1]), y_smooth)
jj, Z, ii = np.meshgrid(np.arange(x_smooth.shape[0]),   # j, k, i
                      z,
                      np.arange(x_smooth.shape[1]))

XT.append(X)
YT.append(Y)
ZT.append(Z)


# establish fort.12
txt = bin_fort12(XT, YT, ZT)
print(txt)


# plot
for q in range(5):
    ax.scatter(XT[q], YT[q], ZT[q])

plt.show()

