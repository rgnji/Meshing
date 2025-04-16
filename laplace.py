#========== solve laplace equation ==========
def solve_grid_laplace_sor(x, y, omega=1.5, tol=1e-6, max_iter=10000):
    x_new = x.copy()
    y_new = y.copy()

    ny, nx = x.shape
    count = 0
    for _ in range(max_iter):
        max_residual = 0.0
        for i in range(1, ny - 1):
            for j in range(1, nx - 1):
                # SOR update for x
                x_old = x_new[i, j]
                x_sor = 0.25 * (x_new[i+1, j] + x_new[i-1, j] + x_new[i, j+1] + x_new[i, j-1])
                x_new[i, j] = (1 - omega) * x_old + omega * x_sor
                res_x = abs(x_new[i, j] - x_old)

                # SOR update for y
                y_old = y_new[i, j]
                y_sor = 0.25 * (y_new[i+1, j] + y_new[i-1, j] + y_new[i, j+1] + y_new[i, j-1])
                y_new[i, j] = (1 - omega) * y_old + omega * y_sor
                res_y = abs(y_new[i, j] - y_old)

                # 更新最大殘差
                max_residual = max(max_residual, res_x, res_y)

        count += 1
        if max_residual < tol:
            break

    return x_new, y_new, count