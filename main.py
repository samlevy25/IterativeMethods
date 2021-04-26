def hasDominantDiagonal(matrix):
    for i in range(len(matrix)):
        pivot, sum = abs(matrix[i][i]), 0
        for k in range(len(matrix)):
            sum = sum + abs(matrix[i][k])
        if sum - pivot > pivot:
            return False
    return True


def swap_row(matrix, i, j):
    temp = matrix[i]
    matrix[i] = matrix[j]
    matrix[j] = temp


def pivoting(matrix, b):
    for i in range(len(matrix)):
        for j in range(i + 1, len(matrix)):
            if abs(matrix[i][i]) < abs(matrix[j][i]):
                swap_row(matrix, i, j)
                swap_row(b, i, j)


def jacobiMethod(matrix, b):
    Xr, Yr, Zr, condition, count, epsilon = 0, 0, 0, 1, 0, 0.00001

    while condition > epsilon:
        count += 1
        Xr1 = (b[0] - matrix[0][1] * Yr - matrix[0][2] * Zr) / matrix[0][0]
        Yr1 = (b[1] - matrix[1][0] * Xr - matrix[1][2] * Zr) / matrix[1][1]
        Zr1 = (b[2] - matrix[2][1] * Yr - matrix[2][0] * Xr) / matrix[2][2]
        condition = abs(Xr1 - Xr)
        Xr, Yr, Zr = round(Xr1, 6), round(Yr1, 6), round(Zr1, 6)
        print("{0} Z = {1}, Y = {2}, X = {3}".format(count, Zr, Yr, Xr))


def gaussSeidelMethod(matrix, b):
    Xr, Yr, Zr, condition, count, epsilon = 0, 0, 0, 1, 0, 0.00001

    while condition > epsilon:
        count += 1
        Xr_1 = Xr
        Xr = (b[0] - matrix[0][1] * Yr - matrix[0][2] * Zr) / matrix[0][0]
        Yr = (b[1] - matrix[1][0] * Xr - matrix[1][2] * Zr) / matrix[1][1]
        Zr = (b[2] - matrix[2][1] * Yr - matrix[2][0] * Xr) / matrix[2][2]
        condition = abs(Xr - Xr_1)
        Xr, Yr, Zr = round(Xr, 6), round(Yr, 6), round(Zr, 6)
        print("{0} Z = {1}, Y = {2}, X = {3}".format(count, Zr, Yr, Xr))


def driver(matrix, b):
    pivoting(matrix, b)
    if not hasDominantDiagonal(matrix):
        print("\nThe Diagonal matrix is not Dominant.")
    print("\n********Jacobi Method*********")
    jacobiMethod(matrix, b)
    print("\n****Gauss Seidel Method****")
    gaussSeidelMethod(matrix, b)


A = [[4, 2, 0], [2, 10, 4], [0, 4, 5]]
B = [2, 6, 5]

driver(A, B)
