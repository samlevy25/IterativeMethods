def precision(num):
    epsilon = 1.0
    while (1.0 + 0.5 * epsilon) != 1:
        epsilon = 0.5 * epsilon
    if abs(num) <= epsilon:
        return 0
    else:
        return num


def swap_row(A, i, j):
    temp = A[i]
    A[i] = A[j]
    A[j] = temp


def pivoting(matrix, b):
    for i in range(len(matrix)):
        for j in range(i + 1, len(matrix)):
            if matrix[i][i] < matrix[j][i]:
                swap_row(matrix, i, j)
                swap_row(b, i, j)


def matrixMultiply(matrix_1, matrix_2):
    result = [[0] * len(matrix_2) for i in range(len(matrix_1))]  # Initializing the result matrix (with zeros)
    for i in range(len(matrix_1)):
        for j in range(len(matrix_2)):
            for k in range(len(matrix_2)):
                result[i][j] += precision(matrix_1[i][k] * matrix_2[k][j])
    return result


def createMatrix(matrix, size, val):
    if val == 0:
        matrix = [[0 for i in range(size)] for j in range(size)]
    elif val == 1:
        matrix = [[int(i == j) for i in range(size)] for j in range(size)]
    return matrix


def calcDet(matrix):
    size = len(matrix)
    det = 0
    if size == 2:
        det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
        return det
    minor = createMatrix(matrix, size - 1, 0)
    for k in range(len(matrix)):
        i, j = 0, 0
        while i < size:
            if i != k:
                minor[j] = matrix[i][1:]
                j += 1
            i += 1
        det += matrix[k][0] * ((-1) ** k) * calcDet(minor)
    return det


def invertMatrix(matrix):
    determinant = calcDet(matrix)
    if len(matrix) == 2:
        return [[matrix[1][1] / determinant, -1 * matrix[0][1] / determinant],
                [-1 * matrix[1][0] / determinant], matrix[0][0] / determinant]

    inverse = []
    for i in range(len(matrix)):
        inverseRow = []
        for j in range(len(matrix)):
            minor = [row[:j] + row[j + 1:] for row in (matrix[:i] + matrix[i + 1:])]
            inverseRow.append(((-1) ** (i + j)) * calcDet(minor))
        inverse.append(inverseRow)
    inverse = list(map(list, zip(*inverse)))
    for i in range(len(inverse)):
        for j in range(len(inverse)):
            inverse[i][j] = inverse[i][j] / determinant
    return inverse


def DLU(matrix):  # create D, L and U matrix
    D = []
    L = []
    U = []
    createMatrix(D, len(matrix), 0)
    createMatrix(L, len(matrix), 0)
    createMatrix(U, len(matrix), 0)
    for i in range(len(matrix)):
        D[i][i] = matrix[i][i]
        for n in range(i):
            L[i][n] = matrix[i][n]
        for m in range(i, len(matrix)):
            U[i][m] = matrix[i][m]
    return D, L, U


A = [[4, 2, 0], [2, 10, 4], [0, 4, 5]]
B = [2, 6, 5]


def jacobiMethod(matrix, b):
    Xr, Yr, Zr, condition, count = 0, 0, 0, 1, 0
    while condition > 0.001:
        count += 1
        Xr1 = (b[0] - matrix[0][1] * Yr - matrix[0][2] * Zr)/matrix[0][0]
        Yr1 = (b[1] - matrix[1][0] * Xr - matrix[1][2] * Zr)/matrix[1][1]
        Zr1 = (b[2] - matrix[2][1] * Yr - matrix[2][0] * Xr)/matrix[2][2]
        condition = abs(Xr1 - Xr)
        Xr, Yr, Zr = round(Xr1, 6), round(Yr1, 6), round(Zr1, 6)
        print("{0} Z = {1}, Y = {2}, X = {3}".format(count, Zr, Yr, Xr))


jacobiMethod(A, B)
