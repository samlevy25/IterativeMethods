def precision(num):
    epsilon = 1.0
    while (1.0 + 0.5 * epsilon) != 1:
        epsilon = 0.5 * epsilon
    if abs(num) <= epsilon:
        return 0
    else:
        return num


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


def jacobiMethod(D, L, U, b, x):
    x = matrixMultiply(-invertMatrix(D), L + U) * x + invertMatrix(D) * b
