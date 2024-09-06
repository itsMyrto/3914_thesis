from math import ceil

def sieve_of_eratosthenes(B: int) -> list[int]:
    primes = [True] * (B + 1)
    primes[0] = primes[1] = False
    root = ceil(pow(B, 0.5))

    for i in range(2, root):
        if primes[i]:
            for j in range(i * i, B + 1, i):
                primes[j] = False

    return [i for i in range(0, len(primes)) if primes[i]]


def legendre_symbol(n: int, p: int) -> int:
    return pow(n, (p - 1) // 2, p)


def tonelli_shanks(n: int, p: int) -> tuple[int, int]:
    S = 0
    Q = p - 1
    while Q % 2 == 0:
        Q //= 2
        S += 1

    if S == 1:
        R = pow(n, (p+1)//4, p)
        return R, p-R

    z = 2
    while legendre_symbol(z, p) != p-1:
        z += 1

    M = S
    c = pow(z, Q, p)
    t = pow(n, Q, p)
    R = pow(n, (Q+1)//2, p)

    while True:
        if t % p == 1:
            break

        i = 1
        while pow(t, pow(2, i), p) != 1:
            i += 1

        b = pow(c, 2 ** (M-i-1), p)
        M = i
        c = (b*b) % p
        t = (t*b*b) % p
        R = (R*b) % p

    return R, p-R


def gcd(x: int, y: int) -> int:
    if x < y:
        tmp = x
        x = y
        y = tmp

    if y == 0:
        return x

    while y != 0:
        t = y
        y = x % t
        x = t
    return x


def get_bit(num: int, pos: int) -> int:
    mask = 1 << pos
    bit = mask & num
    if bit > 0:
        return 1
    else:
        return 0


def gauss_elimination(matrix: list[int], m: int, n: int) -> tuple[list[int], list[int], list[int]]:
    free_columns = []
    pivot_columns = []
    topmost_row = 0

    for column in range(m):
        # Checking if the pivot is not in the topmost row
        if get_bit(matrix[topmost_row], m - 1 - column) != 1:
            # If it's not, we are searching for the row with the leading 1 in the specific column
            found = False
            for row in range(topmost_row + 1, n):
                if get_bit(matrix[row], m - 1 - column) == 1:
                    # print("Pivot found at: ", row, "and is changed with ", topmost_row)
                    found = True
                    temp = matrix[topmost_row]
                    matrix[topmost_row] = matrix[row]
                    matrix[row] = temp
                    break

            # if no pivot found it means that the column has only zeros, it's a free column
            if not found:
                # print("Pivot not found, all zeros in ", column)
                free_columns.append(column)
                continue

        # print("Pivot found in ", column, "at", topmost_row)

        for row in range(n):
            if row == topmost_row:
                continue
            else:
                if get_bit(matrix[row], m - 1 - column) == 1:
                    # print("We are XORing ", row, "with", topmost_row)
                    matrix[row] ^= matrix[topmost_row]
                    # print("with result ", matrix[row])

        pivot_columns.append(column)
        topmost_row += 1
        # print("Next topmost row is: ", topmost_row)

        if topmost_row == n:
            break


    return matrix, pivot_columns, free_columns


def find_nullspace(matrix: list[int], pivot_cols: list[int], free_cols: list[int], m: int) -> list[int]:
    nullspace = [0] * len(free_cols)
    pos = 0

    for free_col in free_cols:
        basis_vec = 1 << (m - free_col - 1)

        rows = len(matrix) if len(matrix) < free_col else free_col

        for i in range(rows):
            if matrix[i] & basis_vec:
                nullspace[pos] ^= (1 << (m - pivot_cols[i] - 1))

        nullspace[pos] ^= basis_vec
        pos += 1

    return nullspace

def compute_square_root(num: int, factor_base: list[int]) -> int:
    square_num = 1

    for p in factor_base:
        factor = p * p
        while num % factor == 0:
            num //= factor
            square_num *= p

    return square_num
