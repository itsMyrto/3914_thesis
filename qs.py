from math import ceil, log2, exp, sqrt
import time

def sieve_of_eratosthenes(B: int) -> list[int]:
    """
    Generating all primes up to a bound B using the sieve of eratosthenes
    pseudocode and examples can be found in https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes

    :param B: an integer, where B>1
    :return: all prime numbers from 2 through B
    """
    primes = [True] * (B + 1)
    primes[0] = primes[1] = False

    root = ceil(pow(B, 0.5))

    for i in range(2, root):
        if primes[i]:
            for j in range(i * i, B + 1, i):
                primes[j] = False

    return [i for i in range(0, len(primes)) if primes[i]]


def legendre_symbol(n: int, p: int) -> int:
    """
    This function uses Euler's criterion to check if n is a quadratic residue modulo a prime number
    n is a quadratic residue modulo p if and only if n^((p-1)/2) ≡ 1 (mod p).

    :param n: The number to factor.
    :param p: The prime number modulo which we are testing n.
    :return: 1 if n is a quadratic residue modulo prime
    """
    return pow(n, (p - 1) // 2, p)



def generate_factor_base(n: int, B: int) -> list[int]:
    """
    This function generates the factor base which contains only the prime numbers up to B
    that are quadratic residues.

    :param n: The number to factor
    :param B: The bound B to find primes
    :return: A list of integers, the prime numbers that create the factor base
    """
    factor_base = []
    primes = sieve_of_eratosthenes(B)

    for p in primes:
        if legendre_symbol(n, p) == 1:
            factor_base.append(p)

    return factor_base



def tonelli_shanks(n: int, p: int) -> tuple[int, int]:
    """
    In this function the tonelli-shanks algorithm is implemented to solve for r
    in a congruence of the form r^2 ≡ n (mod p), where p is a prime, and it finds a square root
    of n modulo p
    pseudocode and examples can be found in https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm

    :param n: An integer of Z/pZ such that solutions to the congruence r^2 ≡ n (mod p) exists
    :param p: A prime
    :return: The two roots r and p-r
    """
    S = 0
    Q = p - 1
    while Q % 2 == 0:
        Q = Q // 2
        S += 1

    if p % 4 == 3:
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


def filter_b_smooth_candidates(candidates: list[int], vals: list[int], factor_base: list[int]) -> tuple[list[int], list[int]]:
    """
    Filters the given list of candidate numbers, returning only those that are B-smooth.

    :param candidates: A list of numbers that are likely to be B-smooth
    :param vals: The corresponding x's
    :param factor_base: A list of prime numbers up to B
    :return: A list of numbers that are B-smooth and the values x
    """
    b_smooth_numbers = []
    x = []
    for i in range(len(candidates)):
        num = candidates[i]
        for p in factor_base:
            while num % p == 0:
                num //= p
        if num == 1:
            b_smooth_numbers.append(candidates[i])
            x.append(vals[i])
    return b_smooth_numbers, x


def get_b_smooth_numbers(n: int, factor_base: list[int], sieve_start: int, sieve_end: int) -> tuple[list[int], list[int]]:
    """
    This function generates the B-smooth numbers using the logarithm approximation.
    The steps:
    1. Calculate the values of Q(xi) = (square_root(n) + x)^2 - n and the logarithm of Q(xi), for xi [M,N] with M<N and M,N>0, which is the sieving subinterval
    2. For each prime number in the factor base (excluding 2), compute the roots of x^2 ≡ n (mod p) using the Tonelli-Shanks algorithm
    3. For each root, subtract from the logarithms of Q(xi) the value log(p) (p the prime), for xi = r + i*p, where r root and i integer between 1 and N
    4. Set a threshold and for each Q(xi) value test if is less than the threshold. If so, the Q(xi) value is considered to be B-smooth candidate
    5. After collecting all the B-smooth candidates, filter them and keep only the B-smooths

    This method speeds up by far the computational time, approximation is all we need
    explanation and example can be found in https://www.cs.umd.edu/users/gasarch/COURSES/456/F19/lecqs/lecqs.pdf in slides 169 and after
    as well as in thesis https://dspace.cvut.cz/bitstream/handle/10467/94585/F8-DP-2021-Vladyka-Ondrej-DP_Vladyka_Ondrej_2021.pdf?sequence=-1&isAllowed=y

    :param n: The number to be factored
    :param factor_base: An list of prime numbers from 2 up to B
    :param sieve_start: The starting of the sieving subinterval
    :param sieve_end: The ending of the sieving subinterval
    :return: The values B-smooth Q(xi) and the values xi
    """
    q_xi = []
    q_xi_log = []
    n_root = ceil(sqrt(n))

    for x in range(sieve_start, sieve_end):
        num = abs((x + n_root) ** 2 - n)
        q_xi.append(num)
        q_xi_log.append(round(log2(num)))


    for p in factor_base[1:]:
        r1, r2 = tonelli_shanks(n, p)
        logp = log2(p)
        for r in (r1, r2):
            pos = (r - n_root - sieve_start) % p
            while pos < sieve_end - sieve_start:
                q_xi_log[pos] -= logp
                pos += p

    threshold = 10

    b_smooth_candidates = []
    x = []

    for i in range(len(q_xi_log)):
        if q_xi_log[i] < threshold:
            b_smooth_candidates.append(q_xi[i])
            x.append(i + sieve_start + n_root)

    return filter_b_smooth_candidates(b_smooth_candidates, x, factor_base)


def build_exponent_matrix(smooth_numbers: list[int], factor_base: list[int]) -> list[int]:
    """
    This function builds the matrix that contains all the exponent vectors (in GF2) of smooth numbers
    In order to save space and work with large numbers, the binary rows of the matrix are saved as integers.

    :param smooth_numbers: A lists containing all the smooth numbers
    :param factor_base: A list with all the primes
    :return: The exponent matrix as a 1D list of integers
    """
    smooths = smooth_numbers[::-1]
    bit_vector = [0] * len(factor_base)
    pos = 0
    for p in factor_base:
        for smooth_number in smooths:
            num = smooth_number
            counter = 0
            while num % p == 0:
                num //= p
                counter += 1
            counter %= 2
            bit_vector[pos] = bit_vector[pos] << 1
            bit_vector[pos] = bit_vector[pos] | counter
        pos += 1

    return bit_vector


def get_bit(num: int, pos: int) -> int:
    """
    This functions finds the bit of an integer number in a specific position
    It creates a mask, which is a number with 1 in the pos position, and then it uses the
    bitwise and (&) operator to find if the bit is 1 or zero.

    :param num: An integer number
    :param pos: The position we want to find the bit
    :return: 1 or 0
    """
    mask = 1 << pos
    bit = mask & num
    if bit > 0:
        return 1
    else:
        return 0


def gcd(x: int, y: int) -> int:
    """
    This function implements the greatest common divisor of two numbers algorithm

    :param x: The first integer
    :param y: The second integer
    :return: The greatest common divisor of x and y
    """
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


def gauss_elimination(matrix: list[int], m: int, n: int) -> tuple[list[int], list[int], list[int]]:
    """
    This function implements gauss elimination in GF(2).
    The flow of the algorithm is:
    (1) Starts searching for a pivot in a column. If the pivot is not in the topmost row then it swaps the two rows.
        If there is no pivot then it continues to the next column
    (2) If a pivot is found, it performs XOR with all the rows that have 1 in the specific column

    :param matrix: The exponent matrix
    :param m: The number of smooth numbers
    :param n: The number of the factors in the factor base
    :return: The matrix after gauss elimination in rref form as a 1D list, the pivot columns and the free variable columns
    """
    topmost_row = 0
    free_cols = []
    pivot_cols = []
    for col in range(m):

        if get_bit(matrix[topmost_row], col) != 1:
            flag = False
            for row in range(topmost_row + 1, n):
                if get_bit(matrix[row], col) == 1:
                    temp = matrix[topmost_row]
                    matrix[topmost_row] = matrix[row]
                    matrix[row] = temp
                    flag = True
                    break

            if not flag:
                free_cols.append(col)
                continue

        for row in range(0, n):
            if row == topmost_row:
                continue
            else:
                if get_bit(matrix[row], col) == 1:
                    matrix[row] = matrix[topmost_row] ^ matrix[row]

        pivot_cols.append(col)
        topmost_row += 1


    return matrix, pivot_cols, free_cols


def null_space(matrix: list[int], pivot_cols: list[int], free_cols: list[int]) -> list[int]:
    """
    This function finds the null space of a matrix in rref form
    For every free column, it creates a basis factor with the only free variable corresponding to the column being 1
    Then, for every pivot row it performs the & operator to check if the pivot row has 1 in the free variable position
    If so, it adds the basis vector with the specific pivot row (XOR in GF2)

    :param matrix: A 1D list, the matrix in rref
    :param pivot_cols: A list containing the columns were the pivots are found in the matrix
    :param free_cols: A list containing the columns were the free variables are found in the matrix
    :return:
    """
    nullspace = [0] * len(free_cols)
    pos = 0
    for free_col in free_cols:

        basis_vec = 1 << free_col
        rows = len(matrix) if len(matrix) < free_col else free_col

        for i in range(rows):

            if matrix[i] & basis_vec:
                nullspace[pos] ^= (1 << pivot_cols[i])

        nullspace[pos] ^= basis_vec
        pos += 1

    return nullspace


def compute_square_root(num: int, factor_base: list[int]) -> int:
    """
    This function returns the square root of a number over the given factor base

    :param num: An integer to find the square root whose factors are in the factor base
    :param factor_base: A list of integers representing the factor base
    :return: The square root
    """
    square_num = 1
    for p in factor_base:
        while num % (p * p) == 0:
            num //= (p * p)
            square_num *= p
    return square_num


def quadratic_sieve(n: int):
    """
    This function implements the quadratic sieve algorithm. Steps:
    (1) Setting the bounds L, B and the subinterval length
    (2) Generate the relations
    (3) Build the exponent matrix
    (4) Perform Gauss Elimination
    (5) Find the null space of the matrix
    (6) Search solutions

    :param n: The number to factorize (odd number that is not a square)
    :return: The factors, if found
    """

    L = round(exp(sqrt((log2(n)*log2(log2(n))))))
    B = round(exp(0.5 * sqrt((log2(n)*log2(log2(n))))))
    # B = round(L ** 0.3) * 3

    SIEVE_LEN = 50000

    factor_base = generate_factor_base(n, B)
    smooth_numbers = []
    x = []

    sieve_start = 0
    sieve_end = min(SIEVE_LEN, L)
    RELATIONS = len(factor_base) + round(0.2 * len(factor_base))
    # RELATIONS = len(factor_base) + 10

    while len(smooth_numbers) < RELATIONS:

        smooth_nums, x_vals = get_b_smooth_numbers(n, factor_base, sieve_start, sieve_end)

        smooth_numbers.extend(smooth_nums)
        x.extend(x_vals)

        sieve_start = sieve_end
        sieve_end = min(sieve_end + SIEVE_LEN, L)


    matrix = build_exponent_matrix(smooth_numbers, factor_base)

    m = len(smooth_numbers)
    s = len(factor_base)

    matrix, pivots, free = gauss_elimination(matrix, m, s)

    nullspace = null_space(matrix, pivots, free)

    for basis_vec in nullspace:

        x_val = 1
        y = 1
        for i in range(m):
            pos = 1 << i
            if basis_vec & pos:
                y *= smooth_numbers[i]
                x_val *= x[i]

        y = compute_square_root(y, factor_base)

        if (x_val % n) != (y % n) and (x_val % n) != (((-1) * y) % n):
            factor1 = gcd(x_val - y, n)
            if factor1 != 1:
                factor2 = n // factor1
                return factor1, factor2

    return None, None



# -------------------------- MAIN PROGRAM -------------------------------

number_to_be_factored = 3744843080529615909019181510330554205500926021947

start = time.time()
f1, f2 = quadratic_sieve(number_to_be_factored)
end = time.time()

print("Time needed in seconds: ", end - start, " to process a number with ", len(str(number_to_be_factored)), " digits")

if f1 == f2 is None:
    print("No solution found")
else:
    print("The factors of ", number_to_be_factored, "are:", f1, "*", f2)
