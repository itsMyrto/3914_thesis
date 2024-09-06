import time
from math import ceil, log2, sqrt, log, e, exp
from utils import sieve_of_eratosthenes, tonelli_shanks, legendre_symbol, find_nullspace, gauss_elimination, compute_square_root, gcd

def get_bounds(n: int) -> tuple[int, int, int]:
    bound = sqrt((log(n)*log(log(n))))
    digits = len(str(n))
    if digits < 10:
        return round(exp(bound)), round(exp(bound)), 1000
    if digits < 40:
        return round(exp(bound) ** 0.6), round(exp(bound)), 10000
    else:
        return round(exp(bound) ** 0.5), round(exp(bound)), 50000

def generate_factor_base(n: int, B: int) -> list[int]:
    factor_base = []
    primes = sieve_of_eratosthenes(B)

    for p in primes:
        if legendre_symbol(n, p) == 1:
            factor_base.append(p)

    return factor_base

def precompute_logs(factor_base: list[int]) -> list[float]:
    logs = []
    for factor in factor_base:
        logs.append(log2(factor))
    return logs

def precompute_roots(factor_base: list[int], n: int) -> list[tuple[int, int]]:
    roots = []
    for p in factor_base:
        roots.append(tonelli_shanks(n, p))
    return roots

def is_b_smooth(n: int, factor_base: list[int]) -> bool:
    for p in factor_base:
        while n % p == 0:
            n //= p
    if n == 1 or n == -1:
        return True
    return False

def factor_base_product(factor_base: list[int]) -> int:
    prod = 1
    for p in factor_base:
        prod *= p
    return prod

def is_b_smooth2(num: int, prod: int) -> bool:
    while True:
        g = gcd(num, prod)
        if g > 1:
            counter = 0
            n = num
            while n % g == 0:
                counter += 1
                n = n // g
            r = num // pow(g, counter)
            if r == 1:
                return True
            num = r
        else:
            return False


def construct_b_smooths(n: int, sieve_start: int, sieve_end: int, step: int, prime_logs: list[float], factor_base: list[int], roots: list[tuple[int, int]], prod: int, threshold: int) -> tuple[list[int], list[int]]:
    Q = []
    Q_log = []
    n_root = ceil(sqrt(n))

    for x in range(sieve_start, sieve_end, step):
        value = abs((x + n_root) ** 2 - n)
        Q.append(value)
        Q_log.append(round(log2(value)))

    for i in range(1, len(factor_base)):
        current_prime = factor_base[i]
        current_prime_log = prime_logs[i]
        for r in roots[i]:
            pos = (r - n_root - sieve_start) % current_prime
            if step == -1:
                if pos < sieve_start - sieve_end:
                    Q_log[pos] -= current_prime_log
                    pos = (-1) * (pos - current_prime)
                for j in range(pos, sieve_start - sieve_end, current_prime):
                    Q_log[j] -= current_prime_log
            elif step == 1:
                for j in range(pos, sieve_end - sieve_start, current_prime):
                    Q_log[j] -= current_prime_log

    Q_smooth = []
    x_vals = []
    for i in range(len(Q_log)):
        if Q_log[i] < threshold:
            if is_b_smooth2(Q[i], prod):
                Q_smooth.append(step * Q[i])
                x_vals.append(step * i + n_root + sieve_start)

    return Q_smooth, x_vals

def build_exponent_matrix(smooth_numbers: list[int], factor_base: list[int]) -> list[int]:
    matrix = [0] * (len(factor_base) + 1)
    pos = 0

    for smooth_number in smooth_numbers:
        matrix[pos] <<= 1
        if smooth_number < 0:
            matrix[pos] |= 1
        else:
            matrix[pos] |= 0

    pos += 1

    for p in factor_base:
        for smooth_number in smooth_numbers:
            num = abs(smooth_number)
            counter = 0
            while num % p == 0:
                counter += 1
                num //= p
            counter %= 2
            matrix[pos] <<= 1
            matrix[pos] |= counter
        pos += 1

    return matrix


def quadratic_sieve(n: int):

    B, M, SIEVE_LEN = get_bounds(n)

    sieve_start = 0
    sieve_end = min(SIEVE_LEN, M)

    smooth_numbers = []
    x_values = []

    factor_base = generate_factor_base(n, B)
    factor_base_logs = precompute_logs(factor_base)
    roots = precompute_roots(factor_base, n)
    product = factor_base_product(factor_base)

    RELATIONS = len(factor_base) + 10

    while len(smooth_numbers) < RELATIONS:

        temp1, temp2 = construct_b_smooths(n, sieve_start, sieve_end, 1, factor_base_logs, factor_base, roots,  product, 20)
        smooth_numbers.extend(temp1)
        x_values.extend(temp2)

        temp1, temp2 = construct_b_smooths(n, -sieve_start, -sieve_end, -1, factor_base_logs, factor_base, roots, product, 20)
        smooth_numbers.extend(temp1)
        x_values.extend(temp2)

        sieve_start = sieve_end
        sieve_end = min(sieve_end + SIEVE_LEN, M)

        if sieve_start == sieve_end:
            break

    m = len(smooth_numbers)

    s = len(factor_base) + 1

    matrix = build_exponent_matrix(smooth_numbers, factor_base)

    matrix, pivots, free = gauss_elimination(matrix, m, s)

    nullspace = find_nullspace(matrix, pivots, free, m)

    for basis_vec in nullspace:
        x_val = 1
        y = 1
        for i in range(m):
            pos = 1 << (m - i - 1)
            if basis_vec & pos:
                y *= smooth_numbers[i]
                x_val *= x_values[i]

        y = compute_square_root(y, factor_base)

        if x_val % n != y % n and x_val % n != ((-1) * y) % n:
            factor1 = gcd(x_val - y, n)
            if factor1 != 1:
                factor2 = n // factor1
                return abs(factor1), abs(factor2)

    return None, None