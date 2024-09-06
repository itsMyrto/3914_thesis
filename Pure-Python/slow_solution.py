# def b_bound(n: int) -> int:
#     factor = sqrt((log(n)*log(log(n))))
#     if len(str(n)) < 10:
#         return round(exp(factor))
#     elif len(str(n)) < 40:
#         return round(exp(factor) ** 0.6)
#     else:
#         return round(exp(factor) ** 0.5)
#
#
# def m_bound(n: int) -> int:
#     factor = sqrt(log(n) * log(log(n)))
#     factor = pow(e, factor)
#     return round(factor)
#
#
# def sieve_len(n: int) -> int:
#     size = len(str(n))
#     if size < 10:
#         return 1000
#     elif size < 30:
#         return 10000
#     else:
#         return 50000
#
# def generate_factor_base(n: int, B: int) -> list[int]:
#     factor_base = []
#     primes = sieve_of_eratosthenes(B)
#
#     for p in primes:
#         if legendre_symbol(n, p) == 1:
#             factor_base.append(p)
#
#     return factor_base
#
# def filter_candidates(candidates: list[int], x_vals: list[int], factor_base: list[int]) -> tuple[list[int], list[int]]:
#     b_smooths = []
#     x_values = []
#
#     for i in range(len(candidates)):
#         candidate = candidates[i]
#
#         for p in factor_base:
#             while candidate % p == 0:
#                 candidate //= p
#
#         if candidate == 1 or candidate == -1:
#             b_smooths.append(candidates[i])
#             x_values.append(x_vals[i])
#
#     return b_smooths, x_values
#
# def evaluate_polynomial(n: int, n_root: int, sieve_start: int, sieve_end: int, step: int) -> tuple[list[int], list[int]]:
#     q = []
#     q_log = []
#
#     for x in range(sieve_start, sieve_end, step):
#         q_val = abs((x + n_root) ** 2 - n)
#         q.append(q_val)
#         q_log.append(round(log2(q_val)))
#
#     return q, q_log
#
#
# def construct_negative_b_smooths(n: int, factor_base: list[int], sieve_start: int, sieve_end: int, threshold: int) -> tuple[list[int], list[int]]:
#     n_root = ceil(sqrt(n))
#     threshold_candidates = []
#     x_values = []
#
#     candidates, logs = evaluate_polynomial(n, n_root, sieve_start, sieve_end, -1)
#
#     for p in factor_base[1:]:
#         r1, r2 = tonelli_shanks(n, p)
#         logp = log2(p)
#         for r in (r1, r2):
#             pos = (r - n_root - sieve_start) % p
#             if pos < sieve_start - sieve_end:
#                 logs[pos] -= logp
#                 pos -= p
#                 pos *= (-1)
#                 while pos < sieve_start - sieve_end:
#                     logs[pos] -= logp
#                     pos += p
#
#     for i in range(len(logs)):
#         if logs[i] < threshold:
#             threshold_candidates.append((-1) * candidates[i])
#             x_values.append(-i + sieve_start + n_root)
#
#     return filter_candidates(threshold_candidates, x_values, factor_base)
#
#
# def construct_positive_b_smooths(n: int, factor_base: list[int], sieve_start: int, sieve_end: int, threshold: int) -> tuple[list[int], list[int]]:
#     n_root = ceil(sqrt(n))
#     threshold_candidates = []
#     x_values = []
#
#     candidates, logs = evaluate_polynomial(n, n_root, sieve_start, sieve_end, 1)
#
#     for p in factor_base[1:]:
#         r1, r2 = tonelli_shanks(n, p)
#         logp = log2(p)
#         for r in (r1, r2):
#             pos = (r - n_root - sieve_start) % p
#             while pos < sieve_end - sieve_start:
#                 logs[pos] -= logp
#                 pos += p
#
#     for i in range(len(logs)):
#         if logs[i] < threshold:
#             threshold_candidates.append(candidates[i])
#             x_values.append(i + sieve_start + n_root)
#
#     return filter_candidates(threshold_candidates, x_values, factor_base)
#
#
# def build_exponent_matrix(smooth_numbers: list[int], factor_base: list[int]) -> list[int]:
#     matrix = [0] * (len(factor_base) + 1)
#     pos = 0
#
#     for smooth_number in smooth_numbers:
#         matrix[pos] <<= 1
#         if smooth_number < 0:
#             matrix[pos] |= 1
#         else:
#             matrix[pos] |= 0
#
#     pos += 1
#
#     for p in factor_base:
#         for smooth_number in smooth_numbers:
#             num = abs(smooth_number)
#             counter = 0
#             while num % p == 0:
#                 counter += 1
#                 num //= p
#             counter %= 2
#             matrix[pos] <<= 1
#             matrix[pos] |= counter
#         pos += 1
#
#     return matrix
#
#
# def quadratic_sieve(n: int):
#
#     M = m_bound(n)
#     B = b_bound(n)
#
#     SIEVE_LEN = sieve_len(n)
#
#     sieve_start = 0
#     sieve_end = min(SIEVE_LEN, M)
#
#     smooth_numbers = []
#     x_values = []
#
#     factor_base = generate_factor_base(n, B)
#
#     RELATIONS = len(factor_base) + 10
#
#     while len(smooth_numbers) < RELATIONS:
#
#         temp1, temp2 = construct_positive_b_smooths(n, factor_base, sieve_start, sieve_end, 10)
#         smooth_numbers.extend(temp1)
#         x_values.extend(temp2)
#
#         temp1, temp2 = construct_negative_b_smooths(n, factor_base, -sieve_start, -sieve_end, 10)
#         smooth_numbers.extend(temp1)
#         x_values.extend(temp2)
#
#         sieve_start = sieve_end
#         sieve_end = min(sieve_end + SIEVE_LEN, M)
#
#         if sieve_start == sieve_end:
#             break
#
#     m = len(smooth_numbers)
#     s = len(factor_base) + 1
#
#     matrix = build_exponent_matrix(smooth_numbers, factor_base)
#
#     matrix, pivots, free = gauss_elimination(matrix, m, s)
#
#     nullspace = find_nullspace(matrix, pivots, free, m)
#
#     for basis_vec in nullspace:
#         x_val = 1
#         y = 1
#         for i in range(m):
#             pos = 1 << (m - i - 1)
#             if basis_vec & pos:
#                 y *= smooth_numbers[i]
#                 x_val *= x_values[i]
#
#         y = compute_square_root(y, factor_base)
#
#         if (x_val % n) != (y % n) and (x_val % n) != (((-1) * y) % n):
#             factor1 = gcd(x_val - y, n)
#             if factor1 != 1:
#                 factor2 = n // factor1
#                 return abs(factor1), abs(factor2)
#
#     return None, None
