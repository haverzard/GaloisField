import functools


def get_sign(x, ignore=False):
    if x < 0:
        return "- "
    elif ignore:
        return ""
    return "+ "


def is_prime(n):
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 6
    while pow(i, 2) <= n:
        if n % (i - 1) == 0 or n % (i + 1) == 0:
            return False
        i += 6
    return True


def compute_poly(x, poly):
    total = 0
    for i in range(len(poly) - 1, -1, -1):
        total = total * x + poly[i]
    return total


def check_irr(poly, primes):
    if poly[0] != functools.reduce(lambda x, y: x * y, primes):
        return False
    for p in primes:
        if not compute_poly(p, poly) or not compute_poly(-p, poly):
            return False
    return True


def egcd(a, b):
    if a % b == 0:
        return 0
    mem = [0, 1]
    while b != 1:
        t = mem[1]
        mem[1] = mem[0] - t * (a // b)
        mem[0] = t
        t = b
        b = a % b
        a = t
    return mem[1]
