#!/usr/bin/python3

import functools


def get_sign(x, ignore=False):
    """
    Get +/- sign from integer x
    """
    if x < 0:
        return "- "
    elif ignore:
        return ""
    return "+ "


def is_prime(n):
    """
    Check if n is prime

    Try to divide with 6n + 1 & 6n + 5
    Else are covered by 2 & 3
    """
    # Ignore negative, 0 & 1
    if n <= 1:
        return False
    # 2 & 3 are primes
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


def check_irr(poly, primes, px):
    """
    Check irreducibility of a polynom with prime factors of `a`
    """
    if poly[0] != functools.reduce(lambda x, y: x * y, primes):
        return False
    for p in primes:
        if not (poly.compute(p) % px) or not (poly.compute(-p) % px):
            return False
    return True


def egcd(a, b):
    """
    Compute extended euclidean for `a` and `b`

    Pre-condition: a > b
    """
    if a % b == 0:
        return (None, None)
    mem = [0, 1, 0, None]
    while b != 1:
        t = mem[1]
        mem[1] = mem[0] - t * (a // b)
        mem[0] = t

        if mem[3] is None:
            mem[3] = 1
        else:
            t = mem[3]
            mem[3] = mem[2] - t * (a // b)
            mem[2] = t

        t = b
        b = a % b
        a = t
    return (mem[3], mem[1])
