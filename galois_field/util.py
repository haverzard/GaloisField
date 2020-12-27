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

    Using fermat's little theorem
    a^(n-1) = 1 + k*n
    """
    # Ignore negative, 0 & 1
    if n <= 1:
        return False
    # 2 & 3 are primes
    if n <= 3:
        return True
    return pow(2, n - 1, n) == 1


def check_irr(fe):
    """
    Check irreducibility of a polynom with prime factors of `a`

    Using euler's theorem
    a^phi = 1 + k*p
    """
    return fe.pow(fe.ff.p ** fe.ff.m - 1).is_one()


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
