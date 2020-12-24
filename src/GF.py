#!/usr/bin/python3
#
# GF - Galois Field in Python3
# by haverzard (https://github.com/haverzard)

from util import get_sign, is_prime, compute_poly, check_irr, egcd
from exceptions import FFOperationException, PrimeFieldNoFitException
from fast_polynom import FastPolynom


class GF:
    """
    Galois Field class implementation - GF(p^m)

    Attributes:
    p - prime number
    m - positive integer (default: 1)
    irr - tuple consists of [0] irreducible polynom & [1] prime (used for extension field)
          set irr[0] with FastPolynom object
          set irr[1] as None to ignore reducibility check
    prime_check - enable/disable prime check for optimization purpose (default: True)

    Supports:
    - Prime Field (m = 0)
    - Extension Field (m != 1)
    """

    def __init__(self, p, m=1, irr=(None, None), prime_check=True):
        """
        Init finite field
        1. Check if p is prime
        2. Check if m is positive integer
        3. Check if type of field
        4. Check polynom's degree if field is extension
           and irreducibility of the polynom if primes are set
        """
        if prime_check:
            assert is_prime(p), "p must be prime"
        assert m > 0, "m must be positive"
        self.p = p
        self.m = m

        if m != 1:
            assert (
                irr[0].get_max_degree() == m
            ), "irreducible polynom is not for the finite field"

            irr[0].broadcast_modulo(p)
            if irr[1]:
                assert len(irr[1]), "make sure it's array"
                for x in irr[1]:
                    assert x == 1 or is_prime(x), "not prime"

                assert check_irr(irr[0], irr[1]), "polynom is reducible"
        self.irr = irr

    def __eq__(self, x):
        return self.p == x.p and self.m == x.m and self.irr[0] == x.irr[0]

    def __str__(self):
        """
        Transform finite field in universal format:
        `GF(p^m)[X] / F(x)`
        with F(x) is the irreducible polynom
        """
        p = self.p
        m = self.m
        return "GF({}{}){}".format(
            p,
            "^" + str(m) if m else "",
            "[X] / " + str(self.irr[0]) if self.irr[0] else "",
        )


class FFElement:
    """
    Galois Field Element/Member class implementation
    Data structure are all in polynomials (using FastPolynom object)

    Attributes:
    ff - GF object
    container - FastPolynom object

    Supports:
    - Addition
    - Subtraction
    - Multiplication
    - Inverse
    - Division
    - Modulo
    """

    def __init__(self, ff, container=None):
        self.ff = ff
        if ff.m:
            if container:
                self.container = container

                self.container.broadcast_modulo(ff.p)

                if self.container.get_max_degree() > ff.m:
                    if ff.m != 1:
                        self._fit()
                    else:
                        raise PrimeFieldNoFitException()
            else:
                self.container = FastPolynom()

    def _div(self, pol1, pol2, p, copy=True):
        d1 = pol1.get_max_degree()
        d2 = pol2.get_max_degree()
        if d1 < d2:
            return [0], pol1
        if d1 == 0 or d2 == 0:
            return pol1, [0]

        keys_pol2 = pol2.get_keys(rev=True)
        num = pol1.deepcopy()
        divs = FastPolynom()
        key = pol2[d2]
        while d1 >= d2:
            x = num[d1]
            assert x % key == 0, "bad reduce"
            fac = x // key
            divs[d1 - d2] = fac
            if fac:
                for j in keys_pol2:
                    num[d1 - d2 + j] -= (pol2[j] * fac) % p
                    num[d1 - d2 + j] %= p
            d1 = num.get_max_degree()
        return (divs, num)

    def __add__(self, x):
        try:
            assert self.ff == x.ff, "x is not in the same finite field"
            res = FFElement(self.ff)
            for i in range(self.ff.m):
                res.container[i] = (self.container[i] + x.container[i]) % self.ff.p
            return res
        except AttributeError:
            raise FFOperationException("+", "x is not FFElement object?")

    def __sub__(self, x):
        try:
            assert self.ff == x.ff, "x is not in the same finite field"
            res = FFElement(self.ff)
            for i in range(self.ff.m):
                res.container[i] = (self.container[i] - x.container[i]) % self.ff.p
            return res
        except AttributeError:
            raise FFOperationException("-", "x is not FFElement object?")

    def __mul__(self, x):
        try:
            assert self.ff == x.ff, "x is not in the same finite field"
            deg = self.ff.m
            p = self.ff.p
            res = FastPolynom()
            for i in range(deg):
                for j in range(deg):
                    res[i + j] += (self.container[i] * x.container[j]) % p
                    res[i + j] %= p
            if ff.m != 1:
                _, res = self._div(res, self.ff.irr[0], p)
            return FFElement(self.ff, res)
        except AttributeError:
            raise FFOperationException("*", "x is not FFElement object?")

    def __floordiv__(self, x):
        try:
            assert self.ff == x.ff, "x is not in the same finite field"
            if self.ff.m == 1:
                res = FastPolynom()
                res[0] = self.container[0] // x.container[0]
            else:
                res, _ = self._div(self.container, x.container, self.ff.p)
            return FFElement(self.ff, res)
        except AttributeError:
            raise FFOperationException("//", "x is not FFElement object?")

    def __mod__(self, x):
        try:
            assert self.ff == x.ff, "x is not in the same finite field"
            if self.ff.m == 1:
                res = FastPolynom(self.ff)
                res.container[0] = self.container[0] % x.container[0]
            else:
                d = self.container.get_max_degree()
                pre = {}
                a = 2
                _, x1 = self._div(FastPolynom({1: 1}), x.container, self.ff.p)
                pre[1] = FFElement(self.ff, x1)
                _, x2 = self._div(FastPolynom({2: 1}), x.container, self.ff.p)
                pre[2] = FFElement(self.ff, x2)
                while a * 2 <= d:
                    _, x2 = self._div(
                        (pre[a] * pre[a]).container, x.container, self.ff.p
                    )
                    a *= 2
                    pre[a] = FFElement(self.ff, x2)
                res = FFElement(self.ff)
                for i in self.container.get_keys(rev=True):
                    temp = FFElement.gen_one(self.ff)
                    j = 1
                    while i:
                        if i & 1:
                            temp *= pre[j]
                        i >>= 1
                        j *= 2
                    res += temp
            return res
        except AttributeError:
            raise FFOperationException("%", "x is not FFElement object?")

    @staticmethod
    def _is_one(a):
        return a.container.get_max_degree() == 0

    @staticmethod
    def _is_zero(a):
        return a.container.get_max_degree() == -1

    @staticmethod
    def gen_one(ff):
        return FFElement(ff, FastPolynom({0: 1}))

    @staticmethod
    def _egcd(a, b):
        if FFElement._is_zero(a % b):
            raise Exception("a & b must be co-prime")
        mem = [FFElement.gen_zero(a.ff), FFElement.gen_one(a.ff)]
        while not FFElement._is_one(b):
            t = mem[1]
            mem[1] = mem[0] - t * (a // b)
            mem[0] = t
            t = b
            b = a % b
            a = t
        return mem[1]

    def inverse(self):
        try:
            if self.ff.m == 1:
                res = FFElement(self.ff)
                res.container[0] = egcd(self.ff.p, self.container[0])[1]
            else:
                irr = FFElement(self.ff, self.container)
                irr.container = self.ff.irr[0]

                res = FFElement._egcd(irr, self)
            return res
        except AttributeError:
            raise FFOperationException("^-1", "Something went wrong?")

    def _fit(self):
        try:
            x = FFElement.gen_zero(self.ff)
            x.container = self.ff.irr[0]
            self.container = (self % x).container
        except AttributeError:
            raise FFOperationException("fit", "Something went wrong?")

    def __str__(self):
        return "{}: {}".format(self.ff, self.container)


if __name__ == "__main__":
    ff = GF(2, 3, ([1, 1, 0, 1], [1]))
    print(FFElement(ff, [1, 1, 1]) * FFElement(ff, [0, 0, 1]))
    print(FFElement(ff, [1, 1, 1]).inverse())
    print(FFElement(ff, FastPolynom({0: 1, 1: 1, 2: 1})).inverse())
    ff = GF(2, 1)
    print(FFElement(ff, FastPolynom({0: 1})).inverse())
