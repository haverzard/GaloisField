from util import get_sign, is_prime, compute_poly, check_irr
from exceptions import FFOperationException
from fast_polynom import FastPolynom


class GF:
    def __init__(self, p, m=1, irr=(None, None), prime_check=True):
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
        p = self.p
        m = self.m
        elements = ""
        if self.irr:
            poly = self.irr[0]
            e = poly[-1]
            if e:
                elements = "{}{}x^{} ".format(get_sign(e, True), e if e != 1 else "", m)
            for i in range(m - 1, 0, -1):
                e = poly[i]
                if e:
                    elements += "{}{}x^{} ".format(
                        get_sign(e, not elements), e if e != 1 else "", i
                    )
            e = poly[0]
            if e or not elements:
                elements += "{}{}".format(get_sign(e, not elements), e)

        return "GF({}{}){}".format(p, "^" + str(m) if m else "", " | " + elements)


class FFElement:
    def __init__(self, ff, container=None, fit=False, optimize=True):
        self.ff = ff
        self.optimize = optimize
        if ff.m:
            if container:
                self.container = container

                self.container.broadcast_modulo(ff.p)

                if self.container.get_max_degree() > ff.m:
                    self._fit()
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
        return self.__add__(x)

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
            _, x = self._div(res, self.ff.irr[0], p)
            return FFElement(self.ff, x)
        except AttributeError:
            raise FFOperationException("*", "x is not FFElement object?")

    def __floordiv__(self, x):
        try:
            assert self.ff == x.ff, "x is not in the same finite field"
            div, _ = self._div(self.container, x.container, self.ff.p)
            return FFElement(self.ff, div)
        except AttributeError:
            raise FFOperationException("//", "x is not FFElement object?")

    def __mod__(self, x):
        try:
            assert self.ff == x.ff, "x is not in the same finite field"
            _, rem = self._div(self.container, x.container, self.ff.p)
            rem = rem[: self.ff.m]
            return FFElement(self.ff, rem + [0 for _ in range(self.ff.m - len(rem))])
        except AttributeError:
            raise FFOperationException("%", "x is not FFElement object?")

    @staticmethod
    def _is_one(a):
        return a.container.get_max_degree() == 0

    @staticmethod
    def _is_zero(a):
        return a.container.get_max_degree() == -1

    @staticmethod
    def gen_zero(ff):
        return FFElement(ff, FastPolynom())

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
            x = FFElement(self.ff, self.container)
            x.container = self.ff.irr[0]

            inv = FFElement._egcd(x, self)
            return inv
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
