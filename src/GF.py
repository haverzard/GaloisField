from util import get_sign, is_prime, compute_poly, check_irr
from exceptions import FFOperationException


class GF:
    def __init__(self, p, m=1, irr=(None, None)):
        assert is_prime(p), "p must be prime"
        assert m > 0, "m must be positive"
        self.p = p
        self.m = m

        if m != 1:
            assert (
                len(irr[0]) == m + 1
            ), "irreducible polynom is not for the finite field"
            assert len(irr[1]), "make sure it's array"

            for x in irr[1]:
                assert x == 1 or is_prime(x), "not prime"

            for i in range(len(irr[0])):
                irr[0][i] %= p

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
    def __init__(self, ff, container=None):
        self.ff = ff
        if ff.m:
            if container:
                assert len(container) == ff.m, "not in the finite field"
                self.container = container
                for i in range(ff.m):
                    self.container[i] %= ff.p
            else:
                self.container = [0 for _ in range(ff.m)]

    def _length(self):
        i = self.ff.m - 1
        while i >= 0 and self.container[i] != 1:
            i -= 1
        return i

    def _div(self, pol1, pol2, p, copy=True):
        d1 = len(pol1)
        while d1 > 0 and pol1[d1 - 1] != 1:
            d1 -= 1
        d2 = len(pol2)
        while d2 > 0 and pol2[d2 - 1] != 1:
            d2 -= 1
        if d1 < d2:
            return [0], pol1
        if d1 == 0 or d2 == 0:
            return pol1, [0]
        d3 = d1 - d2 + 1

        divs = [0 for _ in range(d3)]
        num = pol1[:d1]
        divisor = pol2[:d2]
        key = divisor[-1]
        for i in range(d1 - 1, d2 - 2, -1):
            x = num[i]
            assert x % key == 0, "bad reduce"
            fac = x // key
            divs[-d1 + i] = fac
            if fac:
                for j in range(d2):
                    num[i - j] -= (divisor[-j - 1] * fac) % p
                    num[i - j] %= p
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
            res = [0 for _ in range(deg * 2)]
            for i in range(deg):
                for j in range(deg):
                    res[i + j] += (self.container[i] * x.container[j]) % p
                    res[i + j] %= p
            _, x = self._div(res, self.ff.irr[0], p)
            x = x[: self.ff.m]
            return FFElement(self.ff, x + [0 for _ in range(self.ff.m - len(x))])
        except AttributeError:
            raise FFOperationException("+", "x is not FFElement object?")

    def __floordiv__(self, x):
        try:
            assert self.ff == x.ff, "x is not in the same finite field"
            div, _ = self._div(self.container, x.container, self.ff.p)
            return FFElement(self.ff, div + [0 for _ in range(self.ff.m - len(div))])
        except AttributeError:
            raise FFOperationException("//", "x is not FFElement object?")

    def __mod__(self, x):
        try:
            assert self.ff == x.ff, "x is not in the same finite field"
            _, rem = self._div(self.container, x.container, self.ff.p)
            rem = rem[: self.ff.m]
            return FFElement(self.ff, rem + [0 for _ in range(self.ff.m - len(rem))])
        except AttributeError:
            raise FFOperationException("//", "x is not FFElement object?")

    @staticmethod
    def _is_one(a):
        return a._length() == 0

    @staticmethod
    def _is_zero(a):
        return a._length() == -1

    @staticmethod
    def gen_zero(ff):
        return FFElement(ff, [0 for _ in range(ff.m)])

    @staticmethod
    def gen_one(ff):
        return FFElement(ff, [1] + [0 for _ in range(ff.m - 1)])

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

    def __str__(self):
        elements = ""
        e = self.container[-1]
        if e:
            elements = "{}{}x^{} ".format(
                get_sign(e, True), e if e != 1 else "", self.ff.m - 1
            )
        for i in range(self.ff.m - 2, 0, -1):
            e = self.container[i]
            if e:
                elements += "{}{}x^{} ".format(
                    get_sign(e, not elements), e if e != 1 else "", i
                )
        e = self.container[0]
        if e or not elements:
            elements += "{}{}".format(get_sign(e, not elements), e)

        return "{}: {}".format(self.ff, elements)


if __name__ == "__main__":
    ff = GF(2, 3, ([1, 1, 0, 1], [1]))
    print(FFElement(ff, [1, 1, 1]) * FFElement(ff, [0, 0, 1]))
    print(FFElement(ff, [1, 1, 1]).inverse())
