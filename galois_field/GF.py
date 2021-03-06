#!/usr/bin/python3

from .util import get_sign, is_prime, check_irr, egcd
from .exceptions import FFOperationException, PrimeFieldNoFitException
from .fast_polynom import FastPolynom


class GF:
    """
    Galois Field class implementation - GF(p^m)

    Galois field is a field consists of range defined by a prime and its power (positive integer),
    and finite elements within that range.

    We call galois field with power equal to 1 as prime field, because it's just a prime.
    Elements of prime field are [0...p-1].

    We call galois field with power not equal to 1 as extension field.
    Elements of extension field are polynomials with maximum degree of (p-1).
    Extension field uses `prime` polynomial (or irreducible polynomial) instead for
    limiting the elements in the field (and also prime `p`).

    Attributes:
        p - prime number
        m - positive integer (default: 1)
        irr - irreducible polynom as FastPolynom object

    Usage:
        GF(2, 3, FastPolynom({0:1, 1:1, 3:1}))
    """

    def __init__(self, p, m=1, irr=None, prime_check=True):
        """
        Init finite field
        1. Check if p is prime
        2. Check if m is positive integer
        3. Check if type of field
        4. Check polynom's degree if field is extension and irreducibility of the polynom

        Arguments:
            p - prime number
            m - positive integer (default: 1)
            irr - irreducible polynom as FastPolynom object
            prime_check - enable/disable prime check for optimization purpose (default: True)

        Usage:
            GF(2, 3, FastPolynom({0:1, 1:1, 3:1}))
        """
        if prime_check:
            assert is_prime(p), "p must be prime"
        assert m > 0, "m must be positive"
        self.p = p
        self.m = m

        # Extension field case
        if m != 1:
            if irr:
                assert (
                    irr.get_max_degree() == m
                ), "irreducible polynom is not for the finite field"

                # Make irreducible polynom's elements are within [0...p-1]
                irr.broadcast_modulo(p)

                self.irr = irr
                # Check irreduciblity
                assert check_irr(
                    FFElement(self, FastPolynom({1: 1}))
                ), "polynom is reducible"
            # else:
            #     irr =
        self.irr = irr

    def __eq__(self, x):
        """
        Check if x is also defining the same field

        Arguments:
            x - GF object

        Usage:
            GF(2, 3, FastPolynom({0:1, 1:1, 3:1}))
            == GF(2, 3, FastPolynom({0:1, 1:1, 3:1}))
        """
        return self.p == x.p and self.m == x.m and self.irr == x.irr

    def __ne__(self, x):
        return not (self == x)

    def __str__(self):
        """
        Transform finite field in universal format:
        `GF(p^m)[X] / F(x)`
        with F(x) is the irreducible polynom

        Usage:
            str(GF(2, 3, FastPolynom({0:1, 1:1, 3:1})))
        """
        p = self.p
        m = self.m
        return "GF({}{}){}".format(
            p,
            "^" + str(m) if m != 1 else "",
            "[X] / " + str(self.irr) if self.irr else "",
        )


class FFElement:
    """
    Galois Field Element/Member class implementation

    Element is defined using polynomial data structure (FastPolynom object).
    We must make sure that the element is indeed within the finite field.
    Or we can try to transform it by fitting it into the defined finite field.

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

    Usage:
        FFElement(GF(2, 3, FastPolynom({0:1, 1:1, 3:1})))
    """

    def __init__(self, ff, container=None):
        """
        Init finite field element
        1. Check if GF object `ff` is defined
        2. Check if container is defined
        3. If container is undefined, create empty polynom
        4. If container is defined, check container polynom's maximum degree
        5. If container is not within the field,
            if prime field - raise PrimeFieldNoFitException
            if extension field - try to fit the container (using modulo)

        Arguments:
            ff - GF object
            container - FastPolynom object (default: None)
                        Set to None to create empty polynom

        Usage:
            FFElement(GF(2, 3, FastPolynom({0:1, 1:1, 3:1})))
        """
        self.ff = ff
        if ff.m:
            if container:
                self.container = container

                self.container.broadcast_modulo(ff.p)

                if self.container.get_max_degree() >= ff.m:
                    if ff.m != 1:
                        self._fit()
                    else:
                        raise PrimeFieldNoFitException()
            else:
                self.container = FastPolynom()

    def __eq__(self, x):
        return self.ff == x.ff and self.container == x.container

    def __ne__(self, x):
        return not (self == x)

    def _div(self, pol1, pol2, p, copy=True):
        """
        *** Do not use ***
        Intermediate division operator for 2 polynom
        Make sure their elements are within [0...p-1]

        Arguments:
            pol1 - FastPolynom object (numerator)
            pol2 - FastPolynom object (divisor)
            p - prime number
            copy - bool, copy `pol1` instead of changing `pol1` (default: True)
        """
        d1 = pol1.get_max_degree()
        d2 = pol2.get_max_degree()
        if d1 < d2:
            return FastPolynom(), pol1
        if (d1 == pol1.container) or (d2 == 0 and not pol2.container):
            if not pol2.container:
                raise FFOperationException("//", "Division by 0 is unallowed")
            if not pol1.container or pol1.container[0] == 1 or pol2.container[0] == 1:
                return pol1, FastPolynom()

        keys_pol2 = pol2.get_keys(rev=True)
        num = pol1.deepcopy() if copy else pol1
        divs = FastPolynom()
        key = pol2[d2]
        while d1 >= d2:
            x = num[d1]
            if x % key == 0:
                fac = x // key
            else:
                fac = x * egcd(self.ff.p, key)[1] % self.ff.p
            divs[d1 - d2] = fac
            if fac:
                for j in keys_pol2:
                    num[d1 - d2 + j] -= (pol2[j] * fac) % p
                    num[d1 - d2 + j] %= p
            d1 = num.get_max_degree()
        return (divs, num)

    def __add__(self, x):
        """
        Addition between 2 FFElement (same as polynomial addition)

        Arguments:
            x - FFElement object

        Usage:
            FFElement(GF(2, 3, FastPolynom({0:1, 1:1, 3:1})))
            + FFElement(GF(2, 3, FastPolynom({0:1, 1:1, 3:1})), {0: 1})
        """
        try:
            assert self.ff == x.ff, "x is not in the same finite field"
            res = FFElement(self.ff)
            for i in range(self.ff.m):
                res.container[i] = (self.container[i] + x.container[i]) % self.ff.p
            return res
        except AttributeError:
            raise FFOperationException("+", "x is not FFElement object?")

    def __sub__(self, x):
        """
        Subtraction between 2 FFElement (same as polynomial subtraction)

        Arguments:
            x - FFElement object

        Usage:
            FFElement(GF(2, 3, FastPolynom({0:1, 1:1, 3:1})))
            - FFElement(GF(2, 3, FastPolynom({0:1, 1:1, 3:1})), {0: 1})
        """
        try:
            assert self.ff == x.ff, "x is not in the same finite field"
            res = FFElement(self.ff)
            for i in range(self.ff.m):
                res.container[i] = (self.container[i] - x.container[i]) % self.ff.p
            return res
        except AttributeError:
            raise FFOperationException("-", "x is not FFElement object?")

    def __mul__(self, x):
        """
        Multiplicate between 2 FFElement (same as polynomial multiplication)
        but, modulo the result with irreducible polynom in finite field
        so the result will be within the field

        Arguments:
            x - FFElement object

        Usage:
            FFElement(GF(2, 3, FastPolynom({0:1, 1:1, 3:1})), {1: 1})
            * FFElement(GF(2, 3, FastPolynom({0:1, 1:1, 3:1})), {0: 1})
        """
        try:
            assert self.ff == x.ff, "x is not in the same finite field"
            deg = self.ff.m
            p = self.ff.p
            res = FastPolynom()
            for i in range(deg):
                for j in range(deg):
                    res[i + j] += (self.container[i] * x.container[j]) % p
                    res[i + j] %= p
            if self.ff.m != 1:
                _, res = self._div(res, self.ff.irr, p)
            return FFElement(self.ff, res)
        except AttributeError:
            raise FFOperationException("*", "x is not FFElement object?")

    def __truediv__(self, x):
        """
        Division between 2 FFElement by using inverse of x
        and multiplication

        Arguments:
            x - FFElement object

        Usage:
            FFElement(GF(2, 3, FastPolynom({0:1, 1:1, 3:1})), {1: 1})
            / FFElement(GF(2, 3, FastPolynom({0:1, 1:1, 3:1})), {0: 1})
        """
        try:
            assert self.ff == x.ff, "x is not in the same finite field"
            return self * x.inverse()
        except FFOperationException:
            raise FFOperationException("/", "Division by 0 is unallowed")
        except AttributeError:
            raise FFOperationException("/", "x is not FFElement object?")

    def __floordiv__(self, x):
        """
        Polynomial division (not finite field division)

        Arguments:
            x - FFElement object

        Usage:
            FFElement(GF(2, 3, FastPolynom({0:1, 1:1, 3:1})), {1: 1})
            // FFElement(GF(2, 3, FastPolynom({0:1, 1:1, 3:1})), {0: 1})
        """
        try:
            assert self.ff == x.ff, "x is not in the same finite field"
            if self.is_zero():
                raise FFOperationException("//", "Division by 0 is unallowed")
            if self.ff.m == 1 or x.is_integer():
                res = FastPolynom()
                for x in self.container.get_keys():
                    res.container[x] = self.container[x] // x.container[0]
            else:
                res, _ = self._div(self.container, x.container, self.ff.p)
            return FFElement(self.ff, res)
        except AttributeError:
            raise FFOperationException("//", "x is not FFElement object?")

    def __mod__(self, x):
        """
        Polynomial modulo

        Arguments:
            x - FFElement object

        Usage:
            FFElement(GF(2, 3, FastPolynom({0:1, 1:1, 3:1})), {1: 1})
            % FFElement(GF(2, 3, FastPolynom({0:1, 1:1, 3:1})), {0: 1})
        """
        try:
            assert self.ff == x.ff, "x is not in the same finite field"
            if self.ff.m == 1 or x.is_integer():
                res = FastPolynom(self.ff)
                for x in self.container.get_keys():
                    res.container[x] = self.container[x] % x.container[0]
            else:
                d = self.container.get_max_degree()
                pre = {}
                a = 1
                _, _x = self._div(FastPolynom({1: 1}), x.container, self.ff.p)
                pre[1] = FFElement(self.ff, _x)
                while a * 2 <= d:
                    _, _x = self._div(
                        (pre[a] * pre[a]).container, x.container, self.ff.p
                    )
                    a *= 2
                    pre[a] = FFElement(self.ff, _x)
                res = FFElement(self.ff)
                for i in self.container.get_keys(rev=True):
                    temp = FFElement(self.ff, FastPolynom({0: self.container[i]}))
                    j = 1
                    while i:
                        if i & 1:
                            temp *= pre[j]
                            _, a = self._div(
                                temp.container, x.container, self.ff.p, copy=False
                            )
                        i >>= 1
                        j *= 2
                    res += temp
            return res
        except AttributeError:
            raise FFOperationException("%", "x is not FFElement object?")

    def pow(self, m):
        """
        Polynomial power

        Arguments:
            m - power (positive integer)

        Usage:
            FFElement(GF(2, 3, FastPolynom({0:1, 1:1, 3:1})), {1: 1}).pow(3)
        """
        assert m > 0, "m must be positive"
        pre = {}
        a = 1
        pre[1] = self
        while a * 2 <= m:
            t = a
            a *= 2
            pre[a] = pre[t] * pre[t]
        res = FFElement.gen_one(self.ff)
        j = 1
        while m:
            if m & 1:
                res *= pre[j]
            m >>= 1
            j *= 2
        return res

    def is_one(self):
        """
        Check if FFElement is one
        """
        return self.container.get_max_degree() == 0 and self.container[0] == 1

    def is_zero(self):
        """
        Check if FFElement is zero
        """
        return self.container.get_max_degree() == 0 and not self.container.container

    def is_integer(self):
        """
        Check if FFElement is integer
        """
        return self.container.get_max_degree() == 0 and self.container[0] != 0

    @staticmethod
    def gen_one(ff):
        """
        Create finite field element = 1

        Arguments:
            ff - GF object

        Usage:
            FFElement.gen_one(GF(2, 3, FastPolynom({0:1, 1:1, 3:1})))
        """
        return FFElement(ff, FastPolynom({0: 1}))

    @staticmethod
    def egcd(a, b):
        """
        Extended euclidean for polynom

        Arguments:
            a - FFElement object (greater)
            b - FFElement object (lower)
        """
        if (a % b).is_zero():
            raise Exception("a & b must be co-prime")
        mem = [FFElement(a.ff), FFElement.gen_one(a.ff)]
        while not b.is_integer():
            t = mem[1]
            mem[1] = mem[0] - t * (a // b)
            mem[0] = t
            t = b
            b = a % b
            a = t
        return mem[1] * b.inverse()

    def inverse(self):
        """
        Inverse the finite field element

        Arguments:

        Usage:
            FFElement(GF(2, 3, FastPolynom({0:1, 1:1, 3:1})), {1: 1}).inverse()
        """
        try:
            if self.is_zero():
                raise FFOperationException("^-1", "0 doesn't have any inverse")
            if self.is_one():
                return FFElement.gen_one(self.ff)
            if self.ff.m == 1 or self.container.get_max_degree() == 0:
                res = FFElement(self.ff)
                res.container[0] = egcd(self.ff.p, self.container[0])[1] % self.ff.p
            else:
                irr = FFElement(self.ff, self.container)
                irr.container = self.ff.irr

                res = FFElement.egcd(irr, self)
            return res
        except AttributeError:
            raise FFOperationException("^-1", "Something went wrong?")

    def _fit(self):
        """
        *** Do not use ***
        Fit the container into the finite field by using modulo
        of irreducible polynom (for extension field only)

        Arguments:

        """
        try:
            x = FFElement(self.ff)
            x.container = self.ff.irr
            self.container = (self % x).container
        except AttributeError:
            raise FFOperationException("fit", "Something went wrong?")

    def __str__(self):
        """
        Transform finite field element in universal format:
        `GF(p^m)[X] / F(x): P(x)`
        with F(x) is the irreducible polynom and P(x) is the element

        Usage:
            str(FFElement(GF(2, 3, FastPolynom({0:1, 1:1, 3:1})), {1: 1}))
        """
        return "{}: {}".format(self.ff, self.container)
