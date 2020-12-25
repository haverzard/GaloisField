from unittest import TestCase
from galois_field.GF import GF
from galois_field.fast_polynom import FastPolynom


class TestGF(TestCase):
    def test_extension_field(self):
        gf = GF(2, 3, (FastPolynom({0: 1, 1: 1, 3: 1}), [1]))
        self.assertTrue(isinstance(gf, GF))
        self.assertEqual(gf.p, 2)
        self.assertEqual(gf.m, 3)
        self.assertEqual(str(gf), "GF(2^3)[X] / x^3 + x^1 + 1")

    def test_prime_field(self):
        gf = GF(2)
        self.assertTrue(isinstance(gf, GF))
        self.assertEqual(gf.p, 2)
        self.assertEqual(gf.m, 1)
        self.assertEqual(gf.irr[0], None)
        self.assertEqual(str(gf), "GF(2)")

    def test_equality(self):
        gf1 = GF(2)
        gf2 = GF(2)
        self.assertEqual(gf1, gf2)

        gf1 = GF(2)
        gf2 = GF(3)
        self.assertNotEqual(gf1, gf2)

        gf1 = GF(2, 3, (FastPolynom({0: 1, 1: 1, 3: 1}), [1]))
        gf2 = GF(2, 3, (FastPolynom({0: 1, 1: 1, 3: 1}), [1]))
        self.assertEqual(gf1, gf2)

        gf1 = GF(2, 3, (FastPolynom({0: 1, 1: 1, 3: 1}), [1]))
        gf2 = GF(2, 3, (FastPolynom({0: 1, 2: 1, 3: 1}), [1]))
        self.assertNotEqual(gf1, gf2)


class TestGFException(TestCase):
    def test_not_prime(self):
        self.assertRaises(AssertionError, GF, 6, 3)

    def test_not_positive_m(self):
        self.assertRaises(AssertionError, GF, 2, -3)
        self.assertRaises(AssertionError, GF, 2, 0)

    def test_not_irreducible(self):
        self.assertRaises(AssertionError, GF, 2, 3, (FastPolynom({0: 1, 3: 1}), [1]))
