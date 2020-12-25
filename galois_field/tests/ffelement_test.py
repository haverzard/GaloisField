from unittest import TestCase
from galois_field.GF import GF, FFElement
from galois_field.fast_polynom import FastPolynom
from galois_field.exceptions import PrimeFieldNoFitException


class TestFFElement(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestFFElement, self).__init__(*args, **kwargs)
        self.ff = GF(2, 3, (FastPolynom({0: 1, 1: 1, 3: 1}), [1]))
        self.ff2 = GF(2, 3, (FastPolynom({0: 1, 2: 1, 3: 1}), [1]))

    def test_create(self):
        fe = FFElement(self.ff, FastPolynom({1: 1}))
        self.assertTrue(isinstance(fe, FFElement))
        self.assertEqual(str(fe), "GF(2^3)[X] / x^3 + x^1 + 1: x^1")

    def test_equality(self):
        fe1 = FFElement(self.ff, FastPolynom({1: 1}))
        fe2 = FFElement(self.ff, FastPolynom({1: 1}))
        self.assertEqual(fe1, fe2)

        fe1 = FFElement(self.ff, FastPolynom({1: 1}))
        fe2 = FFElement(self.ff, FastPolynom({1: 2}))
        self.assertNotEqual(fe1, fe2)

        fe1 = FFElement(self.ff, FastPolynom({1: 1}))
        fe2 = FFElement(self.ff2, FastPolynom({1: 1}))
        self.assertNotEqual(fe1, fe2)


class TestFFException(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestFFException, self).__init__(*args, **kwargs)
        self.ff = GF(101)

    def test_no_prime_field_fitting(self):
        self.assertRaises(
            PrimeFieldNoFitException, FFElement, self.ff, FastPolynom({1: 1})
        )
