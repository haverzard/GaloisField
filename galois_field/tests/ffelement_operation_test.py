from unittest import TestCase
from galois_field.GF import GF, FFElement
from galois_field.fast_polynom import FastPolynom
from galois_field.exceptions import FFOperationException
import operator


class TestFFOperation(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestFFOperation, self).__init__(*args, **kwargs)
        self.ff = GF(2, 3, (FastPolynom({0: 1, 1: 1, 3: 1}), [1]))
        self.ff2 = GF(101)
        self.ff3 = GF(101, 3, (FastPolynom({0: 10, 1: 6, 3: 1}), [2, 5]))

    def test_addition(self):
        fe_res = FFElement(self.ff, FastPolynom({1: 1}))
        fe1 = FFElement(self.ff, FastPolynom({0: 1, 1: 1, 2: 1}))
        fe2 = FFElement(self.ff, FastPolynom({0: 1, 2: 1}))
        self.assertEqual(fe1 + fe2, fe_res)

        fe_res = FFElement(self.ff2, FastPolynom({0: 30}))
        fe1 = FFElement(self.ff2, FastPolynom({0: 100}))
        fe2 = FFElement(self.ff2, FastPolynom({0: 31}))
        self.assertEqual(fe1 + fe2, fe_res)

    def test_subtraction(self):
        fe_res = FFElement(self.ff, FastPolynom({1: 1}))
        fe1 = FFElement(self.ff, FastPolynom({0: 1, 1: 1, 2: 1}))
        fe2 = FFElement(self.ff, FastPolynom({0: 1, 2: 1}))
        self.assertEqual(fe1 - fe2, fe_res)

        fe_res = FFElement(self.ff2, FastPolynom({0: 80}))
        fe1 = FFElement(self.ff2, FastPolynom({0: 10}))
        fe2 = FFElement(self.ff2, FastPolynom({0: 31}))
        self.assertEqual(fe1 - fe2, fe_res)

    def test_multiplication(self):
        fe_res = FFElement(self.ff, FastPolynom({0: 1, 2: 1}))
        fe1 = FFElement(self.ff, FastPolynom({0: 1, 1: 1, 2: 1}))
        fe2 = FFElement(self.ff, FastPolynom({1: 1}))
        self.assertEqual(fe1 * fe2, fe_res)

        fe_res = FFElement(self.ff2, FastPolynom({0: 7}))
        fe1 = FFElement(self.ff2, FastPolynom({0: 10}))
        fe2 = FFElement(self.ff2, FastPolynom({0: 31}))
        self.assertEqual(fe1 * fe2, fe_res)

    def test_inverse(self):
        fe = FFElement(self.ff, FastPolynom({0: 1, 1: 1, 2: 1}))
        fe_inv = fe.inverse()
        self.assertEqual(fe * fe_inv, FFElement.gen_one(self.ff))

        fe = FFElement(self.ff2, FastPolynom({0: 10}))
        fe_inv = fe.inverse()
        self.assertEqual(fe * fe_inv, FFElement.gen_one(self.ff2))

        fe = FFElement(self.ff, FastPolynom({0: 1}))
        fe_inv = fe.inverse()
        self.assertEqual(fe_inv, FFElement.gen_one(self.ff))

        fe = FFElement(self.ff3, FastPolynom({0: 14, 1: 10}))
        fe_inv = fe.inverse()
        self.assertEqual(fe * fe_inv, FFElement.gen_one(self.ff2))

    def test_division(self):
        fe1 = FFElement(self.ff, FastPolynom({0: 1, 1: 1, 2: 1}))
        fe2 = FFElement(self.ff, FastPolynom({1: 1}))
        fe_res = FFElement(self.ff, FastPolynom({0: 1, 2: 1}))
        self.assertEqual(fe1 * fe2, fe_res)
        self.assertEqual(fe_res / fe1, fe2)

        fe1 = FFElement(self.ff2, FastPolynom({0: 31}))
        fe2 = FFElement(self.ff2, FastPolynom({0: 10}))
        fe_res = FFElement(self.ff2, FastPolynom({0: 7}))
        self.assertEqual(fe1 * fe2, fe_res)
        self.assertEqual(fe_res / fe1, fe2)


class TestFFException(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestFFException, self).__init__(*args, **kwargs)
        self.ff = GF(2, 3, (FastPolynom({0: 1, 1: 1, 3: 1}), [1]))
        self.ff2 = GF(2, 3, (FastPolynom({0: 1, 2: 1, 3: 1}), [1]))
        self.ff3 = GF(101)

    def test_not_in_the_same_field(self):
        fe1 = FFElement(self.ff, FastPolynom({0: 1, 1: 1, 2: 1}))
        fe2 = FFElement(self.ff2, FastPolynom({1: 1}))
        self.assertRaises(AssertionError, operator.add, fe1, fe2)
        self.assertRaises(AssertionError, operator.sub, fe1, fe2)
        self.assertRaises(AssertionError, operator.mul, fe1, fe2)
        self.assertRaises(AssertionError, operator.truediv, fe1, fe2)

    def test_not_between_ffelement(self):
        fe1 = FFElement(self.ff, FastPolynom({0: 1, 1: 1, 2: 1}))
        fe2 = FastPolynom({1: 1})
        self.assertRaises(FFOperationException, operator.add, fe1, fe2)
        self.assertRaises(FFOperationException, operator.sub, fe1, fe2)
        self.assertRaises(FFOperationException, operator.mul, fe1, fe2)
        self.assertRaises(FFOperationException, operator.truediv, fe1, fe2)

    def test_division_by_zero(self):
        fe1 = FFElement(self.ff, FastPolynom({0: 1, 1: 1, 2: 1}))
        fe2 = FFElement(self.ff)
        self.assertRaises(FFOperationException, operator.truediv, fe1, fe2)
