# Galois Field

![galois-field](https://github.com/haverzard/GaloisField/workflows/Galois%20Field/badge.svg)

## Description
Galois Field (GF) is a field contains a finite number of elements. There are 2 types of Galois Field:
1. Prime Field (m = 1)
2. Extension Field (m != 1)

In prime field, elements are integer within [0, p-1] range. Prime field have a prime `p` that limits our value so it will always be within the field.

In extension field, elements can be polynomials with maximum degree of (m-1). Extension field have a prime `p` and prime polynomial (irreducible) that limits our polynomial and its values so it will always be within the field.

## Install
`pip install GaloisField`

## Requirements
- Python 3

## Usage
```
from galois_field.GF import GF, FFElement
from galois_field.fast_polynom import FastPolynom

irr_poly = FastPolynom({0: 1, 1: 1, 3: 1})
ff = GF(2, 3, (irr_poly, [1]))
e1 = FFElement(ff, FastPolynom({0: 1, 1: 1}))
e2 = FFElement(ff, FastPolynom({0: 1, 1: 1}))
e_res = e1 * e2
assert (e2 == e_res / e1)
```