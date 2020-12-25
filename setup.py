import setuptools

long_description = """
Galois Field (GF) is a field contains a finite number of elements. There are 2 types of Galois Field:
1. Prime Field (m = 1)
2. Extension Field (m != 1)

In prime field, elements are integer within [0, p-1] range. Prime field have a prime `p` that limits our value so it will always be within the field.

In extension field, elements can be polynomials with maximum degree of (m-1). Extension field have a prime `p` and prime polynomial (irreducible) that limits our polynomial and its values so it will always be within the field.
"""

setuptools.setup(
    name="GaloisField",
    version="0.1.0",
    author="haverzard",
    author_email="yonatanviody@gmail.com",
    description="Galois Finite Field Implementation in Python 3",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/haverzard/GaloisField",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.5",
)
