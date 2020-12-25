import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

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
