#!/usr/bin/python3


class FFOperationException(Exception):
    """
    Finite Field Element Operation Exception
    """

    def __init__(self, op="", msg=""):
        self.operator = op
        self.message = msg


class PrimeFieldNoFitException(Exception):
    """
    No fitting available in prime field
    """

    def __init__(self):
        self.message = (
            "There is no fitting in prime field. Use the correct polynom's degree 1"
        )


class FPNegativeDegreeNotAllowed(Exception):
    """
    Negative degree is unallowed
    """

    def __init__(self):
        self.message = "Negative degree polynom is not allowed"
