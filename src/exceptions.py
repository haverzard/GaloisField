class FFOperationException(Exception):
    def __init__(self, op="", msg=""):
        self.operator = op
        self.message = msg


class FPNegativeDegreeNotAllowed(Exception):
    def __init__(self):
        self.message = "Negative degree polynom is not allowed"
