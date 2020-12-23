class FFOperationException(Exception):
    def __init__(self, op="", msg=""):
        self.operator = op
        self.message = msg
