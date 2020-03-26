class SpeciesError(Exception):

    def __init__(self, msg=None):
        self.msg = msg

    def __str__(self):
        return str(self.msg)


class KineticsError(Exception):

    def __init__(self, msg=None):
        self.msg = msg


class ConcentrationError(Exception):
    def __init__(self, msg=None):
        self.msg = msg
