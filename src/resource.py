# Mock resource module for Windows compatibility
RUSAGE_SELF = 0

class StructRusage:
    def __init__(self):
        self.ru_maxrss = 0

def getrusage(who):
    return StructRusage()
