# NOTE (2026-10-19 audit): historical legacy script; not part of the current RANO project. Windows mock of the Unix resource module: reports ru_maxrss = 0 by construction (not measured memory). On Linux/macOS it would shadow the stdlib resource module if src/ is on sys.path.
# Mock resource module for Windows compatibility
RUSAGE_SELF = 0

class StructRusage:
    def __init__(self):
        self.ru_maxrss = 0

def getrusage(who):
    return StructRusage()
