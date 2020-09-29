from __future__ import print_function

class colors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

class context:
    def __init__(self, color):
        self.color = color
    def __enter__(self):
        print( self.color, end="" )
        return self
    def __exit__(self, type, value, traceback):
        print( colors.ENDC, end="" )

verbose=True
        

def iprint(str):
    if verbose:
        print( colors.OKGREEN + "[I] " + str   + colors.ENDC )

def wprint(str):
    if verbose:
        print( colors.WARNING + "[W] "   + str + colors.ENDC )

def eprint(str):
    if verbose:
        print( colors.FAIL + "[E] "  + str  + colors.ENDC)

