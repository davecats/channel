from channel import field_as_memmap
import numpy as np
from argparse import ArgumentParser

# parse input
parser = ArgumentParser(description='Read one or more fields and check for NaNs. If a NaN is present, a message is printed on screen. PLEASE NOTICE: this script TREATS INF AS IF IT WERE NAN!')
parser.add_argument('fields', metavar='field.out', type=str, nargs='*', help='One or more fields to be checked for NaNs.')
settings = parser.parse_args()

# do the checking
for field in settings.fields:
    fld = field_as_memmap(field)
    if np.logical_not(np.isfinite(fld)).any():
        print('Field ' + field + ' contains some non numerical values.')