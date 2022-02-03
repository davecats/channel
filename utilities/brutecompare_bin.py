from channel import field_as_memmap
import numpy as np
import matplotlib
matplotlib.use('qt4Agg')
from matplotlib.pyplot import subplots, show
from argparse import ArgumentParser

# parse input
parser = ArgumentParser(description='Compare two fields, both in terms of their header and their velocity field.')
parser.add_argument('fields', metavar=('field_1', 'field_2'), type=str, nargs=2, help='Two fields to be compared')
parser.add_argument('--plot', '-p', action='store_true')
settings = parser.parse_args()
field1 = settings.fields[0]
field2 = settings.fields[1]

# read whole binary file as an array of realss
realz_first = np.fromfile(field1, count=-1, dtype=np.float64, offset=0)
realz_secnd = np.fromfile(field2, count=-1, dtype=np.float64, offset=0)

# get difference
diff = abs(realz_first - realz_secnd)

# prepare field to normalise
nrmlstr = abs(realz_first)
nozero = np.nonzero(nrmlstr)

# compare fields
maxdiff = np.amax(diff)
if not maxdiff == 0:
    print()
    print('Absolute difference')
    print('Average difference:', np.mean(diff) )
    print('Median difference:', np.median(diff) )
    print('Maximum difference:',maxdiff)
    print()
    # TODO: remove zeroes here!
    print('Relative difference (excluding zeroes)')
    print('Reference file is the first one: ', field1)
    print('Average rel difference:', np.mean(diff[nozero]/nrmlstr[nozero]) )
    print('Median relative difference:', np.median(diff[nozero]/nrmlstr[nozero]) )
    print('Maximum rel difference:', np.amax(diff[nozero]/nrmlstr[nozero]) )
    print()
else:
    print('No difference detected in velocity field.')
