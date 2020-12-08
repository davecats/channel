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

# read headers
intez_first = np.fromfile(field1, count=3, dtype=np.int32); realz_first = np.fromfile(field1, count=7, dtype=np.float64, offset=12)
intez_secnd = np.fromfile(field2, count=3, dtype=np.int32); realz_secnd = np.fromfile(field2, count=7, dtype=np.float64, offset=12)

# compare headers
if np.array_equal(intez_first, intez_secnd) and np.array_equal(realz_first, realz_secnd):
    print('No difference detected in header.')
else:
    print('There are difference in the headers.')
    print('Domain size:')
    print('first ->', intez_first)
    print('second ->', intez_secnd)
    print('Other floats:')
    print('first ->', realz_first)
    print('second ->', realz_secnd)

# define function to calculate difference between fields
def abs_cwise(arr):
    #arr = arr.ravel()
    #result = np.zeros((len(arr)))
    #for ii,el in enumerate(arr):
    #    result[ii] = max(abs(el.real), abs(el.imag))
    return np.maximum(abs(arr.real), abs(arr.imag)).ravel()

# read fields
orig = field_as_memmap(field1)
curr = field_as_memmap(field2)

# get difference
diff = abs_cwise(orig-curr)

# compare fields
maxdiff = np.amax(diff)
avgdiff = np.mean(diff)
if not maxdiff == 0:
    print('Average error:',avgdiff)
    print('Maximum error:',maxdiff)
else:
    print('No difference detected in velocity field.')

# if necessary, plot
if settings.plot:
    fig, ax = subplots(num=field1+' vs '+field2)
    ax.plot(diff)
show()
