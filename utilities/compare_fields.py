from channel import read_field
import numpy as np
import sys
import matplotlib
matplotlib.use('qt4Agg')
from matplotlib.pyplot import plot, show

argz = sys.argv[1:]

orig = read_field(argz[0])
curr = read_field(argz[1])

intez_first = np.fromfile(argz[0], count=3, dtype=np.int32)
realz_first = np.fromfile(argz[0], count=7, dtype=np.float64, offset=12)

intez_secnd = np.fromfile(argz[1], count=3, dtype=np.int32)
realz_secnd = np.fromfile(argz[1], count=7, dtype=np.float64, offset=12)

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

def abs_cwise(arr):
    arr = arr.ravel()
    result = np.zeros((len(arr)))
    for ii,el in enumerate(arr):
        result[ii] = max(abs(el.real), abs(el.imag))
    return result

diff = abs_cwise(orig-curr)

plot(diff)

maxdiff = np.amax(diff)
avgdiff = np.mean(diff)

if not maxdiff == 0:
    print('Average error:',avgdiff)
    print('Maximum error:',maxdiff)
else:
    print('No difference detected in velocity field.')

show()
