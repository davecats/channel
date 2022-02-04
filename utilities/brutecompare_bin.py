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
one = np.fromfile(field1, count=-1, dtype=np.float64, offset=0)
two = np.fromfile(field2, count=-1, dtype=np.float64, offset=0)

def diff_analysis(refarr,arr):

    # get difference
    diff = abs(refarr - arr)

    # prepare field to normalise
    nrmlstr = abs(refarr)
    nozero = np.nonzero(nrmlstr)
    if len(nozero[0])==0:
        print('')
        print('Reference array is zero everywhere.')
        print('')
        return
    reldiff = diff[nozero]/nrmlstr[nozero]
    arr_1 = reldiff[np.greater(reldiff,0.01)]
    arr_10 = reldiff[np.greater(reldiff,0.1)]

    # compare fields
    maxdiff = np.amax(diff)
    if not maxdiff == 0:
        print()
        print('Absolute difference')
        print('Average difference:', np.mean(diff) )
        print('Median difference:', np.median(diff) )
        print('Maximum difference:',maxdiff)
        print()
        print('Relative difference (excluding zeroes)')
        print('Reference file is the first one: ', field1)
        print('Number of points used:', len(reldiff))
        print('Average rel difference:', np.mean(reldiff) )
        print('Median relative difference:', np.median(reldiff) )
        print('Maximum rel difference:', np.amax(reldiff) )
        print('Number of pts with more than 1% difference:',len(arr_1),'out of',len(reldiff),'-',len(arr_1)/len(reldiff)*100,'% of total')
        print('Number of pts with more than 10% difference:',len(arr_10),'out of',len(reldiff),'-',len(arr_10)/len(reldiff)*100,'% of total')
        print()
    else:
        print()
        print('No difference detected.')
        print()

diff_analysis(one,two)