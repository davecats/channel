from channel import read_dnsin
import sys
import os
from math import floor

argz = sys.argv[1:]

# generate path to exec and verify its existance
exec_dir = sys.path[0] + '/print_time' # path where script resides + exec name
if not os.path.isfile(exec_dir):
    raise(Exception('missing executable file "print_time". Please run "make" in the utilities folder of your "channel".'))

# get time of selected fields
stream = os.popen(exec_dir + ' ' + ' '.join(argz))
read_timez = stream.read().split('\n')[0]
read_timez = read_timez.split(' ')
read_timez = list(filter(lambda x: x!='', read_timez))

# read dtfield from dns.in
indns = read_dnsin('dns.in')
dtfield = float(indns['dt_field'])

# THIS PRINTS ONLY FILES WITH MISMATCH BETWEEN TIME AND FILE NUMBER
ii = 0
for ff in argz:
    fno = int(''.join(c for c in ff if c.isdigit())) # file number
    rndt = float(read_timez[ii])
    if not fno == int(round(rndt/dtfield)):
        print(ff, "-->", rndt)
    ii += 1