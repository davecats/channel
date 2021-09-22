import numpy as np
from argparse import ArgumentParser
from channel import read_dnsin
from progressbar import progressbar

indict = read_dnsin('dns.in')
nx = indict['nx']
ny = indict['ny']
nz = indict['nz']

# parse input arguments
argpar = ArgumentParser(description='Shift a velocity field in the x direction by half an extimate of the centerline velocity; designed for Poiseuille (symmetric) flows.')
argpar.add_argument('--centerline', '-c', help='Extimate of centerline velocity.', nargs=1, type=float)
argpar.add_argument('--reverse', '-r', help='Reverse a previous shift, so that wall velocity becomes zero.', action='store_true')
argpar.add_argument('file_list', metavar='file.ext', type=str, nargs='+',
                    help='Files on which to apply shift.')
settings = argpar.parse_args()

# basically parse arguments
if settings.reverse:
    u_shift = 0.0
    trgt_shift = -0.0
else:
    uc_ext = settings.centerline[0]
    print('Extimate of centerline velocity: ', uc_ext)
    u_shift = uc_ext/2
    trgt_shift = u_shift # store value of shift for later reference

# trgt_shift is shift with respect to zero-shift situation (still walls)
# u_shift is shift relative to current shift

# retrieve previous shift
prev_shift = 0
if indict['u0'] == indict['un']:
    prev_shift = -indict['u0']
else:
    raise Exception('mismatch between u0 and un in dns.in.')
# calculate differential shift
u_shift -= prev_shift

print('Shifting by ', u_shift)

# calculate header
integer_size = 4 # size (in bytes) of an integer
doubles_to_skip = 7 # doubles saved before the velocity field (in addition to the first three integers)
bytes_in_double = 8 # number of bytes in a double
bytes_to_skip = doubles_to_skip*bytes_in_double + 3*integer_size # number of bytes to skip (which is, before array)

for ff in progressbar(settings.file_list):
    # read array
    v_field = np.memmap(ff, dtype=np.complex128, mode='readwrite', offset=bytes_to_skip, shape=(3, nx+1, 2*nz+1, ny+3))
    for iy in range(ny+3):
        v_field[0,0,nz,iy] -= (u_shift)

# change dns.in
with open('dns.in', 'r') as infile:
    all_lines = infile.readlines() # read lines
    for (ii, line) in enumerate(all_lines):
        if 'u0' in line: # search for line of u0, uN
            all_lines[ii] = str(-trgt_shift) + "\t" + str(-trgt_shift) + "\t" + '! u0, uN' + '\n' # change line
with open('dns.in', 'w') as infile:
    infile.writelines(all_lines)