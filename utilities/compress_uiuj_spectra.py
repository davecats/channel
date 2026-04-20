from argparse import ArgumentParser
import channel as ch
import numpy as np
import os
import gc

parser = ArgumentParser(prog='compress_uiuj_spectra', description='Compresses uiuj_spectra for online distribution.')
parser.add_argument('path_to_uiuj_spectra', type=str, nargs=1)
stngs = parser.parse_args()

bin_file = stngs.path_to_uiuj_spectra[0]
dnsin_file = bin_file.replace('uiuj_spectra.bin','../dns.in')

m = ch.mesh(ch.read_dnsin(dnsin_file))

new_size =  ((m.ny+3)*7) + ((m.ny+3)*((2*m.nz)+1)*10*6) + ((m.ny+3)*10*6) 
old_size =  ((m.ny+3)*7) + ((m.ny+3)*((2*m.nz)+1)*10*6) + ((m.ny+3)*10*6) + (4*((2*m.nz)+1)*(m.nz+1)*(m.ny+3))

old_size_bytes = os.stat(bin_file).st_size
print(f'Old size in bytes: {old_size_bytes}')
print(f'Expected size in bytes: {old_size*8}')
if not (old_size_bytes == (old_size*8)):
    raise Exception('input file is unexpectedly large.')

old_file = np.memmap(bin_file,  dtype=np.float64, mode='r', offset=0, shape=(old_size))
new_file = np.memmap(bin_file.replace('uiuj_spectra.bin','compressed_uiuj_spectra.bin'), dtype=np.float64, mode='w+', shape=(new_size))

new_file[:] = old_file[:new_size]
