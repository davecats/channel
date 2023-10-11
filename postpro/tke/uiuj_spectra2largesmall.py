from argparse import ArgumentParser
import channel as ch
import numpy as np

parser = ArgumentParser(
    prog = 'uiuj_spectra2largesmall',
    description='Converts the output of uiuj_spectra into the one of uiuj_largesmall.'
    )
parser.add_argument('target_dir', metavar='sim_folder/', type=str, nargs=1, help='the directory containing the simulation data')
settings = parser.parse_args()
base_dir = settings.target_dir[0] + '/'; base_dir = base_dir.replace('//','/')

psdfile = base_dir + '/uiuj_spectra/uiuj_spectra.bin'
dnsin = ch.read_dnsin(base_dir + '/dns.in')
beta0 = dnsin['beta0']
m = ch.mesh(dnsin); y = m.y
nx = dnsin['nx']; ny = dnsin['ny']; nz = dnsin['nz']

# read uiuj_spectra
arr_len = (ny+3)*(2*nz+1)*10*6
offs = (ny+3)*7*8
sp_data = np.fromfile(psdfile, dtype=np.float64, count=arr_len, offset=offs).reshape((ny+3, 2*nz+1,10,6))

# read largesmall_settings
with open('largesmall_settings.in','r') as f:
    kz_thr = float(f.readline().split('!')[0])

# get mask for scale separation
kz = beta0 * np.linspace(-nz,nz,2*nz + 1)
mask_true_is_large_1d = np.abs(kz) <= kz_thr
mask_true_is_large_1d = mask_true_is_large_1d.reshape(1, 2*nz+1,1,1)
all_true = np.empty((ny+3, 2*nz+1,10,6),dtype=bool)
all_true[:] = True
mask_true_is_large = np.logical_and(all_true, mask_true_is_large_1d)

# masked arrays
small_sp = np.ma.masked_array(sp_data, mask_true_is_large)
large_sp = np.ma.masked_array(sp_data, np.logical_not(mask_true_is_large))

# reduce (sum over all kz)
small_reduced = np.sum(small_sp, axis=1)
large_reduced = np.sum(large_sp, axis=1)

# allocate and fill uiuj_largesmall array
uiuj_ls = np.zeros((2,ny+3,11,6), dtype = np.float64)
for ils, src in zip([0,1],[large_reduced,small_reduced]):
    uiuj_ls[ils,:,0:3,:] = src[:,0:3,:]
    uiuj_ls[ils,:,5:,:] = src[:,4:,:]

# retrieve info for binary header
header_data = np.empty((4),dtype=np.intc)
with open(psdfile.replace('.bin','.nfo'),'r') as f:
    linez = f.readlines()
for line in linez:
    if 'nmin' in line:
        header_data[0] = int(line.split('nmin')[-1])
    elif 'nmax' in line:
        header_data[1] = int(line.split('nmax')[-1])
    elif 'dn' in line:
        header_data[2] = int(line.split('dn')[-1])
    elif 'nftot' in line:
        header_data[3] = int(line.split('nftot')[-1])

# read mean data
arr_len = (m.ny+3)*7
mean_data = np.fromfile(psdfile, dtype=np.float64, count=arr_len, offset=0).reshape((ny+3, 7))

# prepare to write
# size_of_header = header_data.nbytes; size_of_mean = mean_data.nbytes; size_of_uiuj = uiuj_ls.nbytes
towrite = header_data.tobytes() + mean_data.tobytes() + uiuj_ls.tobytes()

# write to disk
with open('profiles/uiuj_largesmall.bin', 'wb') as of:
    of.write(towrite)