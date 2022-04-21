import numpy as np
from channel import read_dnsin

# read dns.in
infile = read_dnsin('../dns.in')
ny = infile['ny']

# read the file
fromdisk = np.fromfile('instantaneous.bin',dtype=np.float64,sep='').reshape((-1,3,ny+1))
no_files = fromdisk.shape[0]
varhis = -9 * np.ones((no_files,3,ny+3))
varhis[:,:,1:-1] = fromdisk

# write ascii files
np.savetxt('uu_hist.txt', varhis[:,0,:], delimiter='  ')
np.savetxt('vv_hist.txt', varhis[:,1,:], delimiter='  ')
np.savetxt('ww_hist.txt', varhis[:,2,:], delimiter='  ')