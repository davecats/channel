import numpy as np
import tactical as tct
import dns_channel as ch



def read_psd(fdir, **kwargs):

    y_symm = kwargs.get('y_symm', True)
    
    fname = fdir + 'psd.bin' # generate file name

    # fetch nx, ny, nz
    dnsdict = ch.read_dnsin(fdir)
    mesh = ch.mesh(dnsdict)

    # create memmap
    diskacc = np.memmap(fname, mode='r', dtype=np.float64, shape=(6, mesh.ny+1, mesh.nz+1, mesh.nx+1))

    # check file size and allocate memory
    tct.io.size_ram_check(fname)
    all_spectra = np.zeros((6, mesh.ny+1, mesh.nz+1, mesh.nx+1))

    # load
    all_spectra += diskacc

    # average on y (if requested)
    if y_symm:
        all_spectra += diskacc[:,::-1,:,:]
        all_spectra /= 2

    return all_spectra, mesh.kx, mesh.kz, mesh.y[1:-1]