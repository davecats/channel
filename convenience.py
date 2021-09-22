import numpy
from struct import unpack
from math import pi, floor



class mesh():
# Extracts mesh and properties. All mesh properties here calculated are meant
# as in the standard scaling units (which is, the ones used for the simulation).
# Syntax:
# mesh_obj = mesh(dns_dictionary)
# where dns_dictionary can be extracted from dns.in with read_dnsin()

    def __init__(self, dnsdata):

        self.nx = dnsdata['nx']
        self.ny = dnsdata['ny']
        self.nz = dnsdata['nz']
        self.alfa0 = dnsdata['alfa0']
        self.beta0 = dnsdata['beta0']
        self.a = dnsdata['a']

        # find nxd,nzd
        self.nxd=3*self.nx//2-1
        self.nzd=3*self.nz-1

        # set grid
        self.y = numpy.tanh(self.a*(2*numpy.arange(-1,self.ny+2)/self.ny-1))/numpy.tanh(self.a)+1
        self.kx = self.alfa0*numpy.arange(0,self.nx+1)
        self.kz = self.beta0*numpy.arange(0,self.nz+1)

        # calculate domain size
        self.lx_pih = 2 / self.alfa0
        self.lz_pih = 2 / self.beta0

        # calculate resolution
        self.dx = self.lx_pih/(2*self.nx) * pi
        self.dz = self.lz_pih/(2*self.nz) * pi
        self.dyw = abs(self.y[1] - self.y[2]) # watch out for ghost cell!
        idxy = floor((self.ny + 3) / 2)
        self.dyc = abs(self.y[idxy] - self.y[idxy+1])




def read_sim_prop(bin_file, dnsin_file, **kwargs):

    integer_size = kwargs.get('int_size', 4) # size (in bytes) of an integer

    # read number of points from dns.in for checking
    check_np = get_dim_dnsin(dnsin_file)

    # read dimensions from binary file; read other parameters as well
    np = [0 for jj in range(3)]
    with open(bin_file, 'rb') as f:
        for ii in range(len(np)):
            byte = f.read(integer_size)
            np[ii] = int.from_bytes(byte, byteorder = 'little')
        alfa0 = unpack('d', f.read(8))[0]
        beta0 = unpack('d', f.read(8))[0]
        _ = unpack('d', f.read(8))[0] # ni
        a = unpack('d', f.read(8))[0]
        ymin = unpack('d', f.read(8))[0]
        ymax = unpack('d', f.read(8))[0]
        _ = unpack('d', f.read(8)) # factor

    # check correspondence between dns.in and binary file
    if not check_np == np:
        print('Error: "{}" and "{}" indicate different values for the number of points (respectively, {} and {}). Please check that files are correct, otherwise check that the correct length (in bytes) of integer values is being used.'.format(bin_file, dnsin_file, str(np), str(check_np)))

    # calculate y
    fun_y = lambda ii : ymin + 0.5*(ymax-ymin)*(numpy.tanh(a*(2*ii/np[1]-1))/numpy.tanh(a)+0.5*(ymax-ymin))
    y = [fun_y(iy) for iy in range(-1, np[1]+2)]
    y = numpy.asarray(y)

    # arrays of kx, kz
    kx = [jj*alfa0 for jj in range(np[0]+1)]
    kx = numpy.asarray(kx)
    kz = [jj*beta0 for jj in range(-np[2], np[2]+1)]
    kz = numpy.asarray(kz)

    # get rid of ghost cells while returning
    return np, y[1:-1], kx, kz



def read_time(bin_file, **kwargs):
    integer_size = kwargs.get('int_size', 4) # size (in bytes) of an integer
    bytes_in_double = 8
    bytes_to_skip = bytes_in_double*6 + integer_size*3 # skip 3 integers and 6 doubles
    t = numpy.fromfile(bin_file, dtype=numpy.double, count=1, sep='', offset=bytes_to_skip)
    return t[0]



def read_field(bin_file, **kwargs):
# Load a velocity field into memory; ghost cells are omitted.
# usage: v_field, nx, ny, nz = read_field( 'Dati.cart.whatever.out' [, integer_size=4] )
# Arguments in square brackets are optional.

    integer_size = kwargs.get('int_size', 4) # size (in bytes) of an integer

    # read dimensions from binary file
    np = [0 for jj in range(3)]
    with open(bin_file, 'rb') as f:
        for ii in range(len(np)):
            byte = f.read(integer_size)
            np[ii] = int.from_bytes(byte, byteorder = 'little')

    doubles_to_skip = 7 # doubles saved before the velocity field (in addition to the first three integers)
    bytes_in_double = 8 # number of bytes in a double
    bytes_to_skip = doubles_to_skip*bytes_in_double + 3*integer_size # number of bytes to skip (which is, before array)

    # calculate expected array size
    array_size = 3 * (np[0] + 1) * (2*np[2] + 1) * (np[1]+3) # components * (nx+1) * (2nz+1) * (ny+3)

    # read array
    v_field = numpy.fromfile(bin_file, dtype=numpy.complex128, count=-1, sep='', offset=bytes_to_skip)
    
    # check array size
    if not len(v_field) == array_size:
        print('Error: unespected length of field read from file ({} elements, while {} were expected).'.format(len(v_field), array_size))
    
    # reshape array (keep in mind that fortran stores data in reversed order!)
    v_field = numpy.reshape(v_field, (3, np[0]+1, 2*np[2]+1, np[1]+3))

    # get rid of ghost cells while returning
    return v_field[:,:,:,1:-1]



def field_as_memmap(bin_file, **kwargs):
# Accesses a field on disk as a memmap, meaning it is not loaded to memory.
# FOR SAFETY REASONS, ACCESS IS READ ONLY BY DEFAULT.
# THE MEMMAP INCLUDES GHOST CELLS; this prevents unnecessary slicing.
# Usage:    v_field = field_as_memmap( 'Dati.cart.whatever.out' [, integer_size=4] )
# Arguments in square brackets are optional.

    np = numpy.fromfile(bin_file, count=3, dtype=numpy.int32) # read files from disk

    integer_size = kwargs.get('int_size', 4) # size (in bytes) of an integer

    doubles_to_skip = 7 # doubles saved before the velocity field (in addition to the first three integers)
    bytes_in_double = 8 # number of bytes in a double
    bytes_to_skip = doubles_to_skip*bytes_in_double + 3*integer_size # number of bytes to skip (which is, before array)

    # read array
    v_field = numpy.memmap(bin_file, dtype=numpy.complex128, mode='r', offset=bytes_to_skip, shape=(3, np[0]+1, 2*np[2]+1, np[1]+3))
    
    return v_field


def get_dim_dnsin(dnsin_file):
    check_np = [0 for jj in range(3)]
    dnsin = read_dnsin(dnsin_file)
    check_np[0] = dnsin['nx']
    check_np[1] = dnsin['ny']
    check_np[2] = dnsin['nz']
    return check_np


def read_dnsin(dnsin_file):

    dns_dict = {}

    with open(dnsin_file) as f:
        # first line: nx, ny, nz
        rawline = rawline_dnsin(f)
        dns_dict['nx'] = int(rawline[0])
        dns_dict['ny'] = int(rawline[1])
        dns_dict['nz'] = int(rawline[2])
        # second line: alfa0 beta0
        rawline = rawline_dnsin(f)
        dns_dict['alfa0'] = float(rawline[0])
        dns_dict['beta0'] = float(rawline[1])
        # third line: ni
        rawline = rawline_dnsin(f)
        dns_dict['re'] = float(rawline[0])
        # fourth line: a, ymin, ymax
        rawline = rawline_dnsin(f)
        dns_dict['a'] = float(rawline[0])
        dns_dict['ymin'] = float(rawline[1])
        dns_dict['ymax'] = float(rawline[2])
        # fifth line: CPI, CPItype, gamma
        rawline = rawline_dnsin(f)
        if rawline[0].upper() == ".TRUE.":
            dns_dict['cpi_flag'] = True
        else:
            dns_dict['cpi_flag'] = False
        dns_dict['cpi_type'] = int(rawline[1])
        dns_dict['gamma'] = float(rawline[2])
        # sixth line: meanpx, meanpz
        rawline = rawline_dnsin(f)
        dns_dict['meanpx'] = float(rawline[0])
        dns_dict['meanpz'] = float(rawline[1])
        # seventh line: meanflowx meanflowz
        rawline = rawline_dnsin(f)
        dns_dict['meanflowx'] = float(rawline[0])
        dns_dict['meanflowz'] = float(rawline[1])
        # eigth line: u0 uN
        rawline = rawline_dnsin(f)
        dns_dict['u0'] = float(rawline[0])
        dns_dict['un'] = float(rawline[1])
        # nineth line: deltat, cflmax, time
        rawline = rawline_dnsin(f)
        dns_dict['deltat'] = float(rawline[0])
        dns_dict['cflmax'] = float(rawline[1])
        dns_dict['time'] = float(rawline[2])
        # tenth line: dt_field, dt_save, t_max, time_from_restart
        rawline = rawline_dnsin(f)
        dns_dict['dt_field'] = float(rawline[0])
        dns_dict['dt_save'] = float(rawline[1])
        dns_dict['t_max'] = float(rawline[2])
        if rawline[3].upper() == ".TRUE.":
            dns_dict['time_from_restart'] = True
        else:
            dns_dict['time_from_restart'] = False
        # eleventh line: nstep
        rawline = rawline_dnsin(f)
        dns_dict['nstep'] = float(rawline[0])
        
        return dns_dict


def rawline_dnsin(dnsin):
    preprocessed = dnsin.readline().split("!")[0].split()
    preprocessed = [entry for entry in preprocessed if not (entry=="" or entry=="\t")]
    return preprocessed