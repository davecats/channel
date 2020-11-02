import pandas as pd
from dns_channel.dnsdata import read_dnsin



def read(fdir, **kwargs):

    # runtimePanda, mePanda, uvPanda, uuPanda, vvPanda, wwPanda, kkPanda, mkPanda = read(fdir, invert_sss=False, variant=None, head_len=1)
    # Reads Runtimedata and output of uiuj; returns a tuple of Pandas dataframes

    # INPUT DESCRIPTION
    # fidr: string - directory where output of uiuj and Runtimedata is contained
    invert_sss = kwargs.get('invert_sss', False)   # if true, inverts the sign of shear stress on top wall from runtimedata
    variant = kwargs.get('variant', None)   # variant of uiuj used for postprocessing
    head_len = kwargs.get('head_len', 1) # length of header (notice that the column names should be OUTSIDE OF HEADER!)


    # filenames

    filenames = ['mean', 'uv', 'uu', 'vv', 'ww', 'tke', 'mke']
    for ii in range(len(filenames)):
        if not variant == None:
            filenames[ii] += '_' + variant 
        filenames[ii] += '.dat' 


    # read Reynolds stress tensor components from uiuj

    uvPanda = pd.read_csv(fdir + filenames[1], header = head_len, delim_whitespace=True) # reynolds stresses
    uuPanda = pd.read_csv(fdir + filenames[2], header = head_len, delim_whitespace=True) # kinetic energy: uu
    vvPanda = pd.read_csv(fdir + filenames[3], header = head_len, delim_whitespace=True) # kinetic energy: vv
    wwPanda = pd.read_csv(fdir + filenames[4], header = head_len, delim_whitespace=True) # kinetic energy: ww
    kkPanda = pd.read_csv(fdir + filenames[5], header = head_len, delim_whitespace=True) # turbulent kinetic energy
    

    # read mke and runtimedata

    if variant=='large' or variant =='small':
        
        mkPanda = None
        mePanda = None
        runtimePanda = None
    
    else:
        
        mkPanda = pd.read_csv(fdir + filenames[6], header = head_len, delim_whitespace=True) # mean kinetic energy
        mePanda = pd.read_csv(fdir + filenames[0], header = head_len, delim_whitespace=True) # mean field
        
        # read Runtimedata

        runtimePanda = pd.read_csv(fdir + 'Runtimedata', header=None, delim_whitespace=True)
        _, cols = runtimePanda.shape

        if cols == 11: # column names are assigned depending on which program generated Runtimedata
            runtimePanda.columns = ['t', 'Uy_bottom', 'Uy_top', 'Wy_bottom', 'Wy_top', 'Ub', 'Px', 'Ut', 'Pt', 'CFL', 'dt']
        else:
            runtimePanda.columns = ['t', 'Uy_bottom', 'Uy_top', 'Wy_bottom', 'Wy_top', 'Ub', 'Px', 'Ut', 'Pt', 'CFL', 'dt', 'non', 'so'] # FIXME: controlla

        if invert_sss: # if requested by user, invert sign of shear stress at top wall
            runtimePanda['Uy_top'] = -runtimePanda['Uy_top']


    return runtimePanda, mePanda, uvPanda, uuPanda, vvPanda, wwPanda, kkPanda, mkPanda



def read_integrals(fdir, **kwargs):

    # runtimePanda, mePanda, uvPanda, uuPanda, vvPanda, wwPanda, kkPanda, mkPanda = read(fdir, invert_sss=False, variant=None, head_len=1)
    # Reads Runtimedata and output of uiuj; returns a tuple of Pandas dataframes

    # INPUT DESCRIPTION
    # fidr: string - directory where output of uiuj and Runtimedata is contained
    variant = kwargs.get('variant', None)   # variant of uiuj used for postprocessing
    head_len = kwargs.get('head_len', 1) # length of header (notice that the column names should be OUTSIDE OF HEADER!)

    # filenames
    filenames = ['uv', 'uu', 'vv', 'ww', 'tke', 'mke']
    filenames = [filenames[ii] + 'integrals' for ii in range(len(filenames))]
    for ii in range(len(filenames)):
        if not variant == None:
            filenames[ii] += '_' + variant 
        filenames[ii] += '.dat' 

    # read Reynolds stress tensor components from uiuj
    uvPanda = pd.read_csv(fdir + filenames[0], header = head_len, delim_whitespace=True) # reynolds stresses
    uuPanda = pd.read_csv(fdir + filenames[1], header = head_len, delim_whitespace=True) # kinetic energy: uu
    vvPanda = pd.read_csv(fdir + filenames[2], header = head_len, delim_whitespace=True) # kinetic energy: vv
    wwPanda = pd.read_csv(fdir + filenames[3], header = head_len, delim_whitespace=True) # kinetic energy: ww
    kkPanda = pd.read_csv(fdir + filenames[4], header = head_len, delim_whitespace=True) # turbulent kinetic energy 

    # read mke
    if variant=='large' or variant =='small':
        mkPanda = None
    else:
        mkPanda = pd.read_csv(fdir + filenames[5], header = head_len, delim_whitespace=True) # mean kinetic energy

    # create list of panda dataframes and initialise dictionaries
    pandalist = [uvPanda, uuPanda, vvPanda, wwPanda, kkPanda, mkPanda]
    uv = {}; uu = {}; vv = {}; ww = {}; kk = {}; mk = {}
    dictlist = [uv, uu, vv, ww, kk, mk]
    for pnd,dct in zip(pandalist,dictlist):
        if not pnd is None:
            for term in pnd:
                dct.update({term:pnd[term].iat[0]})

    return uv, uu, vv, ww, kk, mk



def readLM(fabel, **kwargs):

    # file directory
    base_dir = kwargs.get('base_dir', './/')
    fdir = base_dir + fabel + '_'

    mePanda = pd.read_csv(fdir + 'mean.dat', header = 70, delim_whitespace=True, engine='python')
    flPanda = pd.read_csv(fdir + 'fluc.dat', header = 73, delim_whitespace=True, engine='python')
    #uvPanda = pd.read_csv(fdir + 'uv.dat', header = 70, delim_whitespace=True, engine='python') # reynolds stresses
    #uuPanda = pd.read_csv(fdir + 'uu.dat', header = 72, delim_whitespace=True, engine='python') # kinetic energy: uu
    #vvPanda = pd.read_csv(fdir + 'vv.dat', header = 72, delim_whitespace=True, engine='python') # kinetic energy: vv
    #wwPanda = pd.read_csv(fdir + 'ww.dat', header = 75, delim_whitespace=True, engine='python') # kinetic energy: ww
    #kkPanda = pd.read_csv(fdir + 'k.dat', header = 72, delim_whitespace=True, engine='python') # kinetic energy

    return mePanda, flPanda#, uvPanda, uuPanda, vvPanda, wwPanda, kkPanda



def get_nfield(fdir):

    with open(fdir + 'uu.dat') as infile:
        fline = infile.readline()

    tmp = ''

    for char in fline:
        if char.isdigit():
            tmp = tmp + char
        else:
            tmp = tmp + ' '

    tmp = tmp.split()

    return int(tmp[-1])



def get_z_threshold(fdir):

    dct = locals()

    with open(fdir + 'largesmall_settings.in', 'r') as infile:
        exec(infile.readline(), globals(), dct)

    return dct['z_threshold']