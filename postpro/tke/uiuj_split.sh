if [[ $1 == "-h" ]]
then
    echo "Uses uiuj (or its derivatives, as uiuj_largesmall) to calculate statistics on a subrange of the whole set of snapshots. All available snapshots are however used to calculate the mean velocity field. Syntax:"
    echo
    echo "   uiuj_split.sh command arguments"
    echo
    echo "Where command can either be:"
    echo "- compute, compute_largesmall"
    echo "- merge, merge_largesmall"
    echo
    echo "The compute command calls uiuj to calculate statistics; the following arguments are required:"
    echo
    echo "   uiuj_split.sh compute nmin nmax with_mean mean_min mean_max [--sh_mem_proc nsmp]"
    echo
    echo "Where nmin, nmax describe the range of snapshots used for calculation of statistics, while mean_min, mean_max is the range of (all available) snapshots used for calculation of the mean velocity field. Flag -sh_mem_proc is used to set the number nsmp of shared memory processes."
    echo "Compute outputs a file named uiuj_nmin_nmax.bin in folder cm_profiles."
    echo "compute_largesmall is also a valid command; it does the same, but it invokes uiuj_largesmall instead. Output folder is cm_largesmall."
    echo
    echo "Command merge is used to merge multiple .bin files found in folder cm_profiles (or in folder cm_largesmall if merge_largesmall is used). Notice that merge must be called from the parent simulation directory (the one containing snapshots, which is the same compute is called from). Syntax is:"
    echo
    echo "   uiuj_split.sh merge"
    echo
    echo "notice that the .bin files need to retain the name assigned by compute, in particular they need to start with uiuj_. The output of merge is file uiuj.bin in folder profiles."
    exit
fi



# determine path of this script and of uiuj
script_dir=$(dirname $0)
uiujexec="${script_dir}/uiuj"
lasmexec="${script_dir}/uiuj_largesmall"


#############
# calculate #
#############

if [[ $1 == "compute" ]] # normal uiuj
then
    
    # store inputs with more readable names
    nmin=$2
    nmax=$3
    mean_min=$5
    mean_max=$6

    # compute
    $uiujexec 1 1 localhost $nmin $nmax 1 --custom_mean $mean_min $mean_max 1 ${@:7}

    # move files
    mv cm_profiles/uiuj.bin cm_profiles/uiuj_${nmin}_${nmax}.bin

    exit

fi

if [[ $1 == "compute_largesmall" ]] # large/small uiuj
then
    
    # store inputs with more readable names
    nmin=$2
    nmax=$3
    mean_min=$5
    mean_max=$6

    # compute
    mpirun --bind-to core --map-by core -report-bindings $lasmexec $nmin $nmax 1 --custom_mean $mean_min $mean_max 1

    # move files
    mv cm_largesmall/uiuj_largesmall.bin cm_largesmall/uiuj_largesmall_${nmin}_${nmax}.bin

    exit

fi



#########
# merge #
#########

if [[ $1 == "merge" ]] # for standard uiuj
then

    filelist=$(find cm_profiles -type f | grep uiuj_ | sed -e ':a; N; $!ba; s/\n/","/g')
    filelist='"'"${filelist}"'"'
    outfilemerge="profiles/uiuj.bin"

    mkdir -p profiles

    python - << ENDPYTHON
import numpy as np
import channel as ch

# fetch number of y points
ny=ch.get_dim_dnsin('./dns.in')[1]

# definition of datatypes
meanterms = np.dtype([('U', np.float64), ('W', np.float64),('Uy', np.float64), ('Wy', np.float64),('Uyy', np.float64), ('Wyy', np.float64),('P', np.float64)])
mketerms = np.dtype([('pump', np.float64), ('produv', np.float64),('prodvw', np.float64), ('ttrsp', np.float64),('vdiff', np.float64), ('dissU', np.float64),('dissW', np.float64),('PHIttrsp', np.float64),('PHIvdiff', np.float64)])
balance = np.dtype([('var', np.float64, (6,)),('prod', np.float64, (6,)),('psdiss', np.float64, (6,)),('ttrsp', np.float64, (6,)),('vdiff', np.float64, (6,)),('pstrain', np.float64, (6,)),('ptrsp', np.float64, (6,)),('PHIttrsp', np.float64, (6,)),('PHIvdiff', np.float64, (6,)),('PHIptrsp', np.float64, (6,))])
# datatype of binary file
bin_uiuj = np.dtype([('nmin', np.int32), ('nmax', np.int32), ('dn', np.int32), ('ntot', np.int32), ('meandata', meanterms, (ny+3,)), ('mkedata', mketerms, (ny+3,)), ('uiujdata', balance, (ny+3,))])

filez = [$filelist]

# allocation for index validity check
beginlist = []
endlist = []

# cumulative allocation
cumulative = np.zeros(1, dtype=bin_uiuj)

for file in filez: # loop over files

    #read data
    current = np.fromfile(file, dtype=bin_uiuj, count=1, sep="", offset=0)

    # generate lists for index validity check
    beginlist.append(current['nmin'])
    endlist.append(current['nmax'])

    # cumulate ntot
    cumulative["ntot"] += current["ntot"]

    # cumulate mketerms that need to be cumulated
    for mketerm in ['produv', 'prodvw', 'ttrsp', 'PHIttrsp']:
        cumulative['mkedata'][:][mketerm] += current['mkedata'][:][mketerm] * current["ntot"]

    # cumulate all uiuj terms
    for uiujterm in ["var","prod","psdiss","ttrsp","vdiff","pstrain","ptrsp","PHIttrsp","PHIvdiff","PHIptrsp"]:
        cumulative['uiujdata'][:][uiujterm][:] += current['uiujdata'][:][uiujterm][:] * current["ntot"]

# index validity check
beginlist.sort(); endlist.sort(); # sort lists
newmin = beginlist.pop(0) # find minimum and delete it
newmax = endlist.pop(-1) # find maximum and delete it
# check for repetitions of newmin or newmax
if newmin in beginlist:
    print('Warning: invalid set of nmin and nmax detected.')
if newmax in endlist:
    print('Warning: invalid set of nmin and nmax detected.')
# check everything else
for ii,st in enumerate(beginlist):
    if ii < (len(beginlist)-1):
        if beginlist[ii]==beginlist[ii+1] or endlist[ii]==endlist[ii+1]:
            print('Warning: invalid set of nmin and nmax detected.')
            break
    if (st-1) != endlist[ii]:
        print('Warning: invalid set of nmin and nmax detected.')
        break

# write dn, nmin, nmax
cumulative['nmin'] = newmin
cumulative['nmax'] = newmax
cumulative['dn'] = 1

# copy remaining mean terms
for meanterm in ["U","W","Uy","Wy","Uyy","Wyy","P"]:
    cumulative['meandata'][:][meanterm] =  current['meandata'][:][meanterm]

# copy remaining mketerms
for mketerm in ["pump", "vdiff", "dissU", "dissW", "PHIvdiff"]:
    cumulative['mkedata'][:][mketerm] =  current['mkedata'][:][mketerm]

# fix cumulative terms (divide by new ntot)
for mketerm in ['produv', 'prodvw', 'ttrsp', 'PHIttrsp']:
    cumulative['mkedata'][:][mketerm] /= cumulative["ntot"]
for uiujterm in ["var","prod","psdiss","ttrsp","vdiff","pstrain","ptrsp","PHIttrsp","PHIvdiff","PHIptrsp"]:
    cumulative['uiujdata'][:][uiujterm][:] /= cumulative["ntot"]

# output file
cumulative.tofile('$outfilemerge')

ENDPYTHON
    
    exit
    
fi


if [[ $1 == "merge_largesmall" ]]
then

    filelist=$(find cm_largesmall -type f | grep uiuj_largesmall_ | sed -e ':a; N; $!ba; s/\n/","/g')
    filelist='"'"${filelist}"'"'
    outfilemerge="profiles/uiuj_largesmall.bin"

    mkdir -p largesmall

    python - << ENDPYTHONLS
import numpy as np
import channel as ch

# fetch number of y points
ny=ch.get_dim_dnsin('./dns.in')[1]

# definition of datatypes
meanterms = np.dtype([('U', np.float64), ('W', np.float64),('Uy', np.float64), ('Wy', np.float64),('Uyy', np.float64), ('Wyy', np.float64),('P', np.float64)])
balance = np.dtype([('var', np.float64, (6,)),('prod', np.float64, (6,)),('psdiss', np.float64, (6,)),('ttrsp', np.float64, (6,)),('tcross', np.float64, (6,)),('vdiff', np.float64, (6,)),('pstrain', np.float64, (6,)),('ptrsp', np.float64, (6,)),('PHIttrsp', np.float64, (6,)),('PHIvdiff', np.float64, (6,)),('PHIptrsp', np.float64, (6,))])
# datatype of binary file
bin_uiuj = np.dtype([('nmin', np.int32), ('nmax', np.int32), ('dn', np.int32), ('ntot', np.int32), ('meandata', meanterms, (ny+3,)), ('luiujdata', balance, (ny+3,)), ('suiujdata', balance, (ny+3,))])

filez = [$filelist]

# allocation for index validity check
beginlist = []
endlist = []

# cumulative allocation
cumulative = np.zeros(1, dtype=bin_uiuj)

for file in filez: # loop over files

    #read data
    current = np.fromfile(file, dtype=bin_uiuj, count=1, sep="", offset=0)

    # generate lists for index validity check
    beginlist.append(current['nmin'])
    endlist.append(current['nmax'])

    # cumulate ntot
    cumulative["ntot"] += current["ntot"]

    # cumulate all uiuj terms
    for uiujterm in ["var","prod","psdiss","ttrsp","tcross","vdiff","pstrain","ptrsp","PHIttrsp","PHIvdiff","PHIptrsp"]:
        cumulative['luiujdata'][:][uiujterm][:] += current['luiujdata'][:][uiujterm][:] * current["ntot"]
        cumulative['suiujdata'][:][uiujterm][:] += current['suiujdata'][:][uiujterm][:] * current["ntot"]

# index validity check
beginlist.sort(); endlist.sort(); # sort lists
newmin = beginlist.pop(0) # find minimum and delete it
newmax = endlist.pop(-1) # find maximum and delete it
# check for repetitions of newmin or newmax
if newmin in beginlist:
    print('Warning: invalid set of nmin and nmax detected.')
if newmax in endlist:
    print('Warning: invalid set of nmin and nmax detected.')
# check everything else
for ii,st in enumerate(beginlist):
    if ii < (len(beginlist)-1):
        if beginlist[ii]==beginlist[ii+1] or endlist[ii]==endlist[ii+1]:
            print('Warning: invalid set of nmin and nmax detected.')
            break
    if (st-1) != endlist[ii]:
        print('Warning: invalid set of nmin and nmax detected.')
        break

# write dn, nmin, nmax
cumulative['nmin'] = newmin
cumulative['nmax'] = newmax
cumulative['dn'] = 1

# copy remaining mean terms
for meanterm in ["U","W","Uy","Wy","Uyy","Wyy","P"]:
    cumulative['meandata'][:][meanterm] =  current['meandata'][:][meanterm]

# fix cumulative terms (divide by new ntot)
for uiujterm in ["var","prod","psdiss","ttrsp","tcross","vdiff","pstrain","ptrsp","PHIttrsp","PHIvdiff","PHIptrsp"]:
    cumulative['luiujdata'][:][uiujterm][:] /= cumulative["ntot"]
    cumulative['suiujdata'][:][uiujterm][:] /= cumulative["ntot"]

# output file
cumulative.tofile('$outfilemerge')

ENDPYTHONLS
    
    exit
    
fi
