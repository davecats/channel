# This script takes as input a log from a channel run, be it from a slurm job
# or from piping the output of channel to a file.
# It calculates the average wall-time for a timestep and the average time needed
# to write a field on disk.

# Syntax:
# python /path/to/this/script.py /path/to/file.log

import sys
from channel import read_dnsin


# STEP 0: input handling and opening file
#########################################

if len(sys.argv) < 2:
    print("Please pass a filename as argument.")
    exit()
filename = sys.argv[1]

with open(filename) as f:
    lines = f.readlines()


# STEP 1: clean header and tail
###############################

# look for beginning
for ii, line in enumerate(lines):
    if 'nsteps =' in line:
        print("Found beginning of time iterations. Deleting unnecessary lines...")
        break
# check if beginning was found
if ii == (len(lines)-1):
    print("Beginning of time iterations was not found.")
# delete unnecessary lines
else:
    lines = lines[(ii+3):]

# look inside lines in reverse
ii = len(lines)
for line in reversed(lines):
    ii -= 1
    if 'End of time loop' in line:
        print("Found end of time iterations. Deleting unnecessary lines...")
        break
# check if beginning was found
if ii == 0:
    print("Beginning of time iterations was not found.")
# delete unnecessary lines
else:
    lines = lines[:ii]

# new line
print()


# STEP 2: writing time
######################

# make a list of lines to remove and a list of lines containing time steps after saving a file.
idx_todelete = []
idx_twrite = []
ii = 0
while ii < len(lines):
    if "Writing" in lines[ii]:
        idx_todelete.append(ii-1)
        idx_todelete.append(ii)
        ii += 1
        if "Writing" in lines[ii]:
            idx_todelete.append(ii)
            ii += 1
        idx_twrite.append(ii)
    ii += 1

avg_twrite = 0
for ii in idx_twrite:
    avg_twrite += float(lines[ii])
avg_twrite /= len(idx_twrite)
print('###')
print('###  On ' + str(len(idx_twrite)) + ' samples, average wall-time (seconds) for writing a field: ', avg_twrite)

# now delete all lines containing "Writing" or containing timestep for writing
for element in idx_twrite:
    idx_todelete.append(element)
lines = [i for j, i in enumerate(lines) if j not in idx_todelete]

# check if now number of remaining line is even
if not (len(lines) % 2) == 0:
    print('Something went wrong. Exiting.')
    exit()


# STEP 3: wall-time for normal timestep
#######################################

# calculate average
no_samples = 0
avg_wtstep = 0
avg_dtsim = 0
for ii, line in enumerate(lines):
    if not (ii % 2) == 0: # if the index is odd, then index corresponds to wall time for step
        avg_wtstep += float(line)
        no_samples += 1
    else: # if the index is even, you can calculate dt_sim
        temp = line.split(' ')
        temp = list(filter(lambda x: x!='', temp))
        avg_dtsim += float(temp[-1])

if not no_samples == (len(lines)/2): # quick consistency check on number of samples
    print('Something went wrong. Exiting.')
    exit()
avg_wtstep /= no_samples
avg_dtsim /= no_samples

# print out result
print('###  On ' + str(no_samples) + ' samples, average wall-time (seconds) for a standard time-step: ', avg_wtstep)
print('###')
print() # newline
print('On ' + str(no_samples) + ' samples, average time-step in simulation units: ', avg_dtsim)


# STEP 4: calculate time for a field
####################################

# get dt_field from dns.in
indns = read_dnsin('dns.in')
dtfield = indns['dt_field']

ts4fld = dtfield / avg_dtsim # no. of steps to write out a field; left purposely as a floating point number
print('Number of timesteps before writing out a field:', round(ts4fld))
wtime_4fld = ts4fld*avg_wtstep + avg_twrite # wall time to write out a field (seconds)

flds_hr = 1/(wtime_4fld/3600)

print() 
print('###')
print('###  Extimated no. fields per hour:', flds_hr)
print('###')
print()



