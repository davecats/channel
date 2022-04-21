from channel import field_as_memmap, read_sim_prop
from argparse import ArgumentParser
import matplotlib.pyplot as plt
from os import mkdir
from os.path import dirname, isdir
from matplotlib.animation import FuncAnimation
import numpy

parser = ArgumentParser(description='Retrieves time history of the mean velocity from a given set of fields and saves it as ASCII in folder "mean_history". Optionally plots the output. If a single file is passed, no output is written to disk - rather, this script defaults to plotting.')
parser.add_argument('bounds', nargs=2, metavar='nfield', type=int, help='Upper and lower bounds for set of fields to be analysed.')
parser.add_argument('--plot_average', '-pa', help='If a list of files is provided, computes and plot the ensemble average.', action='store_true')
parser.add_argument('--plot', '-p', help='Plots and animates the time history of the mean profile in selected files.', action='store_true')
settings = parser.parse_args()

fields = ['Dati.cart.'+str(ii)+'.out' for ii in range(settings.bounds[0],settings.bounds[1]+1)]
np, y, _, _ = read_sim_prop(fields[0], dirname(fields[0])+'dns.in')

# allocate array and read
hist = numpy.zeros((len(fields),np[1]+3))
for ii,fld in enumerate(fields):
    v_field = field_as_memmap(fld)
    hist[ii,:] = v_field[0,0,np[2],:].real

# save to disk
if not len(fields) == 1:
    # check existence of directory; create if necessary
    if not isdir('mean_history'):
        mkdir('mean_history')
    # write out
    numpy.savetxt('mean_history/u_hist.txt', hist, delimiter='  ')

if len(fields) == 1:
    mean = hist[0,1:-1] # this eliminates ghost cells
    fig, ax = plt.subplots()
    ax.plot(y,mean)
    ax.set_title('mean U profile')
    plt.show()
elif settings.plot:
    ii = 0
    mean = hist[ii,1:-1] # this eliminates ghost cells
    fig, ax = plt.subplots()
    graph = ax.plot(y, mean) # first plot
    def one_frame(xx):
        mean = hist[xx,1:-1]
        graph[0].set_data(y, mean)
        ax.set_title(fields[xx])
        return graph
    animation = FuncAnimation(fig, func=one_frame, frames=range(len(fields)), interval=100)
    plt.show()
elif settings.plot_average:
    mean = numpy.sum(hist,axis=0) / len(fields)
    mean = mean[1:-1] # this eliminates ghost cells
    fig, ax = plt.subplots()
    ax.plot(y, mean)
    plt.show()