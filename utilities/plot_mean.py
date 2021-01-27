from channel import field_as_memmap, read_sim_prop
from argparse import ArgumentParser
import matplotlib.pyplot as plt
from os.path import dirname
from matplotlib.animation import FuncAnimation

parser = ArgumentParser(description='Plots mean field of a turbulent field saved on disk. If more than one file is passed, plot is animated.')
parser.add_argument('fields', metavar='Field.ext', type=str, nargs='+', help='One or more fields, whose mean profile is to be plotted.')
settings = parser.parse_args()

fields = settings.fields

if len(fields) == 1:
    np, y, _, _ = read_sim_prop(fields[0], dirname(fields[0])+'dns.in')
    v_field = field_as_memmap(fields[0])
    mean = v_field[0,0,np[2],1:-1].real
    fig, ax = plt.subplots()
    ax.plot(y,mean)
    ax.set_title('mean U profile')
    plt.show()
else:
    np, y, _, _ = read_sim_prop(fields[0], dirname(fields[0])+'dns.in')
    v_field = field_as_memmap(fields[0])
    mean = v_field[0,0,np[2],1:-1].real
    fig, ax = plt.subplots()
    ii = 0
    graph = ax.plot(y, mean) # first plot
    def one_frame(ii):
        np, y, _, _ = read_sim_prop(fields[ii], dirname(fields[ii])+'dns.in')
        v_field = field_as_memmap(fields[ii])
        mean = v_field[0,0,np[2],1:-1].real
        graph[0].set_data(y, mean)
        ax.set_title('mean U profile, field '+str(ii))
        return graph
    animation = FuncAnimation(fig, func=one_frame, frames=range(len(fields)), interval=100)
    plt.show()

# Accesses a field on disk as a memmap, meaning it is not loaded to memory.
# FOR SAFETY REASONS, ACCESS IS READ ONLY BY DEFAULT.
# THE MEMMAP INCLUDES GHOST CELLS; this prevents unnecessary slicing.
# Usage:    v_field = field_as_memmap( 'Dati.cart.whatever.out' [, integer_size=4] )
# Arguments in square brackets are optional.
