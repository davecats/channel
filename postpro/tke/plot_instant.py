import numpy as np
import channel as ch
from math import sqrt
import pandas as pd
from argparse import ArgumentParser

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.animation import FuncAnimation
matplotlib.use('qt5Agg')


# Parse arguments
argpar = ArgumentParser(description='Plot instantaneous statistics previously calculated by uiuj_instant.')
argpar.add_argument('first', metavar='nmin', type=int, nargs=1,
                    help='Number of first field.')
argpar.add_argument('last', metavar='nmax', type=int, nargs=1,
                    help='Number of last field.')
argpar.add_argument('fname', metavar='file', type=str, nargs=1,
                    help='Bare name of the file containing the desired statistics, e.g. "uu" or "uuintegrals".')
argpar.add_argument('column', metavar='col', type=str, nargs=1,
                    help='Column to be selected inside of file, e.g. "var".')
argpar.add_argument('-animate', help='Animate the plot of the selected statistics.', action='store_true')
argpar.add_argument('-scatter', help='Use a scatter plot; ignored by animate.', action='store_true')
settings = argpar.parse_args()
settings.first = settings.first[0]
settings.last = settings.last[0]
settings.fname = settings.fname[0]
settings.column = settings.column[0]


if "integral" in settings.fname:
    stat_is_profile = False
else:
    stat_is_profile = True


# get number of fields
noflds = settings.last - settings.first + 1
print("Number of fields used:", noflds)


if stat_is_profile: # if statistic is a profile

    # get viscous units
    #meshdata = ch.read_dnsin("dns.in")
    #re = meshdata['re']
    #mean = pd.read_csv('instant_profiles/mean.' + str(settings.nmin) + '.dat', header=1, delim_whitespace=True)
    #dudy = mean['Uy'].iat[0]
    #tauw = 1/re * dudy
    #utau = sqrt(tauw)
    #retau = re * utau
    #print("ReTau:", retau)
    #eta = 1/retau

    # generate colors
    clrmap = cm.get_cmap('copper_r')
    clrs = clrmap(np.linspace(0, 1, noflds))

    # prepare axes
    fig, ax = plt.subplots(figsize=(6,5))
    ax.set_title('instantaneous ' + settings.fname + ' ' + settings.column)
    ax.set_xlabel('y')

    # loop over files
    if settings.animate:
        ii = 0
        target_fname = 'instant_profiles/' + settings.fname + '.' + str(settings.first + ii) + '.dat'
        targ_panda = pd.read_csv(target_fname, header = 1, delim_whitespace=True)
        yseries = targ_panda['y']
        statistic = targ_panda[settings.column]
        graph = ax.plot(yseries, statistic) # first plot
        def one_frame(ii):
            # open desired filename
            target_fname = 'instant_profiles/' + settings.fname + '.' + str(settings.first + ii) + '.dat'
            targ_panda = pd.read_csv(target_fname, header = 1, delim_whitespace=True)
            yseries = targ_panda['y']
            statistic = targ_panda[settings.column]
            graph[0].set_data(yseries, statistic)
            ax.set_title('instantaneous ' + settings.fname + ' ' + settings.column + ', field ' + str(ii+1))
            return graph
        animation = FuncAnimation(fig, func=one_frame, frames=range(noflds), interval=100)
    else:
        for ii in range(noflds):

            # open desired filename
            target_fname = 'instant_profiles/' + settings.fname + '.' + str(settings.first + ii) + '.dat'
            targ_panda = pd.read_csv(target_fname, header = 1, delim_whitespace=True)
            yseries = targ_panda['y']
            statistic = targ_panda[settings.column]

            if settings.scatter:
                ax.scatter(yseries,statistic, color=clrs[ii])
            else:
                ax.plot(yseries,statistic, color=clrs[ii])

else: # if statistic is not a profile

    # prepare axes
    fig, ax = plt.subplots(figsize=(6,5))
    ax.set_title('instantaneous ' + settings.fname + ' ' + settings.column)
    ax.set_xlabel('field')

    # prepare data
    statistic = []
    flds = np.arange(settings.first, (settings.last+1))

    # get data
    for ii in range(noflds):

            # open desired filename
            target_fname = 'instant_profiles/' + settings.fname + '.' + str(settings.first + ii) + '.dat'
            targ_panda = pd.read_csv(target_fname, header = 1, delim_whitespace=True)
            statistic.append(targ_panda[settings.column].iat[0])
    
    if settings.scatter:
        ax.scatter(flds,statistic)
    else:
        ax.plot(flds,statistic)

plt.show()