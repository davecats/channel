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

argpar = ArgumentParser(description='Plot the variance history previously calculated by var_history.')
argpar.add_argument('--animate', help=   '''Animate the plot of instantaneous variances of selected components.
                                            For instance, "--animate uv" creates the animation for components u
                                            and v, while "--animate w" only animates the w component.''')
settings = argpar.parse_args()

# read postpro
cumul = np.fromfile("var_history/cumulative.bin", dtype=np.float64)
insta = np.fromfile("var_history/instantaneous.bin", dtype=np.float64)

# create dictionary for dnsdata
meshdata = ch.read_dnsin("dns.in")
re = meshdata['re']
ny = meshdata['ny']

# generate mesh data
m = ch.mesh(meshdata)
m.y = m.y[1:-1] # remove ghost cells

# get viscous units
with open('var_history/dudy.txt') as f:
    dudy = float(f.readline())
    #dudy = 500
tauw = 1/re * dudy
utau = sqrt(tauw)
retau = re * utau
print("ReTau:", retau)
eta = 1/retau
# switch to viscous units
cumul /= tauw
insta /= tauw
m.y /= eta
 
# reshape arrays
cumul = cumul.reshape((-1,3,ny+1))
insta = insta.reshape((-1,3,ny+1))

# get number of fields
noflds = np.size(cumul, axis=0)
print("Number of fields used:", noflds)

# generate colors
clrmap = cm.get_cmap('copper_r')
clrs = clrmap(np.linspace(0, 1, noflds))

# cumulative: inverse
invcum = np.zeros_like(cumul)
for ii in range(noflds):
    for jj in range(ii,noflds):
        invcum[ii,:,:] += insta[jj,:,:]
    invcum[ii,:,:] /= noflds - ii

# actual plot, uu
fig, (ax, sax, tax) = plt.subplots(1,3, figsize=(18,5))
fig.suptitle('uu fluctuation')
ax.set_title('instantaneous')
sax.set_title('cumulative')
tax.set_title('cumulative inverse')
for ii in range(noflds):
    ax.plot(m.y,insta[ii,0,:], color=clrs[ii])
    sax.plot(m.y,cumul[ii,0,:], color=clrs[ii])
    tax.plot(m.y,invcum[ii,0,:], color=clrs[ii])

# actual plot, vv
fig, (ax, sax, tax) = plt.subplots(1,3, figsize=(18,5))
fig.suptitle('vv fluctuation')
ax.set_title('instantaneous')
sax.set_title('cumulative')
tax.set_title('cumulative inverse')
for ii in range(noflds):
    ax.plot(m.y,insta[ii,1,:], color=clrs[ii])
    sax.plot(m.y,cumul[ii,1,:], color=clrs[ii])
    tax.plot(m.y,invcum[ii,1,:], color=clrs[ii])

# actual plot, vv
fig, (ax, sax, tax) = plt.subplots(1,3, figsize=(18,5))
fig.suptitle('ww fluctuation')
ax.set_title('instantaneous')
sax.set_title('cumulative')
tax.set_title('cumulative inverse')
for ii in range(noflds):
    ax.plot(m.y,insta[ii,2,:], color=clrs[ii])
    sax.plot(m.y,cumul[ii,2,:], color=clrs[ii])
    tax.plot(m.y,invcum[ii,2,:], color=clrs[ii])

# show all plots
plt.show()

# now animate if requested
if settings.animate:
    
    # first off, define parse components
    cmpnt = []
    for char in settings.animate:
        if char == 'u':
            cmpnt.append(0)
        if char == 'v':
            cmpnt.append(1)
        if char == 'w':
            cmpnt.append(2)

    # now animate
    for it,c in enumerate(cmpnt):
        fig, ax = plt.subplots()
        ii = 0
        graph = ax.plot(m.y, insta[0,c,:]) # first plot
        ax.set_title(settings.animate[it]+settings.animate[it]+' fluctuation, field '+str(1))
        def one_frame(ii):
            graph[0].set_data(m.y, insta[ii,c,:])
            ax.set_title(settings.animate[it]+settings.animate[it]+' fluctuation, field '+str(ii))
            return graph
        animation = FuncAnimation(fig, func=one_frame, frames=range(1,noflds), interval=100)
        plt.show()