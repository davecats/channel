import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import channel as ch
from math import sqrt
import pandas as pd

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
dudy = (25.4095 + 25.4133)/2
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