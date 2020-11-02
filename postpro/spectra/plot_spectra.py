import numpy as np
import h5py
from matplotlib.pyplot import subplots, imshow, colorbar, savefig
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap





###################
# CUSTOM COLORMAP #
###################

# colormap settings
white_length = 7 # as a percentage (%)
white_blending = 4 # as a percentage (%)

# preprocess
factor = 3
no_samples = 100 * factor
white_length *= factor
white_blending *= factor

# generate new colormap
inferno = cm.get_cmap('inferno_r')
newcolors = inferno(np.linspace(0, 1, no_samples - white_length))
white_array = np.ones((white_length,4)) # one row of this reads 1,1,1,1 - which is white in RGBA

# apply blending
if white_blending > 0:
    for ii in range(white_length-white_blending, white_length):
        weight = ii/(white_length+1)
        white_array[ii,:] *= (1-weight)
        white_array[ii,:] += weight*newcolors[0,:]

# generate colormap
newcolors = np.concatenate((white_array,newcolors),axis=0)
inferno_wr = ListedColormap(newcolors)





def plot_premultiplied(all_spectra, component, desired_y, y, kx, kz, **kwargs):

    # unpack input
    save_size = kwargs.get('save_size', (1,1)); w = save_size[0]; h = save_size[1]
    save_fig = kwargs.get('save_fig', False)
    cmp = gcmp(component) # get components the proper way
    y_idx = (np.abs(y - desired_y)).argmin() # identify index of y

    # function specific parameters
    labels = (r'$k_x$', r'$k_z$')
    xlabel_alt = r'$\lambda_x$'
    ylabel_alt = r'$\lambda_z$'
    fig_title = r'$\displaystyle k_xk_z \langle \hat{'+cmp[0]+r'}^\dagger \hat{'+cmp[1]+r'} \rangle$'
    xlog = True
    ylog = True
    save_name = cmp + '_premultiplied_y{}'.format(round(y[y_idx], 3))

    # convert component to index
    idx = get_comp_idx(component)
     
    # select desired spectrum
    spectrum = all_spectra[idx, :, :, :]

    # calculate premultiplier using array broadcasting
    bcast_premult = kz.reshape((-1, 1)) * kx
    bcast_premult = abs(bcast_premult)

    # premultiply
    premultiplied = spectrum * bcast_premult

    

    # prepare figure for plotting
    rfig, (rax,cb_ax) = subplots(ncols=2,figsize=(10,7),gridspec_kw={"width_ratios":[1, 0.05]})
    
    # set some scales to be logarithmic
    if xlog:
        rax.set_xscale("log")
    if ylog:
        rax.set_yscale("log")

    # when plotting, 0 modes are excluded
    pos = rax.pcolormesh(kx[1:], kz[1:], premultiplied[y_idx, 1:, 1:], linewidth=0, rasterized=True,shading='gouraud',cmap=inferno_wr)
    rax.set_xlabel(labels[0])
    rax.set_ylabel(labels[1])

    # add secondary x and y axes with wavelengths
    secaxx = rax.secondary_xaxis('top', functions=(get_wavelength_fwr,get_wavelength_inv))
    secaxx.set_xlabel(xlabel_alt)
    secaxy = rax.secondary_yaxis('right', functions=(get_wavelength_fwr,get_wavelength_inv))
    secaxy.set_ylabel(ylabel_alt)

    # last details for plotting
    rfig.colorbar(pos, cax=cb_ax)
    rax.set_title(fig_title + ', y = {}'.format(round(y[y_idx], 3)))

    # save
    if save_fig:
        fig = plt.figure(111,frameon=False)
        fig.set_size_inches(10*w,10*h)
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)
        if xlog:
            ax.set_xscale("log")
        if ylog:
            ax.set_yscale("log")
        ax.pcolormesh(kx[1:], kz[1:], premultiplied[y_idx, 1:, 1:], linewidth=0, rasterized=True,shading='gouraud',cmap=inferno_wr)
        # save figure
        savefig(save_name+'.png', format='png', bbox_inches='tight', pad_inches=0)
        # save tikz code
        generate_tikz(save_name, save_size, labels, fig_title, ax.get_xlim(), ax.get_ylim(), xlog=xlog, ylog=ylog)
        # close figure, so that it does not get plotted
        plt.close(111)
    
    return rfig, rax





def plot(all_spectra, component, desired_y, y, kx, kz, **kwargs):

    # unpack input
    vmin = kwargs.get('vmin', None)
    vmax = kwargs.get('vmax', None)
    cmp = gcmp(component) # get components the proper way
    y_idx = (np.abs(y - desired_y)).argmin() # identify index of y

    # function specific parameters
    labels = (r'$k_x$', r'$k_z$')
    fig_title = r'$\displaystyle \tilde{'+cmp[0]+r'}^\dagger \tilde{'+cmp[1]+'}$'

    # convert component to index
    idx = get_comp_idx(component)
    
    # select desired spectrum
    spectrum = all_spectra[idx, :, :, :]

    # actually print
    fig, (ax,cb_ax) = subplots(ncols=2,figsize=(10,7),gridspec_kw={"width_ratios":[1, 0.05]})
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    pos = ax.pcolormesh(kx,kz,spectrum[y_idx,:,:], vmin=vmin, vmax=vmax, linewidth=0, rasterized=True,shading='gouraud',cmap=inferno_wr)
    fig.colorbar(pos,cax=cb_ax)
    ax.set_title(fig_title + ', y = {}'.format(round(y[y_idx], 3)))





def plot_cumulative_zy(all_spectra, component, y, kz, **kwargs):

    # unpack input
    save_size = kwargs.get('save_size', (1,1)); w = save_size[0]; h = save_size[1]
    save_fig = kwargs.get('save_fig', False)
    y_symm = kwargs.get('y_symm', True)
    cmp = gcmp(component)

    # function specific parameters
    labels = (r'$k_z$', r'$y$')
    xlabel_alt = r'$\lambda_z$'
    fig_title = r'$\displaystyle k_z\sum_{k_x} \langle \hat{'+cmp[0]+r'}^\dagger\hat{'+cmp[1]+r'} \rangle$'
    xlog = True
    ylog = False
    save_name = cmp + '_cumulative_zy'

    # convert component to index
    idx = get_comp_idx(component)
    
    # select desired spectrum component
    spectrum = all_spectra[idx, :, :, :]

    # sum along x-axis
    spectrum = spectrum.sum(axis=(-1))

    # premultiply (reshape kz for broadcasting)
    premultiplied = spectrum * abs(kz)

    # prepare figure for plotting
    rfig, (rax,cb_ax) = subplots(ncols=2,figsize=(10,7),gridspec_kw={"width_ratios":[1, 0.05]})

    # set some scales to be logarithmic
    if xlog:
        rax.set_xscale("log")
    if ylog:
        rax.set_yscale("log")

    # when plotting, 0 modes are excluded
    pos = rax.pcolormesh(kz[1:], y, premultiplied[:,1:], linewidth=0, rasterized=True,shading='gouraud',cmap=inferno_wr)
    rax.set_xlabel(labels[0])
    rax.set_ylabel(labels[1])

    # add secondary x axis with wavelength
    secaxx = rax.secondary_xaxis('top', functions=(get_wavelength_fwr,get_wavelength_inv))
    secaxx.set_xlabel(xlabel_alt)

    # last details for plotting
    rfig.colorbar(pos, cax=cb_ax)
    rax.set_title(fig_title)
    if y_symm:
        rax.set_ylim([0,1])

    # save
    if save_fig:
        fig = plt.figure(111,frameon=False)
        fig.set_size_inches(10*w,10*h)
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)
        if xlog:
            ax.set_xscale("log")
        if ylog:
            ax.set_yscale("log")
        ax.pcolormesh(kz[1:], y, premultiplied[:,1:], linewidth=0, rasterized=True,shading='gouraud',cmap=inferno_wr)
        if y_symm:
            ax.set_ylim([0,1])
        # save figure
        savefig(save_name+'.png', format='png', bbox_inches='tight', pad_inches=0)
        # save tikz code
        generate_tikz(save_name, save_size, labels, fig_title, ax.get_xlim(), ax.get_ylim(), xlog=xlog, ylog=ylog)
        # close figure, so that it does not get plotted
        plt.close(111)

    return rfig, rax





def plot_cumulative_xy(all_spectra, component, y, kx, **kwargs):

    # unpack input
    save_size = kwargs.get('save_size', (1,1)); w = save_size[0]; h = save_size[1]
    save_fig = kwargs.get('save_fig', False)
    y_symm = kwargs.get('y_symm', True)
    cmp = gcmp(component)

    # function specific parameters
    labels = (r'$k_x$', r'$y$')
    xlabel_alt = r'$\lambda_x$'
    fig_title = r'$\displaystyle k_x\sum_{k_z} \langle \hat{'+cmp[0]+r'}^\dagger\hat{'+cmp[1]+r'} \rangle$'
    xlog = True
    ylog = False
    save_name = cmp + '_cumulative_xy'
        
    # convert component to index
    idx = get_comp_idx(component)
    
    # select desired spectrum component
    spectrum = all_spectra[idx, :, :, :]

    # sum along z-axis
    spectrum = spectrum.sum(axis=(-2))

    # premultiply (reshape kz for broadcasting)
    premultiplied = spectrum * abs(kx)

    # prepare figure for plotting
    rfig, (rax,cb_ax) = subplots(ncols=2,figsize=(10,7),gridspec_kw={"width_ratios":[1, 0.05]})
    
    # set some scales to be logarithmic
    if xlog:
        rax.set_xscale("log")
    if ylog:
        rax.set_yscale("log")

    # when plotting, 0 modes are excluded
    pos = rax.pcolormesh(kx[1:], y, premultiplied[:,1:], linewidth=0, rasterized=True,shading='gouraud',cmap=inferno_wr)
    rax.set_xlabel(labels[0])
    rax.set_ylabel(labels[1])

    # add secondary x axis with wavelength
    secaxx = rax.secondary_xaxis('top', functions=(get_wavelength_fwr,get_wavelength_inv))
    secaxx.set_xlabel(xlabel_alt)
    
    # last details for plotting
    rfig.colorbar(pos,cax=cb_ax)
    rax.set_title(fig_title)
    if y_symm:
        rax.set_ylim([0,1])

    # save
    if save_fig:
        fig = plt.figure(111,frameon=False)
        fig.set_size_inches(10*w,10*h)
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)
        if xlog:
            ax.set_xscale("log")
        if ylog:
            ax.set_yscale("log")
        ax.pcolormesh(kx[1:], y, premultiplied[:,1:], linewidth=0, rasterized=True,shading='gouraud',cmap=inferno_wr)
        if y_symm:
            ax.set_ylim([0,1])
        # save figure
        savefig(save_name+'.png', format='png', bbox_inches='tight', pad_inches=0)
        # save tikz code
        generate_tikz(save_name, save_size, labels, fig_title, ax.get_xlim(), ax.get_ylim(), xlog=xlog, ylog=ylog)
        # close figure, so that it does not get plotted
        plt.close(111)

    return rfig, rax





def gcmp(component):
    label = 'uu'
    component = component.lower()
    if component == 'y' or component == 'vv' or component == 'v':
        label = 'vv'
    elif component == 'z' or component == 'ww' or component == 'w':
        label = 'ww'
    elif component == 'xy' or component == 'yx' or component == 'uv' or component == 'vu':
        label = 'uv'
    elif component == 'xz' or component == 'zx' or component == 'uw' or component == 'wu':
        label = 'uw'
    elif component == 'zy' or component == 'yz' or component == 'wv' or component == 'vw':
        label = 'vw'
    elif not component == 'x' or component == 'uu' or component == 'u':
        print('Error: invalid component. Please pass either "x", "y", "z" as the second argument.')
    return label





def get_comp_idx(component):
    idx = 0
    component = component.lower()
    if component == 'y' or component == 'vv' or component == 'v':
        idx = 1
    elif component == 'z' or component == 'ww' or component == 'w':
        idx = 2
    elif component == 'xy' or component == 'yx' or component == 'uv' or component == 'vu':
        idx = 3
    elif component == 'xz' or component == 'zx' or component == 'uw' or component == 'wu':
        idx = 4
    elif component == 'zy' or component == 'yz' or component == 'wv' or component == 'vw':
        idx = 5
    elif not component == 'x' or component == 'uu' or component == 'u':
        print('Error: invalid component. Please pass either "x", "y", "z" as the second argument.')
    return idx

def get_wavelength_fwr(x):
    wvl = np.zeros_like(x)
    x_is_zero = x==0
    wvl[x_is_zero] = 1e+308
    wvl[~x_is_zero] = 2*np.pi/x[~x_is_zero]
    return wvl

def get_wavelength_inv(x):
    return 2*np.pi/x





def generate_tex_colormap():
    # newcolors is the array to print

    ncolors = newcolors.shape[0]

    text = r'''\pgfplotsset{
    colormap={custom_cmap}{
'''
    for ii in range(ncolors):
        text += '\t\t' 'rgb=(' + str(newcolors[ii,0]) + ',' + str(newcolors[ii,1]) + ',' + str(newcolors[ii,2]) + ')' + '\n'

    text += '''\t}
}'''

    with open('colorbar_settings.tex', 'w') as csets:
        csets.write(text)





def generate_tikz(fname, size, labels, title, xlim, ylim, **kwargs):

    xlog = kwargs.get('xlog', True)
    ylog = kwargs.get('ylog', True)

    width = size[0]; height = size[1]
    xlabel = labels[0]; ylabel = labels[1]


    xmode = 'xmode=log,'
    ymode = 'ymode=log,'

    if not xlog: xmode = '%' + xmode
    if not ylog: ymode = '%' + ymode

    tikz_code = f'''
\documentclass{{standalone}}
\\newcommand{{\\aver}}[1]{{\! \left\langle {{#1}} \\right\\rangle \!}}
\\usepackage{{pgfplots}}
\\usepackage{{pgfplotstable}}
\pgfplotsset{{compat=newest}}
\\usetikzlibrary{{pgfplots.colorbrewer}}
\input{{colorbar_settings}}
\\begin{{document}}
\\begin{{tikzpicture}}[]%
\\begin{{axis}}[enlargelimits=false,
             axis on top,
             width={width}\\linewidth, height={height}\\linewidth,
             xlabel={{{xlabel}}},
             {xmode}
             {ymode}
             ylabel={{{ylabel}}},
             view={{0}}{{90}},
             point meta min=0,point meta max=3.8, colorbar,
             colorbar style={{title={{{title}}},height=(\pgfkeysvalueof{{/pgfplots/parent axis height}}-3\\baselineskip),at={{(1.1,0)}},anchor=south west}},
             xmin={xlim[0]}, xmax={xlim[1]}, ymin={ylim[0]}, ymax={ylim[1]},
             scaled ticks=false,
            ]
\\addplot graphics [xmin={xlim[0]}, xmax={xlim[1]}, ymin={ylim[0]}, ymax={ylim[1]},] {{{fname}.png}}; % add external png
\end{{axis}}
\end{{tikzpicture}}
\end{{document}}
'''

    with open(fname+'.tikz', 'w') as f: # TODO: fixme!
        f.write(tikz_code)