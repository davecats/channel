# copy this script in the root folder of a simulation (folder containing fields)
# and run it from there

import channel.postpro.spectra as ch
import matplotlib.pyplot as plt

all_spectra, kx, kz, y = ch.read_psd('.//')
ch.plot_premultiplied(all_spectra, 'x', 0.2, y, kx, kz)
ch.plot(all_spectra, 'x', 0.2, y, kx, kz)
ch.plot_cumulative_zy(all_spectra, 'x', y, kz, save_fig=True, save_size=(0.8,0.412))
ch.plot_cumulative_xy(all_spectra, 'x', y, kx)
plt.show()
ch.generate_tex_colormap()