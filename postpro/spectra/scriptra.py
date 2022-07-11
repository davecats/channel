# copy this script in the root folder of a simulation (folder containing fields)
# and run it from there

import channel as ch
import channel.postpro.spectra as psd
import matplotlib.pyplot as plt

m = ch.mesh(ch.read_dnsin('../dns.in'))
all_spectra, kx, kz, y = psd.read_psd('psd.bin')
psd.plot_premultiplied(all_spectra, 'x', 0.2, y, kx, kz)
psd.plot(all_spectra, 'x', 0.2, y, kx, kz)
psd.plot_cumulative_zy(all_spectra, 'x', y, kz, m.beta0, save_fig=True, save_size=(0.8,0.412))
plt.show()
psd.generate_tex_colormap()