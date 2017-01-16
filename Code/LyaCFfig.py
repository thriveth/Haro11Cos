#!/usr/bin/env python
# encoding: utf-8

import matplotlib as mpl
mpl.use('qt4agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pandas as pd
import cospaper as cosp
import mpl_tricks as mp
import wesanderson as wa
wes = [farve for farve in wa.colors['Darjeeling2']]# if farve != wa.colors['Darjeeling2'][0]]
cf_color = cosp.paircolors[3]
cf_color2 = cosp.paircolors[7]

wdata = pd.read_pickle('../Products/dfdict.pickle')[1]
vdata = pd.read_pickle('../Products/vflis.pickle')[1]
hdata = pd.read_pickle('../Products/vfhis.pickle')[1]

## OBS! Because we know the result of the AOD for Haro 11 COS,
## we set c1260 as proxy for covfrac.
vdata.MeanProfile = vdata.Si_II_1260_Flam
vdata.MeanErrors = vdata.Si_II_1260_Stddev
## End of dirty hardcoding. FIXME

z = 0.020434136
# print wdata.head()
# print vdata.head()
#
fig = plt.figure(figsize=(5.5, 5.5))

cont = plt.subplot(111, frameon=False)
cont.tick_params(labelleft='off', labelbottom='off')

# ax1 = plt.subplot2grid((3, 1), (0, 0), rowspan=2, )
# ax2 = plt.subplot2grid((3, 1), (2, 0))
ax1 = plt.subplot(211)
ax2 = plt.subplot(413)
ax3 = plt.subplot(414)
ax1.tick_params(labelsize=9)
ax1.label_outer()
ax2.label_outer()
ax2.tick_params(labelsize=9)
ax2.set_yticks([0, .5, 1.])
ax3.tick_params(labelsize=9)
ax3.set_yticks([0, .5, 1.])
ax2.axis((-850, 850, -.2, 1.6))
dfig, daxes = plt.subplots(2)
dax1 = daxes[0]
dax2 = daxes[1]
cont = fig.add_subplot(111, frameon=False)
cont.tick_params(labelleft='off', labelbottom='off', length=0)
cont.set_xlabel(r'$v - v_0$ [km/s]', labelpad=15)
cont.set_ylabel(r'\rm Relative flux', labelpad=18)

twax = plt.twinx(ax2)
twax.set_ylim(1.2, -.6)
twax.set_yticks([0, .5, 1])
twax.set_ylabel(r'$\mathbf{f_c(v)}$', color=cf_color)
twax.tick_params(color=cf_color, labelcolor=cf_color, labelsize=9, length=3)
twax1 = plt.twinx(ax3)
twax1.set_ylim(1.2, -.6)
twax1.set_yticks([0, .5, 1])
twax1.set_ylabel(r'$\mathbf{f_c(v)}$', color=cf_color2)
twax1.tick_params(color=cf_color2, labelcolor=cf_color2, labelsize=9, length=3)
[j.set_linewidth(.6) for j in ax2.spines.itervalues()]
[j.set_linewidth(.6) for j in ax1.spines.itervalues()]
[j.set_linewidth(.6) for j in ax3.spines.itervalues()]
[j.set_linewidth(.6) for j in twax.spines.itervalues()]
[j.set_linewidth(.6) for j in twax1.spines.itervalues()]
cosp.draw_coverfrac_axes(
    vdata, dax1, dax2, ax2,
    thin=.3, legend=False, ms=3.,
    dotcolor=cosp.paircolors[3],
    barcolor=cosp.paircolors[2],
    cf_lims=[-700, 500])
cosp.draw_coverfrac_axes(
    hdata, dax1, dax2, ax3,
    thin=.3, legend=False, ms=3.,
    species='Si_IV',
    dotcolor=cosp.paircolors[7],
    barcolor=cosp.paircolors[6],
    cf_lims=[-700, 500])

lyawave = 1215.668 * (1. + z)
wdata = wdata[wdata['lambda'].between(lyawave-10, lyawave+10)]
restwave = wdata['lambda'] / (1 + z)
lyavels = (restwave - 1215.668) / 1215.668 * sp.constants.c / 1000
interval = lyavels[(lyavels > -3000) & (lyavels < 3000)]
topval2 = wdata['flam'].median()
topval = wdata['flam'].loc[(lyavels > -1600) & (lyavels < 1600)].max()/topval2
topval = max(topval, 1.)
ax1.plot(lyavels, wdata['flam']/topval2, lw=.9, color=cosp.paircolors[1], drawstyle='steps-mid')
#ax1.plot(lyavels, wdata['flam']/topval2, lw=.9, color='darkblue', drawstyle='steps-mid')
mp.fill_between_steps(lyavels.values,  # Hacky-hacky: with white first.
                      (wdata['flam'].values-wdata['stddev'].values)/topval2,
                      (wdata['flam'].values+wdata['stddev'].values)/topval2,
                      facecolor='w', edgecolor='w', ax=ax1, lw=.1)
mp.fill_between_steps(lyavels.values,
                          (wdata['flam'].values-wdata['stddev'].values)/topval2,
                          (wdata['flam'].values+wdata['stddev'].values)/topval2,
                          facecolor=cosp.paircolors[0],
                          edgecolor=cosp.paircolors[0], # '.3',#wes[0],
                          ax=ax1, lw=.2, alpha=.99)
ax1.axis((-1200, 1200, -.2*topval, 1.20*topval))
ax2.axis((-1200, 1200, -.2, 1.60))
ax3.axis((-1200, 1200, -.2, 1.60))
ax1.axvline(x=0., lw=.5, ls='--', color='k', zorder=0)
ax1.axhline(y=0., lw=.5, ls='--', color='k', zorder=0)
ax1.yaxis.set_ticks(ax1.yaxis.get_ticklocs()[sp.where(ax1.yaxis.get_ticklocs() >= 0.)])
ax1.annotate('\sffamily{Haro 11 C}', xy=(.05, .88), xycoords='axes fraction',
             size=13, fontstyle='italic')
ax1.annotate(r'Ly$\alpha$', xy=(.10, .60), xycoords='axes fraction', size=18,
             color=cosp.paircolors[1])#'darkblue')#cosp.paircolors[3])
ax2.annotate(r'$f_C$ (Si \textsc{II})', xy=(.01, .2), xycoords='axes fraction', size=13,
             color=cosp.paircolors[3])
ax2.annotate(r'Avg. LIS', xy=(.75, .2), xycoords='axes fraction', size=13,
             color='gray')
ax3.annotate(r'$f_C$ (Si \scshape{IV})', xy=(.01, .2), xycoords='axes fraction', size=13,
             color=cosp.paircolors[7])
ax3.annotate(r'Avg. Si IV', xy=(.75, .2), xycoords='axes fraction', size=13,
             color='gray')
fig.subplots_adjust(bottom=.11, top=.95, hspace=.02, wspace=.02, left=.10, right=.89)
fig.savefig('../Figs/LyACoverfracs.pdf')
fig.savefig('../Figs/LyACoverfracs.png')
plt.show()

