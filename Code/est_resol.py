#!/usr/bin/env python
# encoding: utf-8

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import mpl_tricks as mp
from COS_Vignetting_Throughput import tps
from scipy.ndimage import rotate
import scipy.interpolate as ipl
from astropy.io import fits

plt.rc('text', usetex=True)
plt.rc('font', **{'family':'serif', 'serif':['Computer Modern Roman'],
                  'sans-serif':'Helvetica'})

mirrored = [2, 6, 7, 12, 13, 14]

acqfiles = {1:  '../Data/lc3401luq_rawacq.fits',}

Apertures = {
    1:  [432, 529, 31.25]}

Cuts = {
   # 1:  [300, 700],
   # 2:  [150, 550],
    1:  None,
   #  4:  None,
   #  5:  [300, 500],
   #  6:  [100, 600],
   #  7:  [400, 600],
   #  8:  None,
   #  9:  None,
   #  10: [80, 600],
   #  11: None,
   #  12: [320,  450],
   #  13: [250,  500],
   #  14: [350,  430]
}


def find_angle(larsno):
    ''' Finds the COS y-axis rotation from in degrees E of N.
    '''
    filename = '../Data/haro11_stack.coarse.ascii'
    for line in open(filename):
        if line.startswith('#ORIENT'):
            angle = float(line.split()[2])
            print angle
    return angle


def circle_mask(image, cent_pix, radius):
    ''' Given an array, a center and a radius,
    generates a Bool array of the same dimension as input,
    which is True inside the radius and False outside.

    Parameters:
    -----------

    image: numpy.ndarray
        Two-dimensional numpy array, the entries of which are
        interpreted as image pixels.

    cent_pix: tuple
        (y, x) coordinates (in numpy column-major order) of the
        circle center.

    radius:  float
        radius in pixels of circle.

    Returns:
    --------

    circ_mask: numpy.ndarray
        Numpy array of datatype bool.
    '''
    centx, centy = cent_pix[1], cent_pix[0]
    x, y = np.arange(image.shape[1]), np.arange(image.shape[0])
    XX, YY = np.meshgrid(x, y)
    circ_mask = np.zeros_like(XX).astype(bool)
    circ_mask[np.sqrt((XX - centx) ** 2 + (YY - centy) ** 2) < radius] = True
    return circ_mask


def vignetting(image, cent_pix, radius):
    '''Emulates vignetting in the COS aperture'''
    arcsecbase = np.linspace(0, 1.25, 6)
    centx, centy = cent_pix[1], cent_pix[0]
    x, y = np.arange(image.shape[1]), np.arange(image.shape[0])
    XX, YY = np.meshgrid(x, y)
    vignette = np.zeros_like(image)
    dist_arcsec = np.sqrt((XX - centx) ** 2 + (YY - centy) ** 2) * 0.04
    throughspline = ipl.splrep(arcsecbase, tps.mean(0))
    vignette[dist_arcsec < 1.25] = dist_arcsec[dist_arcsec < 1.25]
    vig1d = vignette.flatten()
    vig1d = ipl.splev(vig1d, throughspline)
    vignette = vig1d.reshape(vignette.shape)
    vignette[dist_arcsec > 1.25] = 0.
    return vignette


def build_plot(larsno=None):
    fig = plt.figure(figsize=(6, 6), dpi=156)
    ax1 = plt.subplot(221)
    ax2 = plt.subplot(222)
    ax3 = plt.subplot(212)
    axes = [ax1, ax2, ax3]
    # plt.show()
    return fig, axes


def halflight_width(profile):
    fine_grid = np.linspace(0, len(profile), len(profile) * 10 + 1)
    fine_prof = sp.interp(fine_grid, np.arange(len(profile)), profile)
    flux = profile.sum()
    accu, step = 0., 0
    first, second = True, True
    limits = [0, 0]
    print flux
    for num in fine_prof:
        accu += num * 0.1
        step += .1
        if first:
            if (accu > flux * 0.25):
                print 'Flux 1: ', accu / flux
                limits[0] = step
                first = False
        if second:
            if (accu > 0.75 * flux):
                limits[1] = step
                second = False
                print 'Flux 2: ', accu / flux
                break
    limits[0] = np.around(limits[0])
    limits[1] = np.around(limits[1])
    print limits
    width = limits[1] - limits[0]
    return limits


def estimate_resol(larsno):
    COS_arcsec_per_resel = 0.171
    FUV_arcsec_per_pix = 0.04
    FUVframe = '../ImageData/haro_lyac.fits'
    FUVdata = 10 ** fits.getdata(FUVframe)
    rot_angle = find_angle(larsno)
    center = Apertures[larsno][:2]
    center = (center[1], center[0])
    radius = Apertures[larsno][2]
    centx, centy = center[1], center[0]
    mask = np.invert(circle_mask(FUVdata, center, radius))
    FUVorig = FUVdata.copy()
    FUVdata = np.ma.masked_where(mask, FUVdata) #[mask] = np.nan
    FUVdata *= vignetting(FUVdata, center, radius)
    FUVcut = FUVdata[centy-40:centy+40, centx-40:centx+40]
    FUVcut = rotate(FUVcut, rot_angle, reshape=False, order=5)
    cutmask = np.invert(circle_mask(FUVcut, (39.5, 39.5), radius))
    FUVcut = np.ma.masked_array(FUVcut, cutmask, fill_value=np.nan)
    profile = FUVcut.sum(0)
    limits = halflight_width(profile.data)
    width = np.array(limits).astype(float).ptp()
    print 'WIDTH:'
    print width
    width_arcsec = width * 0.04
    width_resels = width_arcsec / 0.171
    print width_arcsec, width_resels
    ###============================================================
    #   PLOTTING
    ###============================================================
    fig, axes = build_plot()
    vMax = np.sqrt(FUVorig).max() * .5
    axes[0].imshow(np.sqrt(FUVorig), interpolation='nearest',
                   cmap='gray', vmax=vMax, origin='lower')
    axes[0].imshow(np.sqrt(FUVdata), interpolation='nearest',
                   cmap='YlGnBu_r', vmax=vMax*1.5, origin='lower')#, cmap='hot')
    axes[0].tick_params(labelbottom='off', labelleft='off')
    if Cuts[larsno] is not None:
        c1, c2 = Cuts[larsno][0], Cuts[larsno][1]
        axes[0].axis((c1, c2, c1, c2))
    axes[1].imshow(np.sqrt(FUVcut), cmap='YlGnBu_r', origin='lower',
        interpolation='nearest')
    axes[1].tick_params(labelbottom='off', labelleft='off')
    axes[1].axis((6, 74, 6, 74))
    axes[2].plot(np.arange(len(profile)), profile, 'k-', drawstyle='steps-mid')# lw=1.6)
    mp.fill_between_steps(np.arange(limits[0], limits[1]+1),
        profile[limits[0]:limits[1]+1], color='gray', alpha=.4)
    axes[2].axis((0, 80, 0, profile.max() * 1.2))
    axes[2].set_ylabel('Flux [arbitrary units]')
    axes[2].set_xlabel(r'Pixel')
    axes[2].annotate('LARS {:02d}'.format(larsno), (0.05, 0.85),
        xycoords='axes fraction', color='k', fontsize='16',
        family='serif', weight='bold')
    S = 'Half-light width: \n Arcsec: {:.02f} \n COS Resels: {:.02f}'\
        .format(width * 0.04, width * 0.04 / 0.171)
    axes[2].annotate(S, (0.05, 0.83), xycoords='axes fraction', va='top')
    [j.set_linewidth(.6) for j in axes[1].spines.itervalues()]
    [j.set_linewidth(.6) for j in axes[2].spines.itervalues()]
    for ax in axes: ax.tick_params(length=2, labelsize=12)
    plt.subplots_adjust(hspace=.04, wspace=.04,
        left=.09, right=.96, top=.97, bottom=.07)
    #plt.savefig('../Figs/Resol.pdf')
    return width * 0.04 / 0.171


def main():
    resdict = {}
    for i in range(1, 2):
        resdict[i] = estimate_resol(i)
    return resdict


if __name__ == '__main__':
    main()
    plt.show()
