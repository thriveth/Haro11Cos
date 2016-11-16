#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from COS_Vignetting_Throughput import tps
import scipy as sp
import mpl_tricks as mp
from scipy.ndimage import rotate
from astropy.io import fits
import astropy.units as u
import scipy.interpolate as ipl
import pyregion as pr

Cuts = {1: None}

UVdata, UVhead = fits.getdata('../ImageData/haro_lyac.fits', header=True)
centerpix = [876.01019, 983.80147]
imangle = (np.arctan(1.0523680571667E-7/6.9436464044738E-6) * u.rad).to(u.deg)
COSangle = 101.39 * u.deg
totangle = imangle + COSangle
print imangle, COSangle, totangle
pixradius = 50.

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
    dist_arcsec = np.sqrt((XX - centx) ** 2 + (YY - centy) ** 2) * 0.025
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
    print 'Limits: ', limits
    width = limits[1] - limits[0]
    return limits


def estimate_resol(larsno):
    COS_arcsec_per_resel = 0.171  # Not actually used anywhere
    FUV_arcsec_per_pix = 0.025
    FUVframe = '../ImageData/haro_lyac.fits'
    FUVdata = fits.getdata(FUVframe)
    rot_angle = totangle.value  # find_angle(larsno)
    # center = Apertures[larsno][:2]
    # center = (center[1], center[0])
    center = centerpix[::-1]  # = (center[1], center[0])
    radius = pixradius  # Apertures[larsno][2]
    centx, centy = center[1], center[0]
    mask = np.invert(circle_mask(FUVdata, center, radius))
    FUVorig = FUVdata.copy()
    FUVdata = np.ma.masked_where(mask, FUVdata) #[mask] = np.nan
    FUVdata *= vignetting(FUVdata, center, radius)
    FUVcut = FUVdata[centy-80:centy+80, centx-80:centx+80]
    FUVcut = rotate(FUVcut, rot_angle, reshape=False, order=5)
    cutmask = np.invert(circle_mask(FUVcut, (79.5, 79.5), radius))
    FUVcut = np.ma.masked_array(FUVcut, cutmask, fill_value=np.nan)
    profile = FUVcut.sum(0)
    limits = halflight_width(profile.data)
    width = np.array(limits).astype(float).ptp()
    print 'WIDTH:', width
    width_arcsec = width * 0.025
    width_resels = width_arcsec / 0.171
    print width_arcsec, width_resels
    ###============================================================
    #   PLOTTING
    ###============================================================
    fig, axes = build_plot()
    vMax = np.sqrt(FUVorig).max() * .5
    axes[0].imshow(np.sqrt(FUVorig), interpolation='nearest',
                   cmap='gray', origin='lower', vmin=1e-23, vmax=5e-9)
    axes[0].imshow(np.sqrt(FUVdata), interpolation='nearest',
                   cmap='YlGnBu_r', origin='lower', vmin=1e-23, vmax=5.e-9)#, cmap='hot')
    # axes[0].imshow(np.sqrt(FUVdata), interpolation='nearest',
    #                cmap='viridis', origin='lower', vmin=1e-23, vmax=5.e-9)#, cmap='hot')
    axes[0].tick_params(labelbottom='off', labelleft='off')
    crcl = plt.Circle((centerpix), radius=50, color='w', lw=0.6, fill=False)
    axes[0].add_artist(crcl)
    axes[0].axis((810, 1130, 810, 1130))
    if Cuts[larsno] is not None:
        c1, c2 = Cuts[larsno][0], Cuts[larsno][1]
        axes[0].axis((c1, c2, c1, c2))
    axes[1].imshow(np.sqrt(FUVcut), cmap='YlGnBu_r', origin='lower',
        interpolation='nearest')
    # axes[1].imshow(np.sqrt(FUVcut), cmap='viridis', origin='lower',
    #     interpolation='nearest')
    axes[1].tick_params(labelbottom='off', labelleft='off')
    axes[1].axis((20, 140, 20, 140))
    crcl2 = plt.Circle((80, 80), radius=50, color='k', lw=2.5, fill=False)
    axes[1].add_artist(crcl2)
    profile /= profile.max()
    print len(profile)
    print limits, type(limits)
    limits = np.array(limits)
    #limits *= 0.025
    xes = np.arange(len(profile)) * 0.025
    xes -= xes.mean()
    print xes.min(), xes.max(), xes[31], xes[131], 'EXES!!'
    # axes[2].plot(np.linspace(-1.25, 1.25, len(profile)), profile,
    #              'k-', drawstyle='steps-mid')# lw=1.6)
    axes[2].plot(xes, profile,
                 'k-', drawstyle='steps-mid')# lw=1.6)
    #mp.fill_between_steps(np.arange(limits[0], limits[1]+1),
    #    profile[limits[0]:limits[1]+1], color='gray', alpha=.4)
    mp.fill_between_steps(xes[limits[0]:limits[1]+1],
        profile[limits[0]:limits[1]+1], color='gray', alpha=.4)
    axes[2].axis((-1.3, 1.30, 0, profile.max() * 1.22))
    axes[2].set_ylabel('Flux [arbitrary units]')
    axes[2].set_xlabel(r'Angle ["]')
    axes[2].annotate('Haro 11 C', (0.05, 0.85),
        xycoords='axes fraction', color='k', fontsize='16',
        family='serif', weight='bold')
    S = 'Half-light width: \n Arcsec: {:.02f} \n COS Resels: {:.02f}'\
        .format(width * 0.025, width * 0.025 / 0.171)
    axes[2].annotate(S, (0.05, 0.83), xycoords='axes fraction', va='top')
    [j.set_linewidth(.6) for j in axes[1].spines.itervalues()]
    [j.set_linewidth(.6) for j in axes[2].spines.itervalues()]
    for ax in axes: ax.tick_params(length=2, labelsize=12)
    plt.subplots_adjust(hspace=.04, wspace=.04,
        left=.10, right=.96, top=.97, bottom=.09)

    plt.savefig('../Figs/EffResol.pdf')
    return width * 0.025 / 0.171


def main():
    resdict = {}
    for i in range(1, 2):
        resdict[i] = estimate_resol(i)
    return resdict


if __name__ == '__main__':
    main()
    plt.show()




#ORIENTAT d -101.39  position angle of image y axis (deg. e of n)
