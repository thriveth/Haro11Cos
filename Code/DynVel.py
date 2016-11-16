#! /usr/bin/env/ python
# -*- coding: utf-8 -*-

import matplotlib as mpl
#mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#mpl.rc('font',**{'family':'serif','serif':['Palatino']})
mpl.rc('text', usetex=True)     # For pretty LaTeX fonts in plots
import matplotlib.gridspec as grisp
import matplotlib.pyplot as plt
import numpy as np
#import os
#import pyfits as pf
import scipy as sp
#import scipy.constants as const
import mpl_tricks as mpt

#spedata = pf.getdata('./NormVIS.fits')
#width = pf.getheader('./NormVIS.fits')
#errdata = pf.getdata('./NormVIS.err.fits')
#errhead = pf.getheader('./NormVIS.err.fits')
#
#xs = np.arange(errdata.shape[0])
#wl = xs*width['CD1_1']+width['CRVAL1']
clickcount = 0
exes = []


def onclick(event):
    global clickcount
    global exes
    print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f' % (
        event.button, event.x, event.y,
        event.xdata, event.ydata)
    if event.button == 1:
        exes.append(event.xdata)
        if((clickcount <= 1) | (clickcount > 3)):
            linecolor = 'green'
        else:
            linecolor = 'red'
        if(clickcount == 1):
            plt.title('Click on left border of line')
        elif(clickcount == 2):
            plt.title('Click on right border of line')
        elif(clickcount > 2):
            plt.title('Click twice on right continuum for calibration')
        plt.plot([event.xdata, event.xdata], [.9, 1.1],
                 color=linecolor, linewidth=2)
        plt.draw()
        clickcount += 1
        #print clickcount
        if((clickcount == 2) | (clickcount == 4) | (clickcount == 6)):
            #print exes[clickcount-2:clickcount]
            plt.plot(exes[clickcount - 2:clickcount],
                     [1, 1], linecolor, linewidth=2)
            plt.draw()
        if(clickcount > 6):
            plt.close()


def vel_bin_size(indict):
    velocities = indict.Velocity
    vb1 = velocities.values.ptp() / float(len(velocities) - 1)
    velbinlist = []
    for i in range(len(velocities) - 1):
        width = velocities[i + 1] - velocities[i]
        velbinlist.append(width)
    velbarr = sp.array(velbinlist)
    return velbarr.mean()


def dynvel(indict, ranges=None, species='Average LIS',
           savefig=True, plot=True):
    ''' Takes velocity-space dataframe dict as input.
    '''
    if plot and not savefig:
        print('Plot is set to False, so no figure can be saved')
    plt.rc('text', usetex=True)
    plt.ioff()
    larsname_ul = indict.columns.name.replace(' ', '_')
    spedata = indict.MeanProfile#.values
    errdata = indict.MeanErrors#.values  # sp.ones_like(spedata)
    wl = indict.Velocity.values
    if ranges is None:
        fig = plt.figure()  # figsize=(3.4, 3.4), dpi=160)
        ax = plt.subplot(111)
        #ax.tick_params(labelsize=8)
        #[j.set_linewidth(.5) for j in ax.spines.itervalues()]
        exes = []
        plt.plot(wl, spedata, 'black', drawstyle='steps-mid')
        plt.axis((-800, 800, -.1, 1.4))
        plt.axhline(y=1., linestyle='--', color='k')
        plt.annotate(r'''Mark piece of left continuum, beginning and end of line,\\
                     and piece of right continuum, 6 clicks in total''',
                     (.5, .9), xycoords='axes fraction',
                     ha='center')
        exes = plt.ginput(6)
        plt.close()
        leftexes = np.asarray(exes[0:2])[:, 0]
        riteexes = np.asarray(exes[4:6])[:, 0]
        lineexes = np.asarray(exes[2:4])[:, 0]
        exesarr = np.asarray(exes)
    else:
        try:
            exesarr = sp.array(ranges).reshape(3, 2)
            leftexes = exesarr[0, :]
            lineexes = exesarr[1, :]
            riteexes = exesarr[2, :]
        except:
            print('Could not make sense of given ranges, exiting')
            return
    spedata_orig = spedata.copy()
    spedata = spedata[indict.Velocity.between(lineexes[0], lineexes[1])]

    line = spedata[indict.Velocity.between(lineexes[0], lineexes[1])]
    lidx = sp.where((wl >= lineexes.min()) & (wl <= lineexes.max()))
    liwl = wl[(wl >= lineexes.min()) & (wl <= lineexes.max())]
    zero = sp.where(sp.absolute(wl) == sp.absolute(wl).min())[0]

    norm = 1.      # Because that is what I did when finding robust lower limits
    width = vel_bin_size(indict)
    area = (norm-line)*width
    #datamin = spedata[indict.Velocity.between(lineexes.min(),
    #                                          lineexes.max())].min()
    datamin = spedata.min()
    datazero = spedata[zero].values#[0]
    minwl = wl[spedata[indict.Velocity.between(lineexes.min(),
                                               lineexes.max())].argmin()]

    accu = 0.
    #Array to be plotted to illustrate how area changes:
    accuaccu = np.zeros(line.shape[0])
    lomark = True
    himark = True
    for i in np.arange(line.shape[0]):
        accu += (norm-line.values[i])*width
        accuaccu[i] = accu
        if lomark==True:
            if accu >= area.sum()*.05:
                lomark = False
                xlo = liwl[i]
                #print 'Percentage at xlo: ', accu/area.sum()
                #print 'Percentage at pixel before: ', \
                #              accuaccu[i-1]/area.sum()
        if((lomark == False) & (himark == True)):
            if accu >= area.sum()*.95:
                himark = False
                xhi = liwl[i]
                #print 'Percentage at xhi: ', accu/area.sum()
                #print 'Percentage at pixel before: ', \
                #                accuaccu[i-1]/area.sum()
    #print xlo, xhi
    ###======================================================================
    #      First moment:
    M1 = (liwl*(norm-line)/errdata[(wl >= lineexes.min()) &
                                (wl <= lineexes.max())]**2).sum()/\
                    ((norm-line)/ errdata[(wl >= lineexes.min()) &
                                        (wl <= lineexes.max())]**2).sum()
    #print 'Line center (AA): \t', M1
    llo = lineexes.min()
    lhi = lineexes.max()
    #DynWidth = (xhi-xlo)/M1*const.c
    DynWidth = xhi - xlo
    if plot:
        #print 'Dynamic velocity:\t', DynWidth
        pltmin = np.mean([leftexes.min(), leftexes.max()])
        pltmax = np.mean([riteexes.min(), riteexes.max()])
        pltmin = -900 # JUST FOR THE PAPER
        pltmax = 700  # JUST FOR THE PAPER
        #print pltmin, pltmax
        gs = grisp.GridSpec(4, 4)
        fig2 = plt.figure(figsize=(3.4, 3.), dpi=160)
        ax1 = plt.subplot(gs[0, :])
        ax1.tick_params(labelsize=8, length=2)
        [m.set_linewidth(.6) for m in ax1.spines.itervalues()]
        #plt.title(r'Dynamic width, SiII 1808')
        plt.plot(liwl, accuaccu/area.sum(), 'black', drawstyle='steps-mid', lw=.8)
        ax1.axhline(y=.5, color='k', lw=.6, ls=':')
        ax1.axvline(x=M1, color='k', lw=.6, ls=':')

        plt.ylabel('Accu. area', fontsize=9)
        plt.axis((pltmin, pltmax, 0, 1))
        plt.axvline(x=llo, ymin=0, ymax=1, color='black', alpha=.5, linestyle=':', lw=.5)
        plt.axvline(x=lhi, ymin=0, ymax=1, color='black', alpha=.5, linestyle=':', lw=.5)
        #plt.figtext(.67, .86, indict.columns.name + '\n' + species, fontsize=9)
        plt.figtext(.67, .86, indict.columns.name + '\n' + 'Si II $\lambda$ 1260', fontsize=9)
        #plt.figtext(.64, .77, species, fontsize=9)
        #gray dots:
        plt.axvline(x=xlo, color='gray', linestyle='--', alpha=.6, lw=.5)
        plt.axvline(x=xhi, color='gray', linestyle='--', alpha=.6, lw=.5)
        ax2 = plt.subplot(gs[1:,:], sharex=ax1)
        ax2.tick_params(labelsize=8, length=2)
        [m.set_linewidth(.6) for m in ax2.spines.itervalues()]
        ax2.yaxis.set_ticks([0, .2, .4, .6, .8, 1.])
        ax2.axvline(x=M1, color='k', lw=.6, ls=':')
        plt.plot(wl, spedata_orig, 'black', drawstyle='steps-mid', lw=.8)
        plt.axhline(y=norm, xmin=0, xmax=1, color='black', linestyle='--', lw=.5)
        plt.axhline(y=0., xmin=0, xmax=1, color='black', linestyle='--', lw=.5)
        #plt.fill_between(liwl, line, norm, color='green', alpha=.4)
        mpt.fill_between_steps(liwl, line, norm, h_align='mid', ax=ax2,
                               color='green', alpha=.4)
        plt.axis((pltmin, pltmax, -.5, 1.38))
        plt.xlabel(r'Velocity [km/s]', size=9)
        plt.ylabel('Relative flux', size=9)  # , labelpad=0)
        # mark limits of selected region
        plt.axvline(x=llo, ymin=0, ymax=1, color='black', alpha=.5, linestyle=':', lw=.5)
        plt.axvline(x=lhi, ymin=0, ymax=1, color='black', alpha=.5, linestyle=':', lw=.5)
        #plt.axvspan(xlo, xhi, facecolor='grey', alpha=.3)
        plt.axvline(x=xlo, color='gray', linestyle='--', alpha=.6, lw=.5)
        plt.axvline(x=xhi, color='gray', linestyle='--', alpha=.6, lw=.5)
        # Bar to mark the velocity width, with caps.
        plt.plot([xlo, xlo], [1.08, 1.12], linestyle='-', color='black', linewidth=.8)
        plt.plot([xhi, xhi], [1.09, 1.11], linestyle='-', color='black', linewidth=.8)
        plt.plot([xlo, xhi], [1.10, 1.10], linestyle='-', color='black', lw=.8)
        # Show v_int
        plt.plot([M1, M1], [-.1, .1], 'k-', lw=0.8)
        plt.annotate(r'$v_{\mathrm{int}}$'+'{0}'.format(''), (M1, -.27),
                     fontsize=10, ha='center', va='center', rotation=90.,
                     backgroundcolor='w')
        # Show v_95%
        plt.plot([xlo, xlo], [-.1, .1], 'k-', lw=0.8)
        plt.annotate(r'$v_{95\%}$'+'{0}'.format(''), (xlo, -.27),
                     fontsize=10, ha='center', va='center', rotation=90.,
                     backgroundcolor='none')
        # Show v_min
        plt.plot([minwl, minwl], [-.1, .1], 'k-', lw=0.8)
        plt.annotate(r'$v_{\mathrm{min}}$'+'{0}'.format(''), (minwl, -.27),
                     fontsize=10, ha='center', va='center', rotation=90.,
                     backgroundcolor='none')
        plt.annotate(r'$W(90\%) ='+"%.3s" % str(DynWidth)+
                    ' \mathrm{km/s}$', ((llo+lhi)/2., 1.16),
                     fontsize=10, ha='center', va='bottom', backgroundcolor='w')
        # Show I/I0_min
        #plt.plot([xlo-180, xlo-80], [datamin, datamin], 'k-', lw=0.8)
        plt.plot([pltmin, pltmin+100], [datamin, datamin], 'k-', lw=0.8)
        plt.annotate(r'$\frac{I}{I_0}_{\mathrm{min}}$', (pltmin+115, datamin),
                     fontsize=10, ha='left', va='center',
                     backgroundcolor='none')
        ax1.tick_params(labelbottom='off')
        # Summarize in box:
        summary = [
            r"\noindent $v_{\mathrm{min}} = "  + str(int(abs(minwl))) +r'~ \mathrm{{ km/s}}$ \\'+
            r"$v_{{\mathrm{int}}} = "  + str(int(abs(M1))) +r'~ \mathrm{ km/s} $\\' +
            r"$v_{{95\%}} = " + str(int(abs(xlo))) +r'~ \mathrm{ km/s}$ \\' +
            r"$\frac{{I}}{{I_0}}_{\mathrm{min}} = " +
            str(sp.around(datamin, 3)) + r"$"]
        plt.figtext(.60, .19, summary[0], fontsize=8,
                    bbox=dict(fc='white', lw=.7,
                              ec='0.4', alpha=.9, boxstyle='round'))
        #plt.figtext(.73, .45, summary[1], )
        #plt.figtext(.73, .4, summary[2], )
        fig2.subplots_adjust(hspace=0., bottom=.13, top=.95)
        if savefig:
            fig2.savefig('DynVel-SiII-{0}.pdf'.format(larsname_ul), facecolor='none')
        fig2.show()
    return {'v_min': minwl, 'v_int': M1, 'v_95': xlo,
            'I/I0_min': datamin, 'I/I0_0': datazero[0], 'delta_v': DynWidth}


