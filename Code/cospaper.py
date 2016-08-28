#! /usr/bin/env python
# -*- coding: utf8 -*-

import scipy as sp
import scipy.constants as con
import scipy.interpolate as ip
import pandas as pd
import pyastrolib.astro as pal
import profiles as pro
import lmfit as lm
import LARS_branding as lb
import brewer2mpl as b2m   # Nice ColorBrewer colors instead of ugly defaults
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams, font_manager, cm
from matplotlib.widgets import SpanSelector, Button
import mpl_tricks as mp
import DynVel as dv
from linelists import COSlines, SDSSlines, MWlines, fikdict, wlsdict, \
        hislines, lislines, SiIIlines, wofflines, colorlist
from uncertainties import ufloat as ufl
from uncertainties import unumpy as unp

#================================================================================
#		Plot details setup.
#================================================================================
# pats = 16  # Plot Axis Title Size
# sfont = 14 # Small font for e.g. legends
csz = 5    # Error bar cap size
# leg_prop = font_manager.FontProperties(size=sfont)
brewername = 'Dark2'
brewernumb = max(
    [int(a) for a in b2m.COLOR_MAPS['Qualitative'][brewername].keys()])
brewercolors = b2m.get_map(
    brewername, 'Qualitative', brewernumb, reverse=False).mpl_colors
paircolors = b2m.get_map('Paired', 'Qualitative', 12).mpl_colors
plt.rc('axes', color_cycle=mp.CB_color_cycle)
# Change the tick fontsize for all subsequent figures
# rc('xtick', labelsize='medium')
# rc('ytick', labelsize='medium')
# rc('lines', linewidth = 2.)
rc('font', family='serif')
# rc('font', size = 14)
#rc('text', usetex=True)
rcParams['axes.formatter.limits'] = (-4,4) # Exponential tick labels outside this log range

def load_spectrum(infile):
    larsno = infile.split('_')[0][-2:]
    thedata = sp.genfromtxt(infile)
    thedata = pd.DataFrame(thedata)
    thedata.columns = ['lambda', 'flam', 'stddev', 'flag']
    thedata = thedata.set_index('lambda', drop=False)
    thedata.columns.name='LARS {0}'.format(larsno)
    return thedata


def rebin(indata, factor=3.):
    '''Rebins the COS data by a given factor '''
    # For convenience:
    wave = indata.index.values.astype('float64')
    flam = indata['flam'].values
    errs = indata['stddev'].values
    newshape = sp.floor(wave.shape[0] / factor)
    newwave = sp.zeros((newshape, 1))
    newflam = sp.zeros_like(newwave)
    newerrs = sp.zeros_like(newwave)
    for newbin in sp.arange(newshape):
        newwave[newbin] = wave[newbin*factor:newbin*factor+factor].mean()
        newflam[newbin] = flam[newbin*factor:newbin*factor+factor].sum()
        newerrs[newbin] = sp.sqrt(
            (errs[newbin*factor:newbin*factor+factor]**2).sum())
    binned_array = sp.hstack([newwave, newflam, newerrs])
    binned_frame = pd.DataFrame(binned_array)
    binned_frame = binned_frame.set_index(0, drop=False)
    binned_frame.index.name = indata.index.name
    binned_frame.columns = ['lambda', 'flam', 'stddev']
    binned_frame.columns.name = indata.columns.name
    return binned_frame


def get_redshift(inname, indata=None):
    ''' Gets redshift of galaxy from the Sloan spectroscopy.
    needs an integer as input.

    Call signature:
    a, b, c = get_redshift(larsno)
    '''
    if inname == 4:
        sloanfile = '../LARS-COS/SDSSspectra/l{:02d}_main.ascii.fitresult.ascii'\
            .format(inname)
    else:
        sloanfile = '../LARS-COS/SDSSspectra/l{:02d}.ascii.fitresult.ascii'\
            .format(inname)
    centroids = {}
    for line in open(sloanfile):
        tmp = line.split()

        if len(tmp) < 1:
            continue

        if '.pos' in tmp[0]:
            ident = tmp[0].split('_')[0][1:]
            row = sp.where(SDSSlines[:, 0] == ident)
            restwl = int(sp.floor(float(SDSSlines[row, 2][0][0])))
            try:
                float(tmp[3])
            except ValueError:
                print 'No valid errors on this line'
                continue
            centroids[tmp[0][4:-4]+' {}'.format(str(restwl))] = \
                    [float(SDSSlines[row, 2][0][0]), float(tmp[1]), float(tmp[3])]

    # Create pandas dataframe of the above.
    centframe = pd.DataFrame.from_dict(centroids).T
    try:
        centframe = centframe.drop(['OII 3727', 'OII 3729'])
    except ValueError:
        pass
    centframe.columns = ['Rest wavelength', 'Best fit', 'stddev']
    centframe['redshift'] = centframe['Best fit'] / centframe['Rest wavelength'] - 1.
    # print centframe
    # Weighted average:
    avg = sp.average(
        centframe['redshift'],
        weights=centframe['stddev']**(-2))
    # Standard deviation:
    stddev = centframe.redshift.std()
    print 'Weighted average: {}'.format(avg)
    print 'Unweighted average: {}'.format(centframe.redshift.mean())
    print "difference: {}".format(abs(avg - centframe.redshift.mean()))
    print 'Std. dev. of redshift: {}'.format(centframe.redshift.std())
    if indata is not None:
        indata.z_value = avg
        indata.z_stdev = centframe.redshift.std()
        return centframe, indata
    else:
        return centframe


def redshift_from_list(indata):
    larsno = indata.columns.name.split('_')[0][-2:]
    larsno = int(larsno)
    # print larsno
    zlist = sp.genfromtxt('redshifts.txt', delimiter='\t')
    zlist = pd.read_csv('redshifts.txt', sep='\t', index_col=0)
    #zlist = pd.DataFrame(zlist).set_index(0)
    zlist.index = zlist.index.astype(int)
    redshift = zlist.ix[larsno].values[0]
    indata.z_value = redshift
    return indata


def mw_unred(data, larsno=None):
    ''' Takes as input a DataFrame as created by the load_spectrum() function.
    If the dataframe doesn't have the dataframe.columns.name attribute set to
    "LARS XX", the LARS galaxy number as an integer is also reuired.
    If the `larsno` keyword is set to something else than None, this value will
    override the one in the dataframe.
    '''
    mw = pd.read_table('LARS-MW.dat', sep=',', index_col=0)
    mwebv = mw['MW E(B-V)']
    # If not explicitly given, set "larsno" to the one given in the dataframe.
    # Raise error if not succesful.
    if larsno is None:
        larsno = data.columns.name.split()[-1]
        try:
            int(larsno)
        except ValueError:
            print 'The `data.columns.name` property is not properly formatted.'
            raise
    ebv = mwebv[int(larsno)]
    print ebv
    data['mw_unred'] = pal.ccm_unred(data['lambda'], data.flam, ebv)
    return data


def _single_gauss_residuals(params, wav, data=None, stddev=None):
    '''To be used by lmfit.
    '''
    ampl = params['ampl'].value
    cent = params['cent'].value
    sigma = params['sigma'].value
    cont = params['cont'].value
    model = cont + pro.gauss(wav, cent, sigma, ampl)
    if data is None:
        return model
    resid = data - model
    if stddev is None:
        return resid
    redresid = resid * stddev ** (-2.)
    return redresid


def fix_errors(indata, zeroval='nan'):
    '''
    zeroval: number or string that is floatable.
        The value to put instead of 0 where this is the value of the error
        column, because division by zero is bad and stupid and ugly.
        Default is 'nan', will be convertet to numpy NaN.
    '''
    zeroix = sp.where(indata.stddev == 0.)
    zeroval = sp.float64(zeroval)
    indata.stddev.iloc[zeroix] = zeroval
    return indata


def wavelength_calibrate(data, larsno=None):
    """ 'Calibrates' wavelengths, so far, by simply measuring LyA-GeoCor by
    maximum flux bin and centroid of a gaussian fit, comparing the two.
    What to do with this knowledge is yet to be decided...
    In fact, no further wavelength calibration is done from the LyaGC line now
    that it has been shown that they don't diverge by more then a shioft to
    (from?) a heliocentric frame."""
    if larsno is None:
        larsno = data.columns.name.split('_')[0][-2:]
    wavl = data['lambda'].values
    flam = data.flam.values
    errs = data.stddev.values
    lyaline = 1215.67
    fitrange = sp.where((wavl > lyaline - 5.) & (wavl < lyaline + 5.))
    top = flam[fitrange].max()
    maxwav = wavl[fitrange][flam[fitrange] == top]
    print maxwav
    print 'Top: ', top
    print "maximum at: {} Aangstrohms".format(maxwav)
    fitguess = pro.gauss(wavl[fitrange], maxwav[0], .3, top)
    # load lmfit Parameters() object
    params = lm.Parameters()
    params.clear()
    params.add('cent', value=lyaline, min=lyaline - 5., max=lyaline + 5.)
    params.add('sigma', value=.4)
    params.add('ampl', value = top)
    params.add('cont', value=sp.median(flam[fitrange]))
    thefit = lm.minimize(
        _single_gauss_residuals,
        params,
        args=(
            wavl[fitrange],
            flam[fitrange],
            errs[fitrange]),
            method='nelder')
    thefit = lm.minimize(
        _single_gauss_residuals,
        params,
        args=(
            wavl[fitrange],
            flam[fitrange],
            errs[fitrange]),
            method='leastsq')
    modelplot = _single_gauss_residuals(
        thefit.params, wavl[fitrange])
    lm.report_errors(params)
    lm.report_errors(thefit.params)
    difference = params['cent'].value - maxwav[0]
    print difference
    print 'Galaxy: LARS {}'.format(larsno)
    print 'Center from max value: {}'.format(maxwav[0])
    print 'Center from best fit: {}'.format(params['cent'].value)
    print 'Difference: {}'.format(difference)
    plt.plot(wavl[fitrange], flam[fitrange], label='data')
    plt.fill_between(wavl[fitrange], errs[fitrange], color='grey', alpha=.6)
    plt.plot(wavl[fitrange], fitguess, label = 'initial model')
    plt.plot(wavl[fitrange], modelplot, label='best fit')
    plt.legend()
    plt.show()
    return


def _linear_resids(params, wave, data, stddev):
    const = params['const'].value
    slope = params['slope'].value
    model = slope * wave + const
    resid = sp.sqrt((data - model)** 2 / stddev ** 2)
    return resid


def _fit_contin(wave, data, errs, galaxy=None, line=None, centroid=None):
    """
    Takes three arrays wave, data and errors, and optional arguments galaxy,
    line and centroid.
    Opens an interactive plot, where the user clicks-and-drags to mark
    regions to include in fit, and fits the selected data to a straight line.
    Returns slope and additive constant.

    call signature:
        const, slope = _fit_contin(wave, data, errs,
                                   galaxy='some string',
                                   line='some string',
                                   centroid=float)

    (More docs to come)
    """

    def _onselect(vmin, vmax):
        ###============================================
        # Set everything between vmin and vmax to false.
        ###============================================
        mask = sp.where((wave > vmin) & (wave < vmax))
        idx[mask] = True
        ax.axvspan(vmin, vmax, color='orange', alpha=.9, zorder=0, picker=True)
        fig.canvas.draw()
        return

    def _onfitbutton(event):
        slope_guess = 0.01
        const_guess = sp.median(data[idx])
        pars = lm.Parameters()
        pars.clear()
        pars.add('const', const_guess)
        pars.add('slope', slope_guess)
        print pars
        result = lm.minimize(
            _linear_resids, pars, args=(wave[idx], data[idx], errs[idx]))
        fitresults['slope'] = pars['slope'].value
        fitresults['const'] = pars['const'].value
        const, slope = pars['const'].value, pars['slope'].value
        ax.plot(wave, wave * slope + const)
        print full_path
        fig.savefig(full_path)
        print('Figure saved as {}'.format(full_path))
        return

    def _on_ok_button(event):
        print fitresults
        plt.close(fig)
        return

    def _on_clr_button(event):
        ax.patches = []
        ax.lines.pop()
        idx[:] = False
        fig.canvas.draw()
        return

    idx = sp.zeros_like(data).astype(bool)

    # Set the initial value of const and slope to NaN's. This will stay the
    # value of these if no fit is made, in which case the data and errors for
    # the normalized lines will also be NaN's.
    fitresults = {'const': sp.nan, 'slope': sp.nan}

    filename = '_'.join([galaxy, line])
    full_path = 'Figures/continuum_fits/{}'.format(filename)

    ###=========================================================================
    #   Make plot.
    ###=========================================================================
    fig = plt.figure(10)
    plt.subplots_adjust(left=.1, right=.89)
    ax = fig.add_subplot(111)
    ax.cla()
    ax.plot(wave, data, 'k', drawstyle='steps-mid', label=line)
    ax.fill_between(wave, errs, color='gray', alpha='.5')
    print sp.median(data)
    ax.set_ylim(bottom=0., top=sp.median(data) * 2.)
    axfit = plt.axes([0.9, 0.84, 0.09, 0.06])
    axclr = plt.axes([0.9, 0.77, 0.09, 0.06])
    ax_ok = plt.axes([0.9, 0.70, 0.09, 0.06])
    fitbu = Button(axfit, 'Fit')
    clrbu = Button(axclr, 'Clear')
    ok_bu = Button(ax_ok, 'OK')
    fitbu.on_clicked(_onfitbutton)
    clrbu.on_clicked(_on_clr_button)
    ok_bu.on_clicked(_on_ok_button)
    if centroid is not None:
        ax.axvline(
            x=centroid, linestyle='--', color='gray', alpha=.6, lw=1.5)
    if line is not None:
        ax.legend(loc='lower left', fancybox=True, shadow=True)
    if galaxy is not None:
        ax.set_title(galaxy)
    s = 'Click and drag to mark ranges to include in fit.'# \n \
    helptext = ax.text(.5, .95,  s,
            bbox=dict(fc='white',
                      ec='0.8', alpha=.9, boxstyle='round'),
            horizontalalignment='center',
            verticalalignment='top',
            transform=ax.transAxes)
    span = SpanSelector(ax, _onselect, 'horizontal',
                        useblit=True, minspan=.5)
    plt.show()
    idx = sp.array([], dtype=int)
    const, slope = fitresults['const'], fitresults['slope']
    return const , slope


def normalize_line(indata, line, overwrite=False, smooth=True, illust=False):
    ''' line: string, of the form 'At ION xxxx', must be in the COSlines table.
    '''
    print indata.z_value, 'normalize_line() 1'
    if not line in wlsdict:
        print 'Line not found.'
        return indata
    if '_'.join(line.split()) in indata.columns:
        if not overwrite:
            print(
                'Normalised version of line {} already present in the dataset {}'
                .format(line, indata.columns.name))
            return indata
    # Find the proper restwl from the line name
    rwl = sp.float64(line.split()[-1])
    matchidx = sp.where(
        sp.absolute(COSlines[:, -1].astype(sp.float64) - rwl) < 1.)
    labwl = COSlines[matchidx, -1].astype(sp.float64)[0][0]
    # If redshift is not attached, do it. Pickle can not save this, so depends
    # on how the dataframe is created.
    print indata.z_value, 'normalize_line()'
    # import ipdb; ipdb.set_trace()  # XXX BREAKPOINT
    try:
        z = indata.z_value
        print 'z was defined!!'
    except AttributeError:
        indata = redshift_from_list(indata)
        z = indata.z_value
    obswl = labwl * (1. + z)
    # Slice the part of the dataframe that lies within the desired ploitrange.
    # Create numpy arrays for these, to send to fitting function.
    plotrange = sp.where((indata['lambda'] > obswl - 20.) &
                         (indata['lambda'] < obswl + 20.))
    wave = indata['lambda'].values[plotrange]
    data = indata['flam'].values[plotrange]
    plotdata = indata['flam'].values[plotrange]
    if smooth:
        try:
            plotdata = indata['smooth'].values[plotrange]
        except:
            print 'No smoothed array exists.'
            pass
    errs = indata['stddev'].values[plotrange]
    # Do the actuall fitting set line in dataframe, return dataframe.
    const, slope = _fit_contin(
        wave, plotdata, errs,
        galaxy=indata.columns.name,
        line=line,
        centroid=obswl)
    # TODO: Maybe save slope, const values for documentation? In a file?
    linear = wave * slope + const
    data /= linear
    if illust:
        plt.figure()
        plt.plot(wave, data)#plot
        plt.show()
    errs /= linear
    linefmt = '_'.join(line.split())
    indata.loc[:, linefmt] = sp.nan
    indata.loc[:, '{}_Stddev'.format(linefmt)] = sp.nan
    indata.loc[:, linefmt].iloc[plotrange] = data
    indata.loc[:, '{}_Stddev'.format(linefmt)].iloc[plotrange] = errs
    return indata


def normalize_lines(indata, lineslist, overwrite=False, smooth=True):
    ''' Wrapper for normalize_line(), basically just iterates through a list
    of lines in the format "At ION xxxxx" and sends it through normalize_line()
    '''
    print indata.z_value, 'normalize_lines'
    # import ipdb; ipdb.set_trace()  # XXX BREAKPOINT
    for i in lineslist:
        redshift = indata.z_value
        indata = indata.dropna(subset=['stddev'])
        indata.z_value = redshift
        indata = normalize_line(
            indata, #.dropna(subset=['stddev']),
            i,
            overwrite=overwrite,
            smooth=smooth)
    return indata


def make_line_profile_plots(
    indict=None, normalize=True, smooth=5, lines='default', ionization=False):
    '''
    lines: 'default', 'his', 'lis', any line in COSlines array.
        If not 'default', 'his' or 'lis', must be a dict of the form:
            {'At Io xxxx': xxxx.xx},
        with At and Io being the atom name and the ionization degree,
        respectively, like e.g. {'C II 1334': 1334.53}
    '''
    if lines == 'default':
        lines = {'O I 1302': 1302.17,
                 'Si II 1190': 1190.42,
                 'C II 1334': 1334.53,
                 'Si III 1206': 1206.50}
    elif lines == 'his':
        lines = {'Si IV 1122': 1122.49,
                 #'Si III 1206': 1206.50,
                 'Si IV 1393': 1393.76,
                 'Si IV 1402': 1402.77}
    elif lines == 'lis':
        lines = {'Si II 1190': 1190.42,
                 'Si II 1193': 1193.29,
                 'Si II 1260': 1260.42,
                 'Si II 1304': 1304.37,
                 'Si II 1526': 1526.72}
    elif lines =='hilo':
        lines = {'Si II 1190': 1190.42,
                 'Si II 1193': 1193.29,
                 'Si II 1260': 1260.42,
                 'Si II 1304': 1304.37,
                 'Si II 1526': 1526.72,
                 'Si III 1206': 1206.50,
                 'Si IV 1122': 1122.49,
                 'Si IV 1393': 1393.76,
                 'Si IV 1402': 1402.77}

    elif type(lines) == dict:  # Assume it is a dict with correct format.
        lines = lines
    else:  # Assume it is just a string with the line name in correct format.
        labwl = sp.float64(lines.split()[-1])#.astype(sp.float64)
        matchidx = sp.where(COSlines[:, -1].astype(sp.float64) - labwl < 1.)
        linedict = {}
        linedict[lines] = COSlines[matchidx, -1].astype(sp.float64)
        lines = linedict
    if ionization:
        lines = {'Si IV 1122': 1122.49,
                 'Si IV 1393': 1393.76,
                 'Si IV 1402': 1402.77,
                 'Si II 1190': 1190.42,
                 'Si II 1193': 1193.29,
                 'Si II 1260': 1260.42,
                 'Si II 1304': 1304.37,
                 'Si II 1526': 1526.72}

    dfdict = {}
    # Make sure colors are consistent between subplots:
    # colors = colorlist[:len(lines.keys())]
    # colors = brewercolors[:len(lines.keys())]
    # colors = b2m.get_map(
    #     brewername, 'Qualitative', len(lines.keys())).mpl_colors
    #colors = mp.CB_color_cycle[0:]
    colors = ['green', 'orange', 'purple',
              'gray', 'blue', 'magenta', 'cyan', 'maroon']
    #colors = mp.CUDCycle
    # colors = b2m.get_map(
    # brewername, 'Qualitative', brewernumb, reverse=False).mpl_colors

    # Do this as the last thing, so we know which lines are actually present!
    colord = dict(zip(lines, colors))

    # Before creating the first figure, make sure all needed line are present
    # normalized and masked.
    if indict is None:
        inlist = []
        iterrange = sp.arange(14)
    else:
        iterrange = sp.array([x - 1 for x in indict.keys()])

    print iterrange  # XXX Debug.

    for j in iterrange:  # sp.arange(14):
        i = j + 1
        if indict is None:
            df = prepare_dataframe(i)
        else:
            df = indict[i]
            try:
                z = df.z_value
            except:
                df = redshift_from_list(df)
                z = df.z_value
        # Normalize:
        print df.z_value
        # import ipdb; ipdb.set_trace()  # XXX BREAKPOINT
        df = normalize_lines(df, lines.keys())  # YAY sane defaults!
        if indict is None:
            inlist.append(df)
        else:
            indict[i] = df
    if indict is None: indict = dict(zip(range(1, 15), inlist))
    for j in iterrange:  # sp.arange(14):  # FIXME: Det her er jo åndssvagt!
        i = j + 1
        # Mask:
        df = indict[i]
        df = set_mask(df, lines.keys(), skip=True)
        indict[i] = df
    # Create figure:
    fig = plt.figure(1, figsize=(12, 10))
    fig.subplots_adjust(hspace=.04, wspace=.02)
    fig.show()
    subplots = 0
    #Now, cycle through the galaxies and do what must be done for each.
    for j in iterrange:  # sp.arange(14):
        i = j + 1
        gal = 'Lars{:02d}'.format(i)
        subplots += 1
        # Make possible to send in dictionary of frames with existing lines so
        # not aving to renormalize lines every single time:
        df = indict[i]
        ax = fig.add_subplot(5, 3, subplots)
        if not sp.mod(subplots, 3).astype(bool):
            ax.tick_params(labelleft='off', labelright='on')
        if not sp.mod(subplots + 1, 3).astype(bool):
            ax.tick_params(labelleft='off', labelright='off')
        # Now make sure to only try to handle lines that are actually there, so
        # the thing doesn't die from one line that didn't work (like in e.g.
        # Lars 11).
        lineslist = []
        for line in lines.keys():
            print 'Going with the following lines {} for {}: '.format(
                lines.keys(), df.columns.names)
            print 'Dataframe columns: ', df.columns
            print line
            if not line.replace(' ', '_') in df.columns:
                continue
            if df[line.replace(' ', '_')].notnull().any():
                lineslist.append(line)
                print 'Added line ', line
        #import pdb; pdb.set_trace()  # XXX BREAKPOINT
        # Now, regrid in velocity space (since the spectra have equal
        # wavelength bins and hence non-equal velocity bins).
        #print df.z_value
        # import ipdb; ipdb.set_trace()  # XXX BREAKPOINT
        print lineslist, 'before veloc regrid'  # DEBUG
        vf = velocity_regrid(df, lineslist=lineslist)
        print 'Velocity regridding done for {}'.format(df.columns.name)
        #print vf
        notempty = True
        if len(vf.columns) == 0:
            notempty = False
        # FIXME: THIS IS BAD CODE. When having two essentially different
        # operations in the same function, parted by an if statement, I need to
        # rethink a thing or two. But it works.
        print lineslist, 'after veloc regrid'  # DEBUG
        if notempty:
            vels = vf.Velocity
            # Take averages if ionization comparison
            if ionization:
                hidata = []
                hierrs = []
                himask = []
                lodata = []
                loerrs = []
                lomask = []
                for col in vf.columns:
                    print col.split('_')
                    if not vf[col].notnull().any():
                        continue
                    if ' '.join(col.split('_')[:-1]) in lislines:
                        if (col.endswith('_Flam')):
                            lodata.append(col)
                        elif col.endswith('_Stddev'):
                            loerrs.append(col)
                        elif col.endswith('_Mask'):
                            lomask.append(col)
                    elif col.startswith('Si_IV'):
                        if col.endswith('_Flam'):
                            hidata.append(col)
                        elif col.endswith('_Stddev'):
                            hierrs.append(col)
                        elif col.endswith('_Mask'):
                            himask.append(col)
                if len(himask) > 0:
                    himaskarr = vf[himask].values
                    himaskarr = sp.invert(himaskarr)
                    hidataarr = vf[hidata].values
                    hierrsarr = vf[hierrs].values
                    himamean = sp.ma.masked_array(
                        hidataarr, mask=himaskarr).mean(axis=1)
                    hiermean = sp.ma.masked_array(
                        hierrsarr, mask=himaskarr)
                    hiermean = sp.sqrt((hiermean**2).sum(axis=1))
                    vf['Hi_Ion_Mean'] = himamean
                    vf['Hi_Ion_smooth'] = sp.convolve(
                        himamean, sp.ones(5)/5., mode='same')
                    vf['Hi_Ion_Mean_Stddev'] = hiermean
                    vf['Hi_Ion_Mask'] = himaskarr.all(axis=1)
                    hilist = ['Hi_Ion_Mean', 'Hi_Ion_smooth',
                              'Hi_Ion_Mean_Stddev', 'Hi_Ion_Mask']
                    maskedarr = sp.ma.masked_array(
                        vf.Hi_Ion_smooth.values,
                        vf.Hi_Ion_Mask.values)
                    complemnt = sp.ma.masked_array(
                        vf.Hi_Ion_smooth.values,
                        sp.invert(vf.Hi_Ion_Mask.values))
                    ax.plot(vf.Velocity.values,
                            maskedarr, ls='-',
                            color=mp.CB_color_cycle[0],
                            label='High ionization average')
                    ax.plot(vf.Velocity.values,
                            complemnt, ls='--',
                            color=mp.CB_color_cycle[0],)
                else:
                    hilist = []

                if len(lomask) > 0:
                    lomaskarr = vf[lomask].values.astype(bool)
                    lomaskarr = sp.invert(lomaskarr)
                    lodataarr = vf[lodata].values
                    loerrsarr = vf[loerrs].values
                    lomamean = sp.ma.masked_array(
                        lodataarr, mask=lomaskarr).mean(axis=1)
                    loermean = sp.ma.masked_array(
                        loerrsarr, mask=lomaskarr)
                    loermean = sp.sqrt((loermean**2).sum(axis=1))
                    vf['Lo_Ion_Mean'] = lomamean
                    vf['Lo_Ion_smooth'] = sp.convolve(
                        lomamean, sp.ones(5)/5., mode='same')
                    vf['Lo_Ion_Mean_Stddev'] = loermean
                    vf['lo_Ion_Mask'] = lomaskarr.all(axis=1)
                    lolist = ['Lo_Ion_Mean', 'Lo_Ion_smooth',
                              'Lo_Ion_Mean_Stddev', 'lo_Ion_Mask']
                    maskedarr = sp.ma.masked_array(
                        vf.Lo_Ion_smooth.values,
                        vf.lo_Ion_Mask.values)
                    complemnt = sp.ma.masked_array(
                        vf.Lo_Ion_smooth.values,
                        sp.invert(vf.lo_Ion_Mask.values))
                    ax.plot(vf.Velocity.values,
                            maskedarr, ls='-',
                            color=mp.CB_color_cycle[1],
                            label='Low ionization average')
                    ax.plot(vf.Velocity.values,
                            complemnt, ls='--',
                            color=mp.CB_color_cycle[1],)
                else:
                    lolist = []
                # TODO: Is a masked array with masked pyt-summed errors
                print 'Masked arrays success'
                lineslist = lolist + hilist
                vf = vf[lineslist+['Velocity']]

            if not ionization:
                print lineslist, 'line 709' # DEBUG
                for k in lineslist:
                    k = k.replace(' ', '_')
                    print k
                    data = vf[k+'_Flam']
                    if k+'_Mask' in vf.columns:
                        mask = vf[k+'_Mask'].values
                    else:
                        mask = sp.ones_like(data).astype(bool).values
                    # Smooth, for clarity of the view only.
                    if smooth is not None:
                        data = sp.convolve(
                            data, sp.ones(smooth)/float(smooth),
                            mode='same')
                    maskedarr = sp.ma.masked_array(data, sp.invert(mask))
                    complemnt = sp.ma.masked_array(data, mask)
                    print 'colord.keys', colord.keys()
                    #import pdb; pdb.set_trace()  # XXX BREAKPOINT
                    ax.plot(vels, maskedarr, ls='-',  drawstyle='steps-mid',
                            color=colord[k.replace('_', ' ')], lw=1.2)
                    ax.plot(vels, complemnt, ls='--', drawstyle='steps-mid',
                            color=colord[k.replace('_', ' ')], lw=1.2, alpha=.7)
                    if j == max(iterrange):  #14:
                        ax.plot(vels, data+20.,
                                drawstyle='steps',
                                color=colord[k.replace('_', ' ')],
                                label = k.replace('_', ' '))
                        fig.canvas.draw_idle()
        ax.set_ylim(bottom=-.2, top=2.)
        ax.set_xlim((-900., 900.))
        ax.axvline(x=0., color='k', linestyle='--', linewidth=1.)
        ax.axhline(y=1., color='k', linestyle='--', linewidth=1.)
        ax.axhline(y=0., color='k', linestyle='-', linewidth=1.)
        ax.annotate('LARS {}'.format(i), (.05, .8), xycoords='axes fraction')
        print 'LARS {}'.format(i)
        if subplots < 12:
            ax.tick_params(labelbottom='off')
        if subplots == 7:
            ax.set_ylabel('Relative Flux', fontsize='large')
        print max(iterrange), 'LARSNO'
        if subplots == max(iterrange) + 1:  # 14:
            ax.set_xlabel('Velocity [km/s]', fontsize='large')
            leg = ax.legend(
                frameon=False,
                fontsize='large',
                loc='center',
                bbox_to_anchor=(1.5, 0.5),
                labelspacing=.2)
            leg.draggable()
        fig.subplots_adjust(left=.07, right=.95, top=.95, )
        fig.canvas.draw()
        plt.draw()

        # Append this datadict (df, NOT vf!!) to dict, for later saving.
        dfdict[i] = df
    return dfdict


def show_spec(lnum=None, inframe=None, smooth=None,
              color='black', show_lines=True):
    '''Docstring goes here'''

    if inframe is None:
        datafr = prepare_dataframe(lnum)
        larsname = datafr.columns.name
    else:
        datafr = inframe
        larsname = datafr.columns.name.replace('_', ' ')
    toplot = 'flam'
    if smooth is not None:
        datafr['smooth'] = sp.convolve(
            datafr.flam.values, sp.ones(5) / 5., mode='same')
        toplot = 'smooth'
    datafr.plot('lambda', toplot, color=color,
                drawstyle='steps', label='LARS {}'.format(lnum))
    datafr.plot('lambda', 'stddev',
                color='gray', lw=2, alpha=.6, ax=plt.gca())
    if show_lines is True:
        for line in MWlines.keys():
            plt.axvline(x=MWlines[line], linestyle='--',
                        color='darkorange', lw=.8)
            plt.annotate(line+' MW', (MWlines[line], 0.9),
                        xycoords=('data', 'axes fraction'),
                        color='darkorange',
                        rotation=90,
                        backgroundcolor='w')
            plt.axvline(x=MWlines[line]*(1+datafr.z_value),
                        linestyle='-.', color='teal', lw=.8, zorder=2)
            plt.annotate(line+larsname, (MWlines[line]*(1+datafr.z_value), 0.9),
                        xycoords=('data', 'axes fraction', ),#zorder=2),
                        color='teal',
                        rotation=90,
                        backgroundcolor='w',
                        zorder=2)
    yval = sp.median(datafr[toplot]) * 2.6
    print sp.median(datafr[toplot].dropna())
    plt.ylim(-.1*yval, yval)
    plt.xlim(datafr['lambda'].dropna().min(), datafr['lambda'].dropna().max())
    plt.show()
    return datafr


def prepare_dataframe(larsno, smooth=5):
    """ Runs the steps most commonly needed. """
    filestr = 'COSspectra1.0/l{:02d}_cos.coarse.ascii'.format(larsno)
    datafr = load_spectrum(filestr)
    datafr = rebin(datafr, 6)
    datafr = fix_errors(datafr)
    datafr = redshift_from_list(datafr)  # get_redshift(larsno, datafr)
    datafr = mw_unred(datafr)
    datafr['smooth'] = sp.convolve(
        datafr.flam,
        sp.ones(smooth) / smooth,
        mode='same')
    return datafr


def load_dfdict(infile='dfdict.pickle'):
    ''' Because pandas cannot save metadata, we need to re-calculate redshifts
    and attach them again.
    '''
    import pickle
    thefile = open(infile, 'read')
    dfdict = pickle.load(thefile)
    thefile.close()
    for i in dfdict.keys():
        dfdict[i] = redshift_from_list(dfdict[i])
        # print 'Loaded {0}'.format(i)
    return dfdict


def velocity_regrid(df, lineslist=None, degree=3, illustrate=False, maskit=True):
    ''' Reads a list of lines which must all be in the dataframe, takes the
    data from these lines and runs an interpolation on them so that they fall
    into the same velocity bins, so that later, we can create the covering
    fraction map that we so eagerly yearn for.
    '''
    ###============================================================
    # First, because I can't be botheres passing a list every time:
    ###============================================================
    if lineslist is None:
        lineslist = ['Si_II_1190', 'Si_III_1206', 'O_I_1302', 'C_II_1334']
    if len(lineslist) == 0:
        return pd.DataFrame()
    if lineslist == 'all':
        lineslist = []
        for line in df.columns:
            if ' '.join(line.split('_')) in wlsdict.keys():
                lineslist.append(line)
        #print lineslist
    elif len(lineslist[0].split('_')) < 2:
        ll = []
        for line in lineslist:
            ll.append(line.replace(' ', '_'))
        lineslist = ll
    print lineslist
    # Then, just because it's convenient:
    try:
        larsnr = int(df.columns.name.split()[-1])
    except ValueError:
        larsnr = 42
    ###==================================================
    # Then: We don't want NaN's in errors. And we're not overwriting anything,
    # anyway, cause we're not returning the input dataframe and we're working
    # on a copy.
    ###==================================================
    redshift = df.z_value            # Redshift preserved.
    df = df.dropna(subset=['stddev'])
    df.z_value = redshift
    # df = redshift_from_list(df)
    ###==================================================
    # Now: Make a velocity vector from the wavelengths of each transition
    # data chunk. Also, plot to show how it's going.
    ###==================================================
    velsdict = {}
    if illustrate:
        # plt.ion()
        fig = plt.figure('velvects')
        ax = fig.add_subplot(111)
    ###==================================================
    # Check if masking arrays exist for all lines:
    ###==================================================
    hasmask = []
    isthere = []
    notnull = []
    print lineslist
    #import pdb; pdb.set_trace()  # XXX BREAKPOINT
    for line in lineslist:
        if not df[line].notnull().any():
            print('No data available for line {}'.format(line))
            continue
        if not line+'_Mask' in df.columns and maskit:
            df = set_mask(df, line)
        hasmask.append(line+'_Mask' in df.columns)
        isthere.append(line in df.columns)
    if not sp.array(isthere).all():
        print 'Not all lines in the line list are present in input data.'
        print 'Quitting.'
        return
    if sp.array(hasmask).any() and not sp.array(hasmask).all():
        print 'Some, but not all, lines have masks.'
        print 'Proceeding without masks.'
    hasmask = sp.array(hasmask).all()

    for line in lineslist:
        if not df[line].notnull().any():
            print 'No non-null data for the line {}'.format(line)
            continue
        ###==================================================
        # First: Find the actual rest wavelength from the line name.
        ###==================================================
        restwl = float(line.split('_')[-1])
        matchidx = sp.where(
            sp.absolute(COSlines[:, -1].astype(sp.float64) - restwl) < 1.)
        labwl = COSlines[matchidx, -1][0][0].astype(sp.float64)
        print labwl
        # import pdb; pdb.set_trace()  # XXX BREAKPOINT
        # Lab wavelength to observed wavelength:
        obswl = labwl * (1. + df.z_value)
        ###==================================================
        # Next: turn wavelengths into velocities by aid of this, and slice the
        # data and errors accordingly.
        ###==================================================
        wavlens = df[line].dropna().index.values
        errs = df[line+'_Stddev'].ix[df[line].dropna().index]
        data = df[line].ix[df[line].dropna().index]
        vels = (wavlens - obswl) / obswl * con.c / 1000.
        vels = sp.vstack([vels, data, errs]).T
        ###==================================================
        # Interpolate Masks, if they exist:
        ###==================================================
        if hasmask:
            mask = df[line+'_Mask'].ix[df[line].dropna().index]
            #print mask.shape, vels.shape
            vels = sp.hstack([vels, mask.reshape((-1, 1))])
        ###==================================================
        # Put this into dict of velocity vectors
        ###==================================================
        velsdict[line] = vels
        if illustrate:
            ax.plot(vels[:, 0], label=' '.join(line.split('_')))
    if illustrate:
        ax.legend()
        plt.show()
    ###==================================================
    #Now, find the line with the coarsest velocity binning.
    ###==================================================
    binmax = 0.
    velmin = -1.e9
    velmax = 1.e9
    binline = ''
    print df.columns.name
    for line in velsdict.keys():
        binning = velsdict[line][:, 0].ptp()
        vmin = velsdict[line][:, 0].min()
        vmax = velsdict[line][:, 0].max()
        if binning > binmax:
            binmax = binning
            binline = line
        if vmax < velmax:
            velmax = vmax
        if vmin > velmin:
            velmin = vmin
    ###==================================================
    # Now cut off the coarsest-gridded one so its limits are those of the
    # narrowest limited of the lines. Do the same for its values, so we can
    # still interpolate and/or plot etc.
    ###==================================================
    wavidx = ((velsdict[binline][:, 0] > velmin) &
              (velsdict[binline][:, 0] < velmax))
    velsdict[binline] = velsdict[binline][wavidx, :].astype(sp.float64)
    linedf = pd.DataFrame(velsdict[binline])
    if hasmask:
        linedf.columns = [
            'Velocity', '{}_Flam'.format(binline),
            '{}_Stddev'.format(binline), '{}_Mask'.format(binline)]
        linedf[binline+'_Mask'] = linedf[binline+'_Mask'].astype(bool)
    else:
        linedf.columns = [
            'Velocity', '{}_Flam'.format(binline),
            '{}_Stddev'.format(binline)]
    linedf.columns.name = df.columns.name
    ###==================================================
    # Cycle through other lines, interpolate them and evaluate in the chosen
    # velocity bins. These will be our new data.
    ###==================================================
    for line in velsdict.keys():
        if line == binline:
            continue
        wave = velsdict[line][:, 0]
        flam = velsdict[line][:, 1]
        errs = velsdict[line][:, 2]
        splam = ip.splev(linedf.Velocity.values, ip.splrep(wave, flam))
        spler = ip.splev(linedf.Velocity.values, ip.splrep(wave, errs))
        linedf[line+'_Flam'] = splam
        linedf[line+'_Stddev'] = spler / sp.sqrt(2.)
        if hasmask:
            mask = velsdict[line][:, 3]
            splask = ip.splev(
                linedf.Velocity.values, ip.splrep(wave, mask, k=1))
            linedf[line+'_Mask'] = sp.rint(splask).astype(bool)

    #print linedf.columns
    linedf.z_value = redshift
    return linedf


def _coverfrac_resid(params, wave, fik, data=None, errs=None, diagnostics=False):
    '''Calculates residuals to minimize for
    covering fraction and column density.
    '''
    fc = params['fc'].value
    colden = params['colden'].value
    tau = -1. * fik * wave * colden / 3.768e14
    factorx = 1. - sp.exp(tau)
    rlflx = 1. - fc * (1. - sp.exp(tau))
    ###==================================================
    # If no data input, return model calculated form input params, good for
    # generating plots etc.:
    ###==================================================
    if diagnostics:
        xvals = fiks * wave
        x = sp.arange(xvals.min()-5., xvals.max()+5., 1000)
        taus = -1. * xvals * colden / 3.768e14
    if data is None:
        return rlflx
    elif errs is None:
        return rlflx - data
    else:
        resids = (rlflx - data) / errs
        return resids


def cover_fract(indict, species='Si_II', illustr=True, minimizer='lmfit'):
    ''' This takes as input a dataframe of lines regridded in velocity space
    (output of the velocity_regrid() function), creates a velocity/fc map from
    this.
    '''
    fiks = []
    vels = []
    lines = []
    flamcols = []
    errscols = []
    masks = []
    covfrac = []
    covfrer = []
    ### ========================================================================
    #  Select the lines to include in coverfraction calculations
    ### ========================================================================
    for line in indict.columns:
        if line.startswith(species):
            if indict[line].notnull().any():
                if line.split('_')[-1] == 'Flam':
                    flamcols.append(line)
                    lines.append(' '.join(line.split('_')[:-1]))
                elif line.split('_')[-1] == 'Stddev':
                    errscols.append(line)
                elif line.split('_')[-1] == 'Mask':
                    masks.append(line)
    for line in lines:
        fiks.append(fikdict[line])
        vels.append(wlsdict[line])
    fiks = sp.array(fiks)
    vels = sp.array(vels)
    ###==================================================
    # Lists for fc and coldens and their errors
    ###==================================================
    fclst = []
    cdlst = []
    cslst = []
    numlines = []
    brutedic = {}
    #print indict
        # plt.figure()
        # indict.plot(x='Velocity', y='Si_II_1260_Flam')
        # plt.show()
    #import ipdb; ipdb.set_trace()  # XXX BREAKPOINT
    fitvels = indict[sp.absolute(indict.Velocity) < 1000.].index
    # print indict.index
    # print indict.columns
    # print indict.Velocity

    print "Length of fitvels: {}".format(len(fitvels))
    # print fitvels
    for v in fitvels:
        pars = lm.Parameters()
        pars.clear()
        pars.add('colden', 1.e10, min=1e8, max=1e22)
        pars.add('fc', .5, min=0., max=1.)
        data = indict[flamcols].ix[v].values
        errs = indict[errscols].ix[v].values
        mask = indict[masks].ix[v].values
        ###=========================
        # Only include the data/errrors that are not masked out
        ###=========================
        data = data[mask]
        errs = errs[mask]
        numlines.append(len(data))  # FIXME: May want to ditch this.

        ### ==================================================
        #     Brute force show the chi^2 over parameter space:
        #     FIXME: this is DIRTY, DIRTY coding!
        ### ==================================================
        if minimizer == 'brute':
            covfrac = sp.linspace(0., 1., 101)
            logN = sp.linspace(10., 20., 01)
            chisqarr = sp.ones((len(covfrac), len(logN)))
            CF, NC = sp.meshgrid(covfrac, logN)
            # With only one data point: Don't even try!
            if len(data) < 2:
                chisq = sp.ones_like(CF) * 1000.
            else:
                Fik = fiks[mask].reshape(1, 1, -1)
                Wav = vels[mask].reshape(1, 1, -1)
                CF = sp.dstack([CF] * Fik.shape[2])
                NC = sp.dstack([NC] * Fik.shape[2])
                NC = 10. ** NC
                tau = -1. * Fik * Wav * NC / 3.768e14
                factorx = 1. - sp.exp(tau)
                rlflx = 1. - CF * (1. - sp.exp(tau))
                resid = (rlflx - data)**2 / errs**2
                chisq = (resid ** 2).sum(axis=2) / len(data)
            brutedic[indict.Velocity.ix[v]] = chisq
            if not 'covfrac' in brutedic.keys():
                brutedic['covfrac'] = covfrac
            if not 'logN' in brutedic.keys():
                brutedic['logN'] = logN

        if minimizer == 'lmfit':
            #print "Creating fclist, cslist with lmfit option"  # DEBUG XXX
            ###=========================
            # Don't try to solve this with not enough data:
            ###=========================
            print len(data)
            if len(data) < 2:
                #print 'Not enough data for this velocity bin'
                fclst.append([sp.nan, sp.nan])
                cdlst.append([sp.nan, sp.nan])
                cslst.append(sp.nan)
                continue
            ###===================================================
            # Can't do reduced chisquares with zero DoF:
            ###===================================================
            if len(data) == 2:
                result = lm.minimize(
                    _coverfrac_resid, pars,
                    args=(vels[mask], fiks[mask]),
                    kws={'data':data, 'errs':errs},
                    method='anneal')
                fclst.append([pars['fc'].value, 10.])
                cdlst.append([pars['fc'].value, 1e20])
                cslst.append(sp.nan)
            ###=========================
            # With more than 2 data points, we're peachy
            ###=========================
            if len(data) > 2:
                result = lm.minimize(
                    _coverfrac_resid, pars,
                    args=(vels[mask], fiks[mask]),
                    kws={'data':data, 'errs':errs},
                    method='anneal')
                result = lm.minimize(
                    _coverfrac_resid, pars,
                    args=(vels[mask], fiks[mask]),
                    kws={'data':data, 'errs':errs},
                    method='leastsq')
                fclst.append([pars['fc'].value, pars['fc'].stderr])
                cdlst.append([pars['colden'].value, pars['colden'].stderr])
                cslst.append(result.chisqr)
    ###===========================================================
    # Create arrays of the lists of fc and colden values and errs.
    ###===========================================================

    print "Now working on: {}".format(indict.columns.name)  # DEBUG
    print 'Minimizer: {}'.format(minimizer)
    if minimizer == 'brute':
        print "Created brutedict for {}".format(indict.columns.name)
        fclst, cdlst, cslst = handle_brute_dict(brutedic)
    fcarr = sp.array(fclst)
    cdarr = sp.array(cdlst)
    # print "fclist: {}, cdlist: {}".format(fclst, cdlst)  # DEBUG XXX
    # print 'Minimizer:', minimizer
    # print fcarr.shape, cdarr.shape, 'before'
    columns = ['CoverFrac', 'CoverFerr', 'ColDens', 'ColDerr',]
    if minimizer == 'brute':
        fcarr = fcarr[:, 1:]
        cdarr = cdarr[:, 1:]
        cslst = cslst[:, 1]
        columns = ['CoverFrac', 'CoverFerr', 'CoverFLoerr',
                   'ColDens', 'ColDerr', 'ColDLoerr']
    if minimizer == 'lmfit':
        cslst = sp.array(cslst)
    outdf = pd.DataFrame(sp.hstack((fcarr, cdarr)))
    #print outdf.columns, columns
    outdf.columns = columns
    outdf['ChiSquare'] = cslst#[:, 1]
    outdf['Velocity'] = indict.ix[fitvels]['Velocity'].values#.astype('str')
    outdf.columns.name = indict.columns.name
    # Add the new results to the original input vel-space DataFrame.
    #print indict
    indict.loc[:, species+'_CovFrac_map'] = sp.nan  # Get correct shape
    indict.loc[:, species+'_CovFrac_err'] = sp.nan
    indict.loc[fitvels, species+'_CovFrac_map'] = fcarr[:, 0]
    indict.loc[fitvels, species+'_CovFrac_err'] = fcarr[:, 1]
    if minimizer == 'brute':
        indict.loc[:, species+'_CovFrac_lo_err'] = sp.nan
        indict.loc[fitvels, species+'_CovFrac_lo_err'] = fcarr[:, 2]
    indict.loc[:, species+'_ColDens_map'] = sp.nan  # Same as above
    indict.loc[:, species+'_ColDens_err'] = sp.nan
    indict.loc[fitvels, species+'_ColDens_map'] = cdarr[:, 0]
    indict.loc[fitvels, species+'_ColDens_err'] = cdarr[:, 1]
    if minimizer == 'brute':
        indict.loc[:, species+'_ColDens_lo_err'] = sp.nan
        indict.loc[fitvels, species+'_ColDens_lo_err'] = cdarr[:, 2]
    indict.loc[:, species+'_CovFrac_Xsq'] = sp.nan
    indict.loc[fitvels, species+'_CovFrac_Xsq'] = cslst
    indict.loc[:, species+'_NumLines'] = sp.nan
    indict.loc[fitvels, species+'_NumLines'] = numlines

    ###==========================================================
    # Create a figure to show all the newly generated nicenesses.
    ###==========================================================
    if illustr:
        plt.figure()
        plt.subplot(111)
        plt.plot(indict.Velocity[fitvels], 1.-fcarr[:, 0], color='b', drawstyle='steps-mid')
        plt.plot(
            indict.Velocity[fitvels], fcarr[:, 1], lw=3.,
            color='b', alpha=.3)
        plt.axis((-1000, 1000, -.2, 1.8))
        plt.annotate(indict.columns.name, (.1, .9), xycoords='axes fraction')
        plt.show()
    if minimizer == 'brute':
        return indict, brutedic
    else:
        return indict


def make_cover_frac_plots(
    indict, lines=SiIIlines, compare_lines=lislines,
    minimizer='lmfit', inveldf=None, return_brutedict=False):
    ''' Takes as input a dictionary of dataframes - not vel frames but the
    full data frames - and for each, calculate the covering fractions, and for
    each, create a subplot showing covering fractions and, for comparison, an
    average of the LIS line profiles.

    Parameters:
    ===========
    indict : dict
        Dictionary made from the make_lines_profile_plots() function and the
        functions called therein.
    lines : list
        List of lines to include in calculations of covering fraction.
        Typically, these will be of the same species. Line name format is
        ['At Io xxxx', ...]
    compare_lines : list
       List of lines to include in the average line profile used for comparison
       to the derived covering fractions. Same format as with `lines`
    minimizer : str
        name of the minimization method to use for calculating covering
        fractions. At the moment, 'lmfit' and 'brute' are supported.
    inveldf : pandas.DataFrame
        Not actually supported at the moment, might be completely removed from
        the code soon.
    '''
    fig = plt.figure()
    species = '_'.join(lines[1].split()[:2])
    cntnr = fig.add_subplot(111)  # cntnr := container axes.
    cntnr.set_frame_on(False)
    cntnr.set_xticks([])
    cntnr.set_yticks([])
    cntnr.set_xticklabels([])
    cntnr.set_yticklabels([])
    cntnr.set_ylabel('Flux / Continuum', labelpad=40)
    cntnr.set_xlabel('Velocity [km/s]', labelpad=20)
    veldfdict = {}
    ultilist = []
    if inveldf is not None:
        veldf = inveldf.copy()

    ultilist = list(set(lines).union(set(compare_lines)))

    for larsno in indict.keys():
        ###====================================================================
        # Load the data, check if the requested lines are actually present in
        # the dataset, and perform velocity regridding of
        # relevant line lists.
        ###====================================================================
        df = indict[larsno]
        # Make sure all lines are present
        df = normalize_lines(df, lineslist=ultilist)  # Ulti right choice?
        # Now, check that at least one of the lines requested for the covering
        # fraction calculations exists. Why not two? Because the later code
        # will take care of that, but it cannot handle an empty list and it
        # doesn't make sense anyway. So here> Empty lists out!
        # TODO: Do the check, don't just announce it.
        checklist = set([line.replace(' ', '_') for line in lines])
        availlist = list(df.columns)
        for line in [a.replace(' ', '_') for a in lines]:  # df.columns:
            if not df[line].notnull().any():
                availlist.pop(availlist.index(line))
                # df = df.drop(line, axis=1)
                # df = df.drop(line + '_Stddev', axis=1)
                # try:
                #     df = df.drop(line + '_Mask', axis=1)
                # except:
                #     contine
        datacols = set(availlist)
        if checklist.isdisjoint(datacols):
            ax = fig.add_subplot(5, 3, larsno)
            ax.axhline(y=0., color='k', ls='-')
            ax.axhline(y=1., color='k', ls=':')
            ax.axvline(x=0., color='k', ls=':')
            ax.axis((-900, 900, -.2, 1.5))
            ax.annotate(df.columns.name, (0.05, 0.8), xycoords='axes fraction')
            continue
        # Make a grid for the input lines
        vf = velocity_regrid(df, lineslist=lines)
        # ...for the lines used for the average profile:
        lf = velocity_regrid(df, lineslist=compare_lines)
        # ...and for the union of the two lists
        ff = velocity_regrid(df, lineslist=ultilist)

        print "\n Did velocity regridding of {} \n \n".format(df.columns.name)

        ###============================================
        # Make Select the subset of the lines that are
        # actually in the DataFrame.
        ###============================================
        toavg = []
        toerr = []
        avmask = []
        # Which lines to include?
        includelines = set(lines).union(set(compare_lines))
        for line in includelines:
            if line.replace(' ', '_')+'_Flam' in lf.columns:
                toavg.append(line.replace(' ', '_')+'_Flam')
                toerr.append(line.replace(' ', '_')+'_Stddev')
                avmask.append(line.replace(' ', '_')+'_Mask')
        ###============================================
        # Mean profile:
        ###============================================
        avarr = ff[toavg].values
        erarr = ff[toerr].values
        maarr = ff[avmask].values
        avmar = sp.ma.masked_array(avarr, sp.invert(maarr))
        ermar = sp.ma.masked_array(erarr, sp.invert(maarr))
        avgpro = avmar.mean(axis=1)
        errpro = sp.sqrt((ermar**2.).sum(axis=1))
        # # XXX DEBUG:
        # fig = plt.figure(1)
        # axes = plt.gca()
        # axes.plot(avgpro, label='new')
        # plt.legend().draggable()
        # import pdb; pdb.set_trace()  # XXX BREAKPOINT
        # # END DEBUG
        ax = fig.add_subplot(5, 3, larsno)
        ax.axhline(y=0., color='k', ls='-')
        ax.axhline(y=1., color='k', ls=':')
        ax.axvline(x=0., color='k', ls=':')
        ax.annotate(df.columns.name, (0.05, 0.8), xycoords='axes fraction')
        if not sp.mod(larsno, 3).astype(bool):
            ax.tick_params(labelleft='off', labelright='on')
        if not sp.mod(larsno + 1, 3).astype(bool):
            ax.tick_params(labelleft='off', labelright='off')
        if larsno < 12:
            ax.tick_params(labelbottom='off')
        ax.axis((-900, 900, -.2, 1.5))
        ###=========================
        # plot average line profile
        ###=========================
        # print lf.Velocity.values
        # print avgpro
        # print len(lf.Velocity.values), len(avgpro)
        #
        print 'AVGPRO IS:', type(avgpro)
        print lf.Velocity.values.shape, avgpro.shape
        avgplot = ax.plot(
            ff.Velocity.values, avgpro,#.values,
            drawstyle='steps',
            color='black',
            label='Average Si II line profile')
        ###=========================
        # Calculate cover fractions:
        ###=========================
        if inveldf is None:  # FIXME: Gets a dict, expects a DF. Make expct DF
            print 'Now calculating covering fractions for {}'\
                    .format(df.columns.name)
            veldf = cover_fract(
                ff[vf.columns], species=species,
                illustr=False, minimizer=minimizer)
            if minimizer=='brute':
                brutedict = veldf[1]
                veldf = veldf[0]
            print 'Shape of covfrac input: ', ff[vf.columns].shape
            print type(veldf)# , veldf
            print 'Shape of coverfrac output: ', veldf.shape
        idx = veldf[veldf[species+'_NumLines'].notnull()].index
        covfrac = 1. - veldf[species+'_CovFrac_map'].ix[idx].values
        coverr = veldf[species+'_CovFrac_err'].ix[idx].values
        numlines = veldf[species+'_NumLines'].ix[idx].values
        chisq = veldf[species+'_CovFrac_Xsq'].ix[idx].values
        uperr = covfrac + coverr
        loerr = covfrac - coverr
        # `if both upper and lower errors exist, utilize these:
        if species+'_CovFrac_lo_err' in veldf.columns:
            loerr = covfrac - veldf[species+'_CovFrac_lo_err'].ix[idx].values
            print 'LOWER ERRORS LARS', larsno, loerr.min(), loerr.max()

        plotvels = veldf.Velocity.ix[idx].values
        print 'Vel, ChiSq: ', plotvels.shape, chisq.shape
        uperr = sp.hstack(
            (uperr.reshape(-1, 1),
             sp.ones_like(uperr.reshape(-1, 1)))).min(axis=1)
        loerr = sp.hstack(
            (loerr.reshape(-1, 1),
             sp.zeros_like(loerr.reshape(-1, 1)))).max(axis=1)
        ###============================================
        # Plot dem coverfracs and uncertainties
        ###============================================
        cfplot = ax.plot(
            plotvels, covfrac,
            'o', color='green',
            mec='green',
            mew=1.,
            markersize=4,
            #alpha=.6,
            label='Si II covering fraction')
        nlplot = ax.scatter(
            plotvels,
            sp.ones_like(plotvels) * -.1,
            c=[['yellow', 'orange', 'red', 'purple', 'blue'][int(i)]
               for i in numlines],
            s=30,
            marker='.', linewidths=0., #cmap='heat',
            label=r'\# of lines included in fit')
        #chisqplot = ax.plot(
        #    plotvels, chisq * .25, '-',
        #    color='green', alpha=.6)
        #print plotvels.shape, uperr.shape, loerr.shape
        uperrplt = plt.bar(
            plotvels, uperr-loerr,
            width=15.6605, bottom=loerr,
            color='gray', linewidth=0.,
            align='center', alpha=.8,
            label=r'$1 \sigma$ range',)
        print larsno, 'Is the number LARSNO'
        if larsno == max(indict.keys()):  # 14:
            leg = ax.legend(
                frameon=False,
                fontsize='medium',
                loc='center',
                bbox_to_anchor=(1.5, 0.5),
                labelspacing=.2)
            leg.draggable()
        # Now, combine the two dataframes to make sure to have ALL the data.
        for col in lf.columns:
            veldf.loc[:, col] = ff[col]  # FIXME: I think this is where it fails.
        print veldf.shape, avgpro.shape  # Take from ff for consistency, but
        veldf.loc[:, 'MeanProfile'] = avgpro    # only values present in lf
        veldf.loc[:, 'MeanErrors'] = errpro
        veldfdict[larsno] = veldf
    fig.subplots_adjust(hspace=.1, wspace=.04)
    fig.savefig('Figures/covfracfig.pdf')
    fig.show()
    if return_brutedict:
        return brutedict, veldfdict
    else:
        return veldfdict


def handle_brute_dict(indict):
    """ Maybe a docstring some day.
    """
    data = indict.copy()
        # Clean the list of keys of unwanted elements:
    if 'covfrac' in data.keys():
        covfrac = data.pop('covfrac')
    if 'logN' in data.keys():
        logn = data.pop('logN')
    x2_min = []
    keys = data.keys()
    keys.sort()
    cofralst = []
    coldelst = []
    chisqlst = []

    # print keys
    # import ipdb; ipdb.set_trace()  # XXX BREAKPOINT

    for key in keys:
        frame = data[key]
        min_val = frame.min()
        minix = sp.where(frame < min_val + 1.)
        pointix = sp.where(frame == min_val)
        x2m = [key, min_val]
        chisqlst.append(x2m)
        # [Wavelength, best fit, min conf, max conf]
        cf = [key,
              covfrac[pointix[1]][0],
              covfrac[minix[1]].max(),
              covfrac[minix[1]].min()]
        logN = [key,
                logn[pointix[0]][0],
                logn[minix[0]].max(),
                logn[minix[0]].min()]
        # Now, if there were only 2 datapoints, make NAN's
        if (frame.min() == 1000.) & (frame.ptp() == 0.):
            # print "There was an 'empty' frame!"
            cf = [key] + [sp.nan] * 3
            logN = [key] + [sp.nan] * 3
        cofralst.append(cf)
        coldelst.append(logN)
    cofralst = sp.array(cofralst)
    coldelst = sp.array(coldelst)
    chisqlst = sp.array(chisqlst)
    # To fix errors to be relative to the measured value,
    # rather than absolute numbers:
    cofralst[:, 3] = cofralst[:, 1] - cofralst[:, 3]
    cofralst[:, 2] -= cofralst[:, 1]
    coldelst[:, 3] = coldelst[:, 1] - coldelst[:, 3]
    coldelst[:, 2] -= coldelst[:, 1]
    return cofralst, coldelst, chisqlst


def view_brute_dicts(indict):
    """ Maybe a docstring one day
    """
    import adhoc_functions as ad
    reload(ad)
    data = indict.copy()
    if 'covfrac' in data.keys():
        covfrac = data.pop('covfrac')
    if 'logN' in data.keys():
        logn = data.pop('logN')
    thething = ad.CubeView(data, covfrac, logn)
    thething.configure_traits()
    return


def draw_coverfrac_axes(veldf, ax1, ax2, ax3, ms=3,
                        species='Si_II', legend=True, thick=.8, thin=.2,
                        dotcolor='red', barcolor='pink', cf_lims=None,
                        **kwargs):
    # veldf = veldf.dropna(subset=[species+'_CovFrac_err'])
    avgprofile = veldf.MeanProfile.values
    avgproerrs = veldf.MeanErrors.values
    velocities = veldf.Velocity.values
    shortvels = velocities
    if cf_lims is not None:
        veldf = veldf[veldf.Velocity.between(cf_lims[0], cf_lims[1])]
        shortvels = veldf.Velocity
    coverfracs = veldf[species+'_CovFrac_map'].values
    #coverfracs = veldf[veldf.Velocity.between(cf_lims[0], cf_lims[1])]\
    #    [species+'_CovFrac_map'].values
    coverfrerr = coverfracs + sp.absolute(veldf[species+'_CovFrac_err'].values)
    columndens = veldf[species+'_ColDens_map'].values
    columnderr = columndens + sp.absolute(veldf[species+'_ColDens_err'].values)
    chisquared = veldf[species+'_CovFrac_Xsq'].values

    #avgprofile = sp.ma.masked_array(avgprofile, veldf.)

    # Make sure annotations get correct.
    if species == 'Si_II':
        phase = 'LIS'
    elif species == 'Si_IV':
        phase = 'HIS'

    if species+'_CovFrac_lo_err' in veldf.columns:
        covfrloerr = coverfracs - veldf[species+'_CovFrac_lo_err'].values
        coldeloerr = columndens - veldf[species+'_ColDens_lo_err'].values
    else:
        covfrloerr = coverfracs - sp.absolute(coverfrerr)
        coldeloerr = columndens - sp.absolute(columnderr)

    chisqplot = ax1.plot(
        shortvels, chisquared, 'k-', label='$\chi^2 / # datapoints', lw=thick)
    ax1.axvline(x=0., color='k', linewidth=.5, linestyle='--')
    ax2.axvline(x=0, color='k', ls='--', lw=.5)
    ax3.axhline(y=0, color='black', lw=.5, zorder=0.)
    ax3.axhline(y=1, color='black', lw=.5, linestyle='--', zorder=0.)
    ax3.axvline(x=0, color='k', ls='--', lw=.5, zorder=0.)
    #ax1.axhline(y=0, color='k', lw=1)

    #print velocities.shape, columnderr.shape, coldeloerr.shape
    if len(velocities) > 2:
        ax1.set_ylim(-.5, chisquared.ptp()*.8)
        shortvels = sp.array(shortvels)
        # print('cospaper.py sez:')
        # print shortvels, type(shortvels)
        # print columnderr, type(columnderr)
        # import pdb; pdb.set_trace()  # XXX BREAKPOINT
        cderr = mp.fill_between_steps(
            shortvels, columnderr, y2=coldeloerr,
            ax=ax2,
            color='lightgray', edgecolor='gray', lw=thin)
        cdplot = ax2.plot(
            shortvels, columndens, label='log$N_{S_{II}}$',
            color='k', drawstyle='steps-mid', lw=thick)
        ax2.set_ylim(columndens.min()-1., columndens.max()+2.0)

        cofrplot = ax3.plot(
            shortvels, 1.-coverfracs, 'o',
            mec=dotcolor, mfc=dotcolor, alpha=.99, mew=.8,
            label=r'$f_c$', zorder=2, ms=ms)
        cferplot = mp.fill_between_steps(
            shortvels, 1-coverfrerr, 1-covfrloerr, ax=ax3,
            facecolor=barcolor, edgecolor='gray', alpha=.7, zorder=2,
            lw=thin)
            #label=r'$f_c: \Delta X^2 < 1$')#edgecolor='red', alpha=.5)
        avprplot = ax3.plot(
            velocities, avgprofile, 'k-', drawstyle='steps-mid', zorder=1,
            label=r'$I/I_{0, \rm{'+phase+'}}$', lw=thick)
        averplot = mp.fill_between_steps(
            velocities, avgprofile+avgproerrs, avgprofile-avgproerrs,
            ax=ax3, #label='Average profile conf. levels',
            color='lightgray', edgecolor='gray', zorder=1, lw=thin)
        larsname = veldf.columns.name
        if legend:
            leg = ax3.legend(
                fancybox=True, shadow=False, loc='lower left',
                fontsize=8)
            leg.set_title(larsname, prop={'size':7})
            leg.draggable()
    # TODO: Proxies to get fill_between in legend.

    return ax1, ax2, ax3


def single_coverfrac_plot(
    veldf, species='Si_II', save='False', sav_format='png',
    dotcolor='red', barcolor='pink'):
    ''' Put docstring here'''
    larsname = veldf.columns.name
    fig = plt.figure('Single N-Fc-X^2 plot, '+larsname,
                     figsize=(3.4, 3.4), dpi=196)
    ax1 = plt.subplot(411)
    ax2 = plt.subplot(412, sharex=ax1)
    ax3 = plt.subplot(212, sharex=ax1)
    ax4 = plt.twinx(ax3)

    draw_coverfrac_axes(veldf, ax1, ax2, ax3, species=species,
                        interpolate=False, dotcolor=dotcolor, barcolor=barcolor)
    fig.subplots_adjust(hspace=0, right=.87, bottom=.11)
    ax1.tick_params(labelbottom='off', labelsize=8)
    ax1.locator_params(nbins=4, axis='y')
    ax2.tick_params(labelbottom='off', labelsize=8)

    ax1.set_ylim(0, 80)
    ax2.set_ylim(9, 16)
    ax3.set_xlim(-900, 600)  # Yeah, hardcoding. A sane default, though.
    ax3.tick_params(labelsize=8)
    ax4.tick_params(labelsize=8)
    ax3.set_ylim(-.2, 1.6)
    ax4.set_ylim(1.2, -.6)
    ax3.yaxis.set_ticks([0., .2, .4, .6, .8, 1.])
    ax4.yaxis.set_ticks([0., .2, .4, .6, .8, 1.])
    ax2.yaxis.set_ticks([10., 11., 12., 13., 14., 15])
    ax4.yaxis.label.set_color(dotcolor)
    ax4.tick_params(axis='y', colors=dotcolor)
    ylo, yup = ax3.get_ylim()
    # ax4.set_ylim(1-ylo, 1-yup)
    # ax3.set_xlim(veldf.Velocity.min()-20, veldf.Velocity.max()+20)
    ax4.set_ylabel('Covering fraction', fontsize=9)
    ax3.set_xlabel('v [km/s]', fontsize=9)
    ax3.set_ylabel('Relative flux', fontsize=9)
    ax2.set_ylabel(
        '$\log N_{\mathrm{'+'{}'.format(species.replace('_', ''))+'}}$',
        fontsize=9)
    ax1.set_ylabel(r'Pseudo-red. $\chi ^2$', fontsize=9)

    if save:
        filename = 'Figures/CoverFracs/{}_{}.{}'.format(
            larsname.replace(' ', '_'), species, sav_format)
        plt.savefig(filename)

    plt.show()
    return


def set_mask(indf, lines, skip=True, mode='overwrite'):
    '''
    Edits the masking array of a line.
    Input must be a dataframe, a list of
    lines of the
    Possible modes:
    'append'
    :    Adds selection to mask (i.e. sets it to False)
    'overwrite'
    :    Flushes old array, starts over from scratch.
    '''
    larsno = indf.columns.name.split()[-1]
    bbox_args = dict(boxstyle='round', fc='w')
    if not type(lines) == list:
        lines = [lines]  # Assume string.

    def _onselect(vmin, vmax):
        ###============================================
        # Set everytging between vmin and vmax to false.
        ###============================================
        idx = indf[(indf['lambda'] > vmin) & (indf['lambda'] < vmax)].index
        indf[i+'_Mask'].ix[idx] = False
        print idx.values.min(), idx.values.max()
        print sp.where(indf[i+'_Mask'] == True)[0].shape
        print len(sp.where(indf[i+'_Mask'] == True)[0])
        ax.axvspan(vmin, vmax, color='darkorange', alpha=.2, zorder=9)
        ax.set_xbound(plmin, plmax)
        ax.set_ybound(0, 2.5)
        fig.canvas.draw()
        return

    def _onbutton(event):
        fig.savefig(filepath)
        print 'Saved'
        return

    for i in lines:
        filepath = 'Figures/masks/{0}{1}.pdf'.format(
            indf.columns.name.replace(' ', '_'), i.replace('_', ' '))
        if len(i.split(' ')) > 1:  # Check for correct string format
            i = i.replace(' ', '_')
        if not i in indf.columns or not indf[i].notnull().any():
            print 'line array doesn not exist or is empty; continuing.'
            continue
        if i+'_Mask' in indf.columns and skip:
            continue
        if not i+'_Mask' in indf.columns or mode=='overwrite':
            indf[i+'_Mask'] = False
            indf[i+'_Mask'].ix[indf[i].dropna().index] = True
        barr = sp.ones(len(indf[i].dropna())).astype(bool)
        fig = plt.figure()
        ax = plt.subplot(111)
        axsave = plt.axes([0.8, 0.8, 0.09, 0.09])
        savb = Button(axsave, 'Save')
        savb.on_clicked(_onbutton)

        if (i+'_Mask' in indf.columns) and (not skip):
            print 'Existing mask for {}'.format(i)
            themask = indf[i+'_Mask'].values
            thelmbd = indf['lambda'].values
            LMBD = sp.ma.masked_array(thelmbd, themask)
            ax.fill_between(LMBD, -1, 3, color='darkorange', alpha=.2)

        plmin = indf[i].dropna().index.values.min()
        plmax = indf[i].dropna().index.values.max()

        for theline in MWlines.keys():
            #if MWlines[theline] < plmin or MWlines[theline] > plmax:
            #    continue
            if theline.replace(' ', '_') == i:
                lc = 'red'
            else:
                lc = 'green'
            ax.axvline(
                x=MWlines[theline], ymin=.5, ymax=.99,
                lw=1.5, color='b', ls=':')
            ax.annotate(
                theline+' MW', (MWlines[theline]-.0, 2.2),
                backgroundcolor='white',
                rotation='vertical', horizontalalignment='center', color='b')
            ax.axvline(
                x=MWlines[theline] * (indf.z_value + 1.),
                ymin=.5, ymax=.99, color=lc, ls='--', lw=1.5)
            ax.annotate(
                theline+' L'+larsno,
                (MWlines[theline] * (indf.z_value + 1.) - .2, 2.2),
                backgroundcolor='white',
                rotation='vertical', horizontalalignment='center', color=lc)

        ax.plot(indf['lambda'], indf[i], color='k', drawstyle='steps-mid')
        print indf['lambda'].shape, indf[i+'_Stddev'].shape
        mp.fill_between_steps(indf['lambda'].values, indf[i+'_Stddev'].values,
                              ax=ax, color='gray', alpha=.6, interpolate=False)
        ax.axhline(y=0, linewidth=1, color='k')
        ax.set_xbound(plmin, plmax)
        ax.set_ybound(0, 2.5)
        span = SpanSelector(ax, _onselect, 'horizontal',
                            useblit=True, minspan=.5)
        #fig.savefig(filepath)
        plt.show()
        print filepath
    return indf


def pickle_df(indf, filename='dfdict.pickle'):
    import pickle
    outfile = open(filename, 'w')
    pickle.dump(indf, outfile)
    outfile.close()
    return


def velwidth_table(vfdict):  # OLD: Use notebook instead.
    velwidict = {}
    for i in vfdict.keys():
        larsname = vfdict[i].columns.name
        wavedict = dv.dynvel(vfdict[i])
        velwidict[larsname] = wavedict

    velwifr = pd.DataFrame.from_dict(velwidict).T
    # print velwifr
    outfr = velwifr.copy()
    outfr['I/I0_min'] = outfr['I/I0_min'].map(lambda x: sp.around(x, 2))
    outfr['I/I0_0'] = outfr['I/I0_0'].map(lambda x: sp.around(x, 2))
    outfr.v_min = outfr.v_min.map(lambda x: sp.around(x))
    outfr.v_95 = outfr.v_95.map(lambda x: sp.around(x))
    outfr.v_int = outfr.v_int.map(lambda x: sp.around(x))
    outfr.delta_v = outfr.delta_v.map(lambda x: sp.around(x))
    return outfr


def fit_flourescent(indf, smooth=None, interactions=1):
    '''
    Function to make simple fits of flourescent lines of normalized specs.
    Lines to look at:
        C II 1334, Si II 1190, Si II 1193, O I 1302, Si II 1306
    Cycles through lines and fits a model to each of them.
    '''
    ###=======================================================================
    #     Functions defining the behaiiour of buttons and span selector.
    ###=======================================================================

    interactions = float(interactions)

    def _onselect(vmin, vmax):
        ax.axvspan(vmin, vmax, color='lightgreen',  zorder=0)
        thidx = sp.where((wave > vmin) & (wave < vmax))
        idx[thidx] = True
        fig.canvas.draw()
        return

    def _on_ok_button(event):
        print fitresults
        output_dict[line] = fitresults
        # plt.close(fig)
        fitresults['LineFik'] = flourlines.Aki.ix[line]
        if not line in lmfit_fitters.keys():
            lmfit_fitters[line] = None
        try:
            fitresults['RedXSq'] = lmfit_fitters[line].redchi
        except:
            fitresults['RedXSq'] = None
        # Calculate velocity relative to systemic:
        #vels = (wavlens - obswl) / obswl * con.c / 1000.
        try:
            pecwl = fitresults['centroid']
            #pecer = fitresults['centroer']
            #pecwl = ufl(pecwl, pecer)
            syswlidx = flourlines[sp.absolute(flourlines['ObsWl'] - pecwl) ==
                sp.absolute(flourlines['ObsWl'] - pecwl).min()].index
            fitresults['FlourLine'] = syswlidx.values[0]
            syswl = flourlines['ObsWl'].ix[syswlidx].values
            offsetvel = (pecwl - syswl) / syswl * con.c / 1000.
            fitresults['VelOffset'] = offsetvel[0]
            # Calculate EW
            EW = fitresults['amp'] * fitresults['sigma'] * sp.sqrt(2 * sp.pi)
            fitresults['EquiWidth'] = EW
            fitresults['FlourFik'] = flourlines.Aki.ix[syswlidx][0]
            print syswlidx, syswl, offsetvel
            fitresults['RelFik'] = (fitresults['LineFik'] / \
                    fitresults['FlourFik'])
            fitresults['ResEW'] = \
                    fitresults['EquiWidth'] * fitresults['RelFik']

            ###================================================================
            #  Experimental: take multiple sscatteringsinto account
            ###================================================================
            #      IS THIS CORRECT?!?
            # ===============================
            # (I think so but needs checking)

            fitresults['ResEW1'] = \
                    fitresults['EquiWidth'] * fitresults['RelFik'] * \
                    (fitresults['LineFik'] /
                        (fitresults['LineFik'] +
                         fitresults['FlourFik'])) ** 1  # (interactions - 1.)

            fitresults['ResEW2'] = \
                    fitresults['EquiWidth'] * fitresults['RelFik'] * \
                    (fitresults['LineFik'] /
                        (fitresults['LineFik'] + fitresults['FlourFik'])) ** 2  # (interactions - 1.)

            fitresults['ResAmp'] = \
                    fitresults['ResEW'] / (fitresults['sigma'] *
                                           sp.sqrt(2. * sp.pi))
            fitresults['ResAmp1'] = \
                    fitresults['ResEW1'] / (fitresults['sigma'] *
                                           sp.sqrt(2. * sp.pi))
            fitresults['ResAmp2'] = \
                    fitresults['ResEW2'] / (fitresults['sigma'] *
                                           sp.sqrt(2. * sp.pi))
            syn_gauss = pro.gauss(wave, center,
                                  fitresults['sigma'].nominal_value,
                                  fitresults['ResAmp'].nominal_value)
            syn_gauss1 = pro.gauss(wave, center,
                                  fitresults['sigma'].nominal_value,
                                  fitresults['ResAmp1'].nominal_value)
            syn_gauss2 = pro.gauss(wave, center,
                                  fitresults['sigma'].nominal_value,
                                  fitresults['ResAmp2'].nominal_value)
            ax.plot(wave, syn_gauss+1., color='green')
            ax.plot(wave, syn_gauss1+1., color='blue')
            ax.plot(wave, syn_gauss2+1., color='purple')
            ax.plot(wave, pltflux-syn_gauss,  color='green',
                    drawstyle='steps-mid', lw=1.2, zorder=1)
            ax.plot(wave, pltflux-syn_gauss1, color='blue',
                    drawstyle='steps-mid', lw=1.2, zorder=1)
            ax.plot(wave, pltflux-syn_gauss2, color='purple',
                    drawstyle='steps-mid', lw=1.2, zorder=1)
            plt.draw()
        except:
            print "Could not compute flourescent velocity offset of line {}"\
                    .format(line)
            raise
        return



    def _on_clr_button(event):
        if not clear_lock[0]:
            ax.patches = []
            ax.lines.pop()
            ax.lines.pop()
            idx[:] = False
            fig.canvas.draw()
            clear_lock[0] = True
        return

    def _flour_model(pars, wave, data, stddev):
        amp = pars['amp'].value
        mu = pars['centroid'].value
        sigma = pars['sigma'].value
        model = pro.gauss(wave, mu, sigma, amp) + 1.
        resid = sp.sqrt((model - data)**2/stddev**2)
        return resid

    def _on_fit_button(event):
        # global result  # Just for debugging!
        fitdata = flux[idx]
        fiterrs = errs[idx]

        mu_guess = sp.median(wave[idx])
        sigma_guess = 1.
        amp_guess = flux[idx].max() - 1.
        fig.canvas.draw()
        pars = lm.Parameters()
        pars.clear()
        pars.add('amp', amp_guess)
        pars.add('sigma', sigma_guess)
        pars.add('centroid', mu_guess)
        result = lm.minimize(
            _flour_model, pars, args=(wave[idx], flux[idx], errs[idx]))
        lmfit_fitters[line] = result
        fitresults['sigma'] = ufl(pars['sigma'].value, pars['sigma'].stderr)
        fitresults['amp'] = ufl(pars['amp'].value, pars['amp'].stderr)
        fitresults['centroid'] = ufl(
            pars['centroid'].value, pars['centroid'].stderr)
        fitted_model = pro.gauss(
            wave, pars['centroid'].value,
            pars['sigma'].value, pars['amp'].value) + 1.
        ax.plot(wave, fitted_model, color='orange', lw=2.5, zorder=6)
        fig.canvas.draw()
        clear_lock[0] = False
        return

    def load_flourlines():
        flourlines = pd.read_table(
            'lislines_NIST.ascii',
            delimiter='|',
            skiprows=[1, 6, 10, 16, 19],
            index_col=0)
        flourlines.columns = flourlines.columns.map(str.strip)
        flourlines.index.name = flourlines.index.name.strip()
        flourlines.index = flourlines.index.map(str.strip)
        flourlines['ObsWl'] = flourlines['RestWl (AA)'] * (1. + indf.z_value)
        #print flourlines
        return flourlines

    def _dummy_func():
        """ Only to get correct folding in Vim,
        ignore and remove when no longer necessary.
        """
        return

    wave = indf['lambda'].values
    fitresults = {}
    output_dict = {}
    clear_lock = [True]
    lmfit_fitters = {}
    flourlines = load_flourlines()
    for i in range(len(SiIIlines)):
        line = SiIIlines[i]
        center = wlsdict[line] * (1 + indf.z_value)
        flux = indf['_'.join(line.split())].values
        fitresults = {}
        # If the smoothing flag is set, smooth the visualized data, but still
        # use the non/smoothed data for fitting.
        if smooth is not None:
            try:
                kernel = sp.ones(smooth) / smooth
            except:
                kernel = sp.ones(5) / 5.
            pltflux = sp.convolve(flux, kernel, mode='same')
        else:
            pltflux = flux
        errs = indf['_'.join(line.split() + ['Stddev'])].values
        idx = sp.zeros_like(flux).astype(bool)

        ###==========
        # Define plot
        fig, ax = plt.subplots(1, 1, figsize=(7, 6), dpi=160)
        # Define axes for buttons:
        axfit = plt.axes([0.9, 0.84, 0.09, 0.06])
        axclr = plt.axes([0.9, 0.77, 0.09, 0.06])
        ax_ok = plt.axes([0.9, 0.70, 0.09, 0.06])
        # define buttons:
        fitbu = Button(axfit, 'Fit')
        clrbu = Button(axclr, 'Clear')
        ok_bu = Button(ax_ok, 'OK')
        # What to do when buttons are clicked:
        fitbu.on_clicked(_on_fit_button)
        clrbu.on_clicked(_on_clr_button)
        ok_bu.on_clicked(_on_ok_button)
        # Plot data & erro)s
        ax.plot(wave, pltflux, color='black', label=line, lw=1.6,
                drawstyle='steps-mid', zorder=4)
        mp.fill_between_steps(wave, pltflux-errs, pltflux+errs,
                              ax=ax, color='lightgrey', zorder=1)
        # Set markers for expected lines at systemic velocity
        for flourline in flourlines.index:
            # print flourline, line, flourline == line, len(flourline), len(line)
            if flourline == line:
                continue
            ax.axvline(
                x=flourlines.ix[flourline]['ObsWl'],
                color='gray',
                ls='--',
                ymin=.6)
            ax.annotate(flourline,
                        xy=(flourlines.ix[flourline]['ObsWl']-0.5, 0.9),
                        xycoords=('data', 'axes fraction'),
                        color='gray',
                        rotation=90,
                        ha='center',
                        )
        ax.axhline(y=0, lw=1., color='k')
        ax.axhline(y=1, lw=1., color='k', ls='--')
        ax.axvline(x=0, lw=1., color='k', ls='--')
        #ax.axvline(x=center, color='gray', ls='--', lw=2.5)
        ax.axvline(x=center, color='maroon', ls='--', lw=1.)
        ax.annotate(line, xy=(center-.5, .9), xycoords=('data', 'axes fraction'),
                    color='maroon', rotation=90)
        ax.legend(loc='lower left')
        ax.axis((center-10, center+10, -.2, 2.))
        span = SpanSelector(ax, _onselect, 'horizontal',
                            useblit=True, minspan=.05)
        fig.subplots_adjust(left=.1, right=.89)
        plt.show()

    return pd.DataFrame.from_dict(output_dict), flourlines  #, lmfit_fitters
