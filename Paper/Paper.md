---
title: "Neutral ISM, Lyman-Alpha and Lyman-continuum in nearby starburst Haro
        11" 
titlefootnote:
  - mark: 1
    text: 'Based on observations Cosmic Origins Spectrograph on the Hubble Space Telescope, program GO 13017, PI Timothy Heckman'

shorttitle: "Neutral ISM, Ly-Alpha and LyC in Haro 11"
author:
  - name: T. Emil Rivera-Thorsen
    altaffil: '2, 3'
  - name: Göran Östlin
    altaffil: '2, 3'
  - name: Matthew Hayes
    altaffil: '2, 3'

altaffiltext:
  - num:  '2'
    text: Department of Astronomy, Oskar Klein Centre, Stockholm University,
           AlbaNova University Centre, SE-106 91 Stockholm, Sweden
  - num:  '3'
    text: Oscar Klein Centre for Cosmoparticle Physics, Department of
          Astronomy, Stockholm University, Stockholm, Sweden

shortauthor: 'T. E. Rivera-Thorsen et al.'
biblio-style: aasjournal
header-includes: \usepackage[caption=false]{subfig}
---


## Abstract:

Star forming galaxies are believed to be the main source of ionizing radiation
responsible for reionizing the early Universe. Direct observations have however
been few and with escape fractions far below what would be necessary to account
for the ionization history of the Universe. Especially in the local Universe,
only a few handfuls of emitters have been observed with typical escape
fractions of a few percent. It seems the mechanisms regulating this escape need
to be ones that are strongly evolving with redshift, if galaxies are to be held
responsible for the reionization. Gas mass, star formation feedback and thus
star formation activity are among the main suspects, known to both regulate
neutral gas coverage and evolve with cosmic time. In this paper, we present a
reanalysis of far-UV HST-COS spectra of the first detected Lyman-continuum
leaker in the local Universe, Haro 11, in which we examine the connections
between Lyman-continuum leakage and Lyman-$\alpha$ line shape and
feedback-influenced neutral ISM properties like kinematics, geometry and column
density. In particular, we discuss the two extremes of an optically thin,
density bounded ISM and a riddled, optically thick, ionization bounded ISM and
how Haro 11 fit into the theoretical predictions made from these models. We
find that the most likely ISM model for Haro 11 knot C is one of clumpy neutral
medium with a highly ionized interclump medium with a combined covering
fraction of near-Unity and a with a residual neutral gas column density of $XX
N_{HI} < XX$, high enough to be optically thick to Lyman-$\alpha$, but low
enough to be at least partly transparent to Lyman-continuum.


Introduction and Observations
-------------------------------------------------------------------------------

### Background


### About Haro 11

Some of the same stuff as in Rivera-Thorsen et al. (in prep.) could be
mentioned, with a different emphasis. Star formation rate should be mentioned
from @Adamo2010. Population from same and from @Micheva2010. Knots named by
@Vader1993, and also mentioned by @Kunth2003. It's an LBA [@Hoopes2007], and an
LAE analog [@Hayes2007; @Leitet2011]. Knot C a strong Ly$\alpha$ emitter
[@Hayes2007], and also a LyC leaker [@Bergvall2006; @Leitet2011]. Kinematics,
it's a merger [@Ostlin1999; @Ostlin2001; @Ostlin2015]. More weight on the
Lyman-$\alpha$ escape thing, and LyC. Mention Ly$\alpha$ escape and the history
here, too? References: 

This observation is first (?) published in @Alexandroff2015, and also in
@Verhamme2015 in which the Ly$\alpha$ profile is dicussed as a .

This observation is part of the sample in @Alexandroff2015, in which details
about observation and data reduction are described in depth; we pint the reader
there for further information about these. The Ly$\alpha$ profile in this
spectrum is discussed in @Verhamme2015. @Heckman2011, @Heckman2015.
@Jaskot2013, @JaskotOey. @Bouwens2012.


#### Figure: HaroCOSApert                                          {#fig:apert}

![](../Figs/Haroslit.pdf){width=3.5in}

Caption: Approximate position of the COS aperture, shown on HST imaging data of
@Hayes2009; @Ostlin2009, encoding UV continuum in green, H$\alpha$ in red, and
continuum subtracted Lyman $\alpha$ in blue. 


In this paper, we analyze HST-COS spectra of knot C in Haro 11. The
observations were done as part of HST program GO 13017, PI Timothy Heckman. The
spectrum has first been published by @Alexandroff2015, a sample paper in which
an attempt is made of 



Analysis 
-------------------------------------------------------------------------------

### Individual lines

Figure [#fig:SingleLines] shows the individual profiles of the transitions
included in our analysis; the upper panel shows transitions of \ion{Si}{2},
lower panel of \ion{Si}{4}. It is plainly visible that the ionization fraction
is high, with the \ion{Si}{4} curves being considerably deeper than the
low-ionized lines.  Looking at the upper panel, \ion{Si}{2} $\lambda \lambda
1304, 1526$ are somewhat shallower than \ion{Si}{2} $\lambda 1260$. The former
two lines have comparable resonator strengths, both about a factor of 10 lower
than that of $\lambda 1260$. It is thus clear that we do not find ourselves in
the optically thin regime, in which the latter line should be correspondingly
around 10 times stronger; on the other hand, it is possible that the two weak
lines are not completely saturated. In the lower panel, the two lines have
resonator strengths within a factor of 2 of each other. They are thus at first
glance consistent with a medium that is not completely opaque, but not with an
optically thin one. The stronger absorption in \ion{Si}{4} reveals a high level
of ionization of the medium covering the central cluster. 


#### Figure: HisLisProf.                                     {#fig:SingleLines}

![](../Figs/HISLISProfiles.pdf){width=3.5in}

Caption: The \ion{Si}{2} (**upper**) and \ion{Si}{4} (**lower**) profiles
included in this study. 

\input{table1}


### $N_{\rm Si}$ and $f_C$                                           {#sec:aod}

Following the method described in @RiveraThorsen2015; we have performed fits
for column density and covering factor in each velocity bin, for both the high-
and low-ionization state. The results are shown in fig. [#fig:WithColDens]. The
upper panels show the pseudo-reduced $\chi^2$ as defined in @RiveraThorsen2015
( $=\chi^2 / (\mathrm{DOF} + 1)$ ) for each bin, middle panels show the
inferred column density in each bin, with surrounding shaded columns showing
the confidence intervals. In the lower panels, the mean LIS line profile is
shown in black with gray shaded uncertainty intervals. On these are overlaid
the best-fit values of $f_C$ as colored dots, with surrounding shaded bars
showing the confidence intervals. We again caution that $f_C$ is the covering
fraction of HI atoms *within the given velocity bin*, and hence only provides a
lower limit for the total, geometric neutral gas covering fraction, since gas
at different velocities generally does not occupy the same projected area.


#### Figure: WithColumnDensity                         {#fig:WithColDens .wide}

![LIS](../Figs/Fc_haro11c_LIS.pdf){#coldenLIS width=40%}
![\ion{Si}{4}](../Figs/Fc_haro11c-SiIV.pdf){#coldenHIS width=40%}

Caption: **Upper panels**: Pseudo-reduced $\chi^2$ as described in
@RiveraThorsen2015. **Middle panels**: Best-fit ion column density with
confidence intervals in shaded gray. **Lower panels**: Mean LIS/\ion{Si}{4}
profile shown as black steps, with inferred $f_C$ shown with yellow dots.
Lighter shaded columns show confidence intervals for both. 



Discussion and conclusions 
-------------------------------------------------------------------------------

### Lyman-$\alpha$ and neutral absorption profiles                {#sec:LISLya}

In fig. [#fig:HisLisLya], we show the neutral and ionized absorption profile as
in fig. [#fig:WithColDens] together with the profile of Ly$\alpha$ on a common
velocity scale. 

#### Figure: HisLisLya                                         {#fig:HisLisLya}

![](../Figs/LyACoverfracs.pdf){width=3.5in}

Caption: **Upper panel**: Ly$\alpha$ profile of Haro 11 C, in approximate units
of the surrounding continuum level. Full line is the measured values smoothed
by a 5 px. flat kernel; surrounding shading encloses the $\pm 1 \sigma$ error
band. **Middle panel**: Black steps show the averaged, LIS line profile,
smoothed by a 5px kernel. Surrounding gray shading denotes the $\pm 1 \sigma$
confidence band.  **Lower panel**: Same as middle panel, but for the
\ion{Si}{4} transitions. 


The picture is what we would expect from a Lyman Continuum leaker: weak neutral
and strong ionized absorption profiles, the latter almost everywhere consistent
with full coverage,  reveal a highly ionized medium in front of the central
cluster. In addition to this, it is interesting to note the close similarity in
shape between the \ion{Si}{2} and \ion{Si}{4} line profiles, indicating that
they likely originate in the same higher density regions. These regions will
likely be ionized on the side facing the central cluster, being photoionized by
this. Also the absorption feature in the Lyman $\alpha$ profile in the upper
panel seems to morphologically follow the shape of the metal lines, indicating
that radiative transfer effects are modest and most Lyman $\alpha$ photons only
scatter once, indicative of a fairly low column density of HI. We find a Lyman
$\alpha$ peak velocity of $v_{\rm peak}^{\mathrm{Ly}\alpha} = 158 \pm 1$ km s⁻¹
w.r.t the H$\alpha$-derived systemic velocity found by @Sandberg2013. This
velocity, as it is also discussed in @Verhamme2015, is just consistent with
their theoretical predictions for a density-bounded, low-column density system,
albeit on the upper limit of their allowed range. 

Lyman $\alpha$ escape is mainly governed by gas at or near systemic velocity.
In the middle panel of left subfigure of fig. [#fig:WithColDens] is shown the
best-fit column density $N_{\rm Si II}$ for each velocity bin. The value at
systemic velocity is $\log_{10}(N_{\rm Si II}) = 12.1 \pm 0.2$. Assuming a
solar Si / O ratio for Haro 11^[Given the ionization potential for Si I is
$\sim 8.5$ and for SiII is around 16, it is justified to assume the majority of
Silicon is singly-ionized FIXME Johannes, was there a reference for this?], we
can use this to estimate the column density of neutral Hydrogen in front of the
light source in the same way is in [Puschnig et al., submitted to MNRAS] as
follows.  @Guseva2012 found a metallicity of Haro 11 C of $10 + \log_{10}(O/H)
= 8.1$, corresponding to $0.25 Z_{\odot}$ when assuming solar relative
abundances. Assuming a solar Si abundance of $12 + \log_{10}(Si /H) = 7.55$
[@Asplund2009], this leads to a Si/H ratio in the neutral medium of $Si / H
\approx 1.4 \cdot 10^{-5}$, leading to a Hydrogen column density ranging from
$9.5 \cdot 10^{16} \leq N_{\rm HI} \leq 2.4 \cdot 10^{17}$ cm^-2^. Since HI
gets opaque to ionizing radiation at $\log N \sim 17$, this range is marginally
consistent with the low- optical depth, density bounded scenario, but the large
majority of the range is inconsistent with this. 

However, the pure riddled ionization bounded scenario is easily ruled out since
the Ly$\alpha$ profile does not have any appreciable component at zero
velocity. We expect a residual neutral fraction to remain in the ionized,
inter-clump phase; a fraction which has a column density high enough to block
Ly$\alpha$ efficiently at line center, but low enough to be at least part
transparent to LyC radiation. In Puschnig et al. (submitted), it was shown for
Tololo 1247, an object with very similar \ion{Si}{2} column densities, that the
sensitivity limit for \ion{Si}{2} is $\sim 1 \cdot 10^{16}$ cm^-2^, which with
the metallicity for Haro 11 C corresponds to $\sim 1.1 \cdot 10^{16}$ cm^-2^.
This leaves around 2-3 orders of magnitude in which the gas has no detectable
\ion{Si}{2} and is optically thick to Ly$\alpha$ and translucent to Lyman
continuum. 


### Neutral gas metallicity

The exact value of the metallicity is however uncertain. The values found from
nebular recombination lines by @Guseva2012 are measured mainly in the central
HII regions around the clusters; and differ by 0.2 dex between knot B and C.
The neutral, outflowing gas could be mixed, or have an unseen LOS distance
component larger than the knot separation, drawing into question which is the
better value to assume for this gas. We base our conclusions on the value found
for knot C, but note that using the metallicity of $12 + \log (O/H) = 8.3$
found for knot B by @Guseva2012, the min, median and max value of inferred
column densities are $6.03 \cdot 10^{16}$, $9.55 \cdot 10^{16}$, and $1.51
\cdot 10^{17}$ cm^-2^. However, the question of the exact metallicity of the
outflowing gas is complicated. Gas closer to the star forming regions is
expected to be more strongly enriched than gas further away, which would imply
that the HI column density is *larger* than inferred from SiII above. To this
can be added the further complication stemming from the merger event that the
galaxy is currently undergoing, which may have mixed gas of different metal
contents thoroughly. 


---
biblio-files: './main.bib' 
---
