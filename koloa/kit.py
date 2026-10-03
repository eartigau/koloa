#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
An analysis kit of a star: everything needed to run, and change, koloa's
analysis of it on any machine, offline, as a starting point (for students
first)

    kit = build(dict(target='AN Sex', files=[dict(path='AN_SEX.rdb')],
                     root='archives', detailed=dict(dace=True)), lang='fr')
    open('kit.tar.gz', 'wb').write(kit)

The kit, a .tar.gz of one folder, koloa_<star>_analysis/:

    analysis.py      the analysis, step by step, every step explained in
                     its comments (in English or in French), calling
                     koloa's routines: the velocities read and put
                     together, a fit of their noise, the FIP in two passes,
                     the Keplerian orbits of the signals with their errors
                     and minimum masses, the acceleration of the star, the
                     folds, the FIP of what is left; its figures as PDF, its
                     numbers in results.json
    README.md        what is in the kit, how to run it, where the data come
                     from (and which may not be public)
    star.yaml        what SIMBAD, the NASA Exoplanet Archive and TESS say
                     of the star (identifiers, position, spectral type,
                     mass, periods of variability, known planets, TOIs);
                     star.json the same, for a Python without PyYAML
    data/files/      the files of velocities given
    data/archives/<star>/   what koloa.gather put on disk for the star:
                     rv/ (DACE, CARMENES DR1, VizieR), manifest.json,
                     target.json, phot/tess.csv

The settings of analysis.py are those of the page when the kit is made
(the archives ticked, the instruments left out, the trend, the FIP's
signals and sweeps, the rotation for a GP, the mass of the star).

Created on 2026-10-03

@author: artigau
"""
import io
import json
import os
import re
import tarfile
import textwrap
import time
from typing import Any, Dict, List, Optional

# =============================================================================
# Define variables
# =============================================================================
#: the comments of the script, by step and language: what each step does
#: and why (they are its teaching)
COMMENTS: Dict[str, Dict[str, str]] = {
    'intro': dict(
        en="""This script analyses the radial velocities of the star with
koloa, from the files of this folder only (no network): it can be run as
it is, and it is meant to be read and changed. Each step is a function,
and the comment above it says what it does and why.

To run it: install koloa (pip install git+https://github.com/eartigau/koloa.git),
then, in this folder: python analysis.py. Its figures go to figures/ (PDF),
its numbers to results.json and to the terminal. With the settings below,
it takes a few minutes to half an hour, most of it in the FIP.

The steps:
 1. the settings (change them here, nothing else needs to change);
 2. the star (star.yaml: its names, its mass, its known planets);
 3. the velocities: the files and the archives, put together;
 4. a first look: the velocities, instrument by instrument;
 5. the noise: a fit without planets, with outliers and jitters;
 6. the FIP: which periods hold a signal, and how sure we are;
 7. the peaks of the FIP, one per family of aliases;
 8. the orbits of the signals: Keplerians, their errors;
 9. their minimum masses;
10. the acceleration of the star (the trend);
11. the folds: the velocities at the phase of each orbit;
12. what is left: the FIP of the residuals.""",
        fr="""Ce script analyse les vitesses radiales de l'étoile avec koloa,
à partir des seuls fichiers de ce dossier (sans réseau) : il peut être
lancé tel quel, et il est fait pour être lu et modifié. Chaque étape est
une fonction, et le commentaire au-dessus dit ce qu'elle fait et pourquoi.

Pour le lancer : installez koloa
(pip install git+https://github.com/eartigau/koloa.git), puis, dans ce
dossier : python analysis.py. Ses figures vont dans figures/ (PDF), ses
nombres dans results.json et dans le terminal. Avec les réglages
ci-dessous, il prend de quelques minutes à une demi-heure, surtout dans le
FIP.

Les étapes :
 1. les réglages (changez-les ici, rien d'autre n'a à changer) ;
 2. l'étoile (star.yaml : ses noms, sa masse, ses planètes connues) ;
 3. les vitesses : les fichiers et les archives, mis ensemble ;
 4. un premier coup d'oeil : les vitesses, instrument par instrument ;
 5. le bruit : un ajustement sans planète, avec valeurs aberrantes et
    gigues ;
 6. le FIP : quelles périodes portent un signal, et avec quelle certitude ;
 7. les pics du FIP, un par famille d'alias ;
 8. les orbites des signaux : des képlériennes, leurs erreurs ;
 9. leurs masses minimales ;
10. l'accélération de l'étoile (la tendance) ;
11. les repliements : les vitesses à la phase de chaque orbite ;
12. ce qui reste : le FIP des résidus."""),
    'imports': dict(
        en="""koloa's routines used here: koloa.data (a series of velocities,
RVData, and merge to put series together), koloa.detailed (reading files
of velocities, read_files, and setting an archive apart from them,
distinct), koloa.gather (the archives on disk), koloa.published (the
velocities published on VizieR), koloa.fit (RVModel: a model of the
velocities and its fit), koloa.fip (the FIP), koloa.aliases (the aliases of
a period), koloa.kepler (Keplerian orbits), koloa.stars (minimum masses)
and koloa.secular (the acceleration of the star); koloa.utils.blas_threads
keeps each chain of the FIP to one thread. matplotlib makes the
figures ('Agg': into files, without a screen).""",
        fr="""Les routines de koloa utilisées ici : koloa.data (une série de
vitesses, RVData, et merge pour mettre des séries ensemble), koloa.detailed
(lire des fichiers de vitesses, read_files, et mettre une archive à part
d'eux, distinct), koloa.gather (les archives sur le disque),
koloa.published (les vitesses publiées sur VizieR), koloa.fit (RVModel : un
modèle des vitesses et son ajustement), koloa.fip (le FIP), koloa.aliases
(les alias d'une période), koloa.kepler (les orbites képlériennes),
koloa.stars (les masses minimales) et koloa.secular (l'accélération de
l'étoile) ; koloa.utils.blas_threads garde chaque chaîne du FIP à un fil.
matplotlib fait les figures ('Agg' : dans des fichiers, sans
écran)."""),
    'settings': dict(
        en="""The settings, those of koloa's page when this kit was made.
Change them here: the rest of the script reads them.""",
        fr="""Les réglages, ceux de la page de koloa quand ce paquet a été
fait. Changez-les ici : le reste du script les lit."""),
    'settings_files': dict(
        en="""The files of velocities (LBL .rdb, csv), relative to this folder,
each with the name of its instrument (None: koloa finds it, from the file
or its columns). An instrument is a set of velocities that share one
zero point: two pipelines of one spectrograph are two instruments.""",
        fr="""Les fichiers de vitesses (LBL .rdb, csv), relatifs à ce
dossier, chacun avec le nom de son instrument (None : koloa le trouve,
d'après le fichier ou ses colonnes). Un instrument est un ensemble de
vitesses qui partagent un point zéro : deux pipelines d'un même
spectrographe sont deux instruments."""),
    'settings_archives': dict(
        en="""The archives of the star gathered by koloa (koloa.gather): DACE
(HARPS, ESPRESSO, NIRPS...), CARMENES DR1, and the velocities published on
VizieR (Keck HIRES, the APF, the Lick Hamilton, HARPS by SERVAL, the star's
papers). True to use them.""",
        fr="""Les archives de l'étoile récupérées par koloa (koloa.gather) :
DACE (HARPS, ESPRESSO, NIRPS...), CARMENES DR1, et les vitesses publiées
sur VizieR (Keck HIRES, l'APF, le Lick Hamilton, HARPS par SERVAL, les
articles de l'étoile). True pour les utiliser."""),
    'settings_exclude': dict(
        en="""Instruments left out of the analysis, by name (['HARPS03'], say):
an old instrument with large errors, or one whose velocities look wrong.""",
        fr="""Les instruments écartés de l'analyse, par leur nom (['HARPS03'],
par exemple) : un vieil instrument aux grandes erreurs, ou un dont les
vitesses semblent fausses."""),
    'settings_trend': dict(
        en="""The trend in time fitted with the signals: 0 none, 1 a straight
line (the acceleration of the star, from a companion too far out to show
its orbit), 2 a parabola (and its change).""",
        fr="""La tendance dans le temps ajustée avec les signaux : 0 aucune,
1 une droite (l'accélération de l'étoile, due à un compagnon trop loin pour
montrer son orbite), 2 une parabole (et sa variation)."""),
    'settings_fip': dict(
        en="""The FIP. KMAX: the largest number of signals it considers at
once. NSWEEP and NBURN: the length of its Markov chains, after NBURN
sweeps thrown away while they settle; more sweeps, smoother FIPs, longer
runs. PMIN and PMAX: the periods searched [days] (PMAX None: twice the
time spanned by the velocities).""",
        fr="""Le FIP. KMAX : le plus grand nombre de signaux qu'il considère
à la fois. NSWEEP et NBURN : la longueur de ses chaînes de Markov, après
NBURN itérations jetées le temps qu'elles se stabilisent ; plus
d'itérations, des FIP plus lisses, des calculs plus longs. PMIN et PMAX :
les périodes cherchées [jours] (PMAX None : deux fois la durée couverte
par les vitesses)."""),
    'settings_threshold': dict(
        en="""A period is a signal when its FIP (of the period or any of its
aliases) is below THRESHOLD: 1 %, a chance of one in a hundred that no
signal is there.""",
        fr="""Une période est un signal quand son FIP (de la période ou d'un
de ses alias) est sous THRESHOLD : 1 %, une chance sur cent qu'aucun
signal ne s'y trouve."""),
    'settings_gp': dict(
        en="""The activity of the star (its spots, rotating with it) makes
signals too, at the rotation period and its harmonics, that are not
planets. FIP_GP says how the FIP takes it: 'banded' (the report's
default), the FIP by bands of periods, each band with a Gaussian process
of the activity only as flexible as the data ask for (koloa.bandfip);
'sho', one Gaussian process for every period, a stochastically driven
damped oscillator at ROTATION (the rotation period [days], when it is
known well); 'none', no GP (a peak at the rotation, or at its half, is
then expected: it is the star's, not a planet's).""",
        fr="""L'activité de l'étoile (ses taches, qui tournent avec elle)
fait aussi des signaux, à la période de rotation et à ses harmoniques, qui
ne sont pas des planètes. FIP_GP dit comment le FIP la prend : 'banded'
(le défaut du rapport), le FIP par bandes de périodes, chaque bande avec
un processus gaussien de l'activité seulement aussi souple que les données
le demandent (koloa.bandfip) ; 'sho', un seul processus gaussien pour
toutes les périodes, un oscillateur amorti forcé à ROTATION (la période
de rotation [jours], quand elle est bien connue) ; 'none', pas de GP (un
pic à la rotation, ou à sa moitié, est alors attendu : c'est celui de
l'étoile, pas d'une planète)."""),
    'settings_periods': dict(
        en="""Periods [days] always fitted, whatever the FIP says of them
(a planet known from its transits, say): [] for none.""",
        fr="""Des périodes [jours] toujours ajustées, quoi qu'en dise le FIP
(une planète connue par ses transits, par exemple) : [] pour aucune."""),
    'settings_nightly': dict(
        en="""True: the analysis runs on the nightly means of each
instrument (see the comment in main()); False: on every exposure.""",
        fr="""True : l'analyse se fait sur les moyennes par nuit de chaque
instrument (voir le commentaire dans main()) ; False : sur chaque
pose."""),
    'settings_mass': dict(
        en="""The mass of the star and its error [solar masses], for the
minimum masses of the planets (None: no masses).""",
        fr="""La masse de l'étoile et son erreur [masses solaires], pour les
masses minimales des planètes (None : pas de masses)."""),
    'read_star': dict(
        en="""2. The star. star.yaml holds what SIMBAD, the NASA Exoplanet
Archive and TESS say of it: its identifiers, position, spectral type and
mass, the periods of variability SIMBAD lists (a rotation, say), its known
planets and its TOIs. It is read as a dict (with PyYAML, or star.json
without it).""",
        fr="""2. L'étoile. star.yaml contient ce que SIMBAD, la NASA Exoplanet
Archive et TESS en disent : ses identifiants, sa position, son type
spectral et sa masse, les périodes de variabilité que SIMBAD liste (une
rotation, par exemple), ses planètes connues et ses TOI. Il est lu comme
un dict (avec PyYAML, ou star.json sans lui)."""),
    'read_files': dict(
        en="""3. The velocities. First the files: read_files reads each (its
time, velocity and error columns are found by name), in BJD - 2400000 and
m/s, an instrument per file or per instrument column; two files of the
same instrument become two instruments (their own zero points).""",
        fr="""3. Les vitesses. D'abord les fichiers : read_files lit chacun
(ses colonnes de temps, de vitesse et d'erreur sont trouvées par leur
nom), en BJD - 2400000 et m/s, un instrument par fichier ou par colonne
d'instrument ; deux fichiers d'un même instrument deviennent deux
instruments (leurs propres points zéro)."""),
    'read_archives': dict(
        en="""Then the archives gathered for the star: rv/all_rv.csv holds
DACE's velocities (one instrument per era of a spectrograph: HARPS03 and
HARPS15 before and after its 2015 upgrade) and CARMENES DR1's (named
CARMENES).""",
        fr="""Puis les archives récupérées pour l'étoile : rv/all_rv.csv
contient les vitesses de DACE (un instrument par époque d'un
spectrographe : HARPS03 et HARPS15 avant et après sa mise à niveau de
2015) et celles de CARMENES DR1 (nommées CARMENES)."""),
    'distinct': dict(
        en="""An archive that has an instrument of the files too (NIRPS in an
LBL file and on DACE) is another pipeline of the same spectra: distinct
names it apart (NIRPS_DACE, its own offset) and leaves out its exposures
within a minute of one of the files (the same spectrum, not counted
twice).""",
        fr="""Une archive qui a aussi un instrument des fichiers (NIRPS dans un
fichier LBL et sur DACE) est un autre pipeline des mêmes spectres :
distinct la nomme à part (NIRPS_DACE, son propre offset) et écarte ses
poses à moins d'une minute d'une pose des fichiers (le même spectre, pas
compté deux fois)."""),
    'read_vizier': dict(
        en="""Then the velocities published on VizieR (rv/published/, see
published.json for their sources): each source keeps its own offset, and
new_spectra leaves out the spectra the other series already have (a
spectrum published twice); enough leaves out an instrument with too few
velocities for an offset of its own.""",
        fr="""Puis les vitesses publiées sur VizieR (rv/published/, voir
published.json pour leurs sources) : chaque source garde son propre
offset, et new_spectra écarte les spectres que les autres séries ont déjà
(un spectre publié deux fois) ; enough écarte un instrument qui a trop peu
de vitesses pour un offset à lui."""),
    'merge': dict(
        en="""merge puts the series together: one series, its instruments
kept apart, each about its own zero point (its median: data.zero_point);
the analysis fits an offset per instrument anyway.""",
        fr="""merge met les séries ensemble : une seule série, ses
instruments gardés à part, chacun autour de son propre point zéro (sa
médiane : data.zero_point) ; l'analyse ajuste de toute façon un offset par
instrument."""),
    'exclude': dict(
        en="""The instruments of EXCLUDE left out (select keeps a part of a
series).""",
        fr="""Les instruments de EXCLUDE écartés (select garde une partie
d'une série)."""),
    'plot_series': dict(
        en="""4. A first look: every velocity, instrument by instrument, each
about its median. Look for what does not belong: an instrument off the
others, a few points far away, a trend.""",
        fr="""4. Un premier coup d'oeil : chaque vitesse, instrument par
instrument, chacun autour de sa médiane. Cherchez ce qui détonne : un
instrument décalé des autres, quelques points très loin, une tendance."""),
    'nightly': dict(
        en="""The nightly means: the exposures of one instrument in one night
averaged into one point (nightly()). A planet does not move in a night;
the noise of the star within a night (oscillations, granulation) is
averaged out, and an outlier is then a whole night.""",
        fr="""Les moyennes par nuit : les poses d'un instrument dans une nuit
moyennées en un point (nightly()). Une planète ne bouge pas en une nuit ;
le bruit de l'étoile dans une nuit (oscillations, granulation) est
moyenné, et une valeur aberrante est alors une nuit entière."""),
    'fit_noise': dict(
        en="""5. The noise. Before looking for planets, how noisy is each
instrument? RVModel fits an offset per instrument, the trend, and a
jitter per instrument (a white noise added to the errors: the activity of
the star, an instrument less stable than its errors say), with a mixture
likelihood: each point is either good, or an outlier drawn from a much
wider distribution, and the fit finds the probability of each
(fit.reliability: one minus that of being an outlier). An outlier then
weighs little, instead of being cut by hand. With several exposures per
visit (NIGHTLY = False), seq_jitter adds a jitter per visit too: a whole
visit can be off (a night of bad weather, a drift of the instrument). The orbits of PLANETS, when
given, are fitted too (their variance is then not taken for noise).
fit() finds the maximum a posteriori from a few starts.""",
        fr="""5. Le bruit. Avant de chercher des planètes, à quel point chaque
instrument est-il bruité ? RVModel ajuste un offset par instrument, la
tendance, et une gigue par instrument (un bruit blanc ajouté aux erreurs :
l'activité de l'étoile, un instrument moins stable que ne le disent ses
erreurs), avec une vraisemblance de mélange : chaque point est soit bon,
soit une valeur aberrante tirée d'une distribution bien plus large, et
l'ajustement trouve la probabilité de chacun (fit.reliability : un moins
celle d'être aberrant). Une valeur aberrante pèse alors peu, au lieu
d'être coupée à la main. Avec plusieurs poses par visite (NIGHTLY =
False), seq_jitter ajoute aussi une gigue par visite : une visite entière
peut être décalée (une nuit de mauvais temps, une dérive de
l'instrument). Les orbites des périodes données (planets) sont
ajustées aussi (leur variance n'est alors pas prise pour du bruit). fit()
trouve le maximum a posteriori à partir de quelques départs."""),
    'run_fip': dict(
        en="""6. The FIP, the false inclusion probability (Hara et al. 2022,
A&A 663, A14):
for each band of frequencies, the probability that no signal has its
period in it. A Markov chain (reversible jump) explores models with 0 to
KMAX Keplerian signals, their periods, their amplitudes, the jitters and
the outliers; the FIP of a band is the fraction of the chain with no
signal there. Unlike a periodogram's false alarm probability, it accounts
for every other signal and for the noise at once, and a FIP of 1 % means a
chance of one in a hundred that the signal is not there.""",
        fr="""6. Le FIP, la probabilité de fausse inclusion (Hara et al. 2022,
A&A 663, A14) : pour chaque bande de fréquences, la probabilité qu'aucun signal
n'y ait sa période. Une chaîne de Markov (à sauts réversibles) explore des
modèles de 0 à KMAX signaux képlériens, leurs périodes, leurs amplitudes,
les gigues et les valeurs aberrantes ; le FIP d'une bande est la fraction
de la chaîne sans signal là. Contrairement à la probabilité de fausse
alarme d'un périodogramme, il tient compte de tous les autres signaux et
du bruit à la fois, et un FIP de 1 % veut dire une chance sur cent que le
signal n'y soit pas."""),
    'inflate': dict(
        en="""The chain samples a single jitter; the instruments are not
equally noisy. inflate_to_fit adds to the errors of each instrument, in
quadrature, the jitter the fit of step 5 found for it beyond the
smallest: the chain then sees each instrument's own noise.""",
        fr="""La chaîne n'échantillonne qu'une gigue ; les instruments ne sont
pas également bruités. inflate_to_fit ajoute aux erreurs de chaque
instrument, en quadrature, la gigue que l'ajustement de l'étape 5 lui a
trouvée au-delà de la plus petite : la chaîne voit alors le bruit propre à
chaque instrument."""),
    'oafip': dict(
        en="""oafip runs two chains in processes of their own (hence the
if __name__ == '__main__' at the end of this script: the processes import
it again). outliers='both': a night or a whole visit can be an outlier.
The result: res.freq (the frequencies), res.fip (the FIP of the period
alone), res.family (of the period or any of its aliases), res.peaks (the
peaks), res.pk (the probability of k signals). blas_threads(1): one thread
of linear algebra per chain, so that the chains do not fight over the
cores.""",
        fr="""oafip lance deux chaînes dans des processus à elles (d'où le
if __name__ == '__main__' à la fin de ce script : les processus le
réimportent). outliers='both' : une nuit ou une visite entière peut être
aberrante. Le résultat : res.freq (les fréquences), res.fip (le FIP de la
période seule), res.family (de la période ou d'un de ses alias),
res.peaks (les pics), res.pk (la probabilité de k signaux).
blas_threads(1) : un fil d'algèbre linéaire par chaîne, pour que les
chaînes ne se disputent pas les coeurs."""),
    'banded': dict(
        en="""The banded FIP: the periods cut into bands; in each, a FIP
without a GP first, then with a GP of the activity more and more flexible
as long as the data ask for it (the evidence of the models, from the top
band down). as_result puts the bands together as one FIP, read as
oafip's.""",
        fr="""Le FIP par bandes : les périodes coupées en bandes ; dans
chacune, un FIP sans GP d'abord, puis avec un GP de l'activité de plus en
plus souple tant que les données le demandent (l'évidence des modèles, de
la bande du haut vers le bas). as_result met les bandes ensemble en un
seul FIP, lu comme celui de oafip."""),
    'peaks': dict(
        en="""7. The peaks, one per family of aliases. Observing from the
ground, at night, a few months a year, a period P shows at its aliases
too: 1 / (1/P + n/1 day), 1 / (1/P + n/year)... A real signal and its
aliases are one family: same_family groups them (within the width of a
peak, 1 / the time spanned), and the FIP that decides is that of the
period or any of its aliases (res.family_containing); the period alone
(res.fip_containing) tells which of the family is the true one.""",
        fr="""7. Les pics, un par famille d'alias. En observant depuis le sol,
de nuit, quelques mois par an, une période P se montre aussi à ses alias :
1 / (1/P + n/1 jour), 1 / (1/P + n/an)... Un vrai signal et ses alias
forment une famille : same_family les regroupe (à la largeur d'un pic
près, 1 / la durée couverte), et le FIP qui décide est celui de la période
ou d'un de ses alias (res.family_containing) ; celui de la période seule
(res.fip_containing) dit laquelle de la famille est la vraie."""),
    'plot_fip': dict(
        en="""The FIP as -log10: 2 is 1 %, 3 is 0.1 %. Black: the period or any
of its aliases; orange, over it: the period alone; dashed: the known
planets;
dotted: FIP = 1 %.""",
        fr="""Le FIP en -log10 : 2 est 1 %, 3 est 0,1 %. En noir : la période
ou un de ses alias ; en orange, par-dessus : la période seule ; en tirets : les planètes
connues ; en pointillé : FIP = 1 %."""),
    'second_pass': dict(
        en="""Two passes. The first fit of the noise had no planets: a planet's
variance went into the jitters, the errors were inflated too much, and
the FIP was too cautious. The second pass fits the noise again with the
signals found (FIP < THRESHOLD) and the known planets, then runs the FIP
again with errors closer to the truth.""",
        fr="""Deux passages. Le premier ajustement du bruit n'avait pas de
planète : la variance d'une planète est allée dans les gigues, les
erreurs ont été trop gonflées, et le FIP était trop prudent. Le second
passage ajuste le bruit à nouveau avec les signaux trouvés (FIP <
THRESHOLD) et les planètes connues, puis refait le FIP avec des erreurs
plus proches de la vérité."""),
    'fit_orbits': dict(
        en="""8. The orbits. RVModel now fits a Keplerian orbit per signal, its
eccentricity free and its period free within its peak, with the offsets,
the trend, the jitters and the outliers. A Keplerian is P (the period), K
(the semi-amplitude), e (the eccentricity), omega (the argument of
periastron) and tp (the time of periastron).""",
        fr="""8. Les orbites. RVModel ajuste maintenant une orbite
képlérienne par signal, son excentricité libre et sa période libre dans
son pic, avec les offsets, la tendance, les gigues et les valeurs
aberrantes. Une képlérienne, c'est P (la période), K (la
semi-amplitude), e (l'excentricité), omega (l'argument du périastre) et
tp (le temps du périastre)."""),
    'errors': dict(
        en="""Their errors: fit.cov is the Laplace covariance of the
parameters at the maximum (the curvature of the likelihood there). Its
parameters are not P, K and e themselves (they are log P, sqrt(K) cos
omega... which behave better): fit.orbits() draws parameters from the
covariance, turns each draw into an orbit (model.orbit), and gives each
quantity as (value, minus, plus), from the 16th and 84th percentiles of
the draws. An MCMC (koloa's report, --mcmc) does better when the
posterior is far from a Gaussian. tc is the time of conjunction (when a
transiting planet transits).""",
        fr="""Leurs erreurs : fit.cov est la covariance de Laplace des
paramètres au maximum (la courbure de la vraisemblance là). Ses
paramètres ne sont pas P, K et e eux-mêmes (ce sont log P, sqrt(K) cos
omega... qui se comportent mieux) : fit.orbits() tire des paramètres dans
la covariance, fait de chaque tirage une orbite (model.orbit), et donne
chaque quantité comme (valeur, moins, plus), d'après les 16e et 84e
centiles des tirages. Une MCMC (le rapport de koloa, --mcmc) fait mieux
quand la loi a posteriori est loin d'une gaussienne. tc est le temps de
la conjonction (quand une planète qui transite passe devant
l'étoile)."""),
    'masses': dict(
        en="""9. The minimum masses. The mass function, f = (m sin i)^3 /
(M + m)^2 = P K^3 (1 - e^2)^1.5 / (2 pi G), gives m sin i, the mass of the
planet times the sine of the inclination of its orbit (the velocities see
only the motion along the line of sight). minimum_mass solves it exactly,
its error from draws of K, P, e and the mass of the star, in Earth,
Neptune and Jupiter masses.""",
        fr="""9. Les masses minimales. La fonction de masse, f = (m sin i)^3 /
(M + m)^2 = P K^3 (1 - e^2)^1.5 / (2 pi G), donne m sin i, la masse de la
planète fois le sinus de l'inclinaison de son orbite (les vitesses ne
voient que le mouvement sur la ligne de visée). minimum_mass la résout
exactement, son erreur tirée de tirages de K, P, e et de la masse de
l'étoile, en masses de la Terre, de Neptune et de Jupiter."""),
    'acceleration': dict(
        en="""10. The acceleration of the star: the slope of the trend, dv/dt
[m/s/yr] (and its change with TREND = 2), with its error
(secular.acceleration). A companion too far out to complete an orbit
shows as one. The velocities also hold the perspective acceleration (the
proper motion of the star times its distance): not removed here.""",
        fr="""10. L'accélération de l'étoile : la pente de la tendance, dv/dt
[m/s/an] (et sa variation avec TREND = 2), avec son erreur
(secular.acceleration). Un compagnon trop loin pour boucler une orbite se
montre ainsi. Les vitesses contiennent aussi l'accélération de
perspective (le mouvement propre de l'étoile fois sa distance) : pas
retirée ici."""),
    'folds': dict(
        en="""11. The folds: each orbit seen alone, the velocities minus the
offsets, the trend and the other orbits, against the phase (the time
modulo the period, 0 at the conjunction, when the planet passes in front
of the star if it transits). A circle marks a night more likely an
outlier than not.""",
        fr="""11. Les repliements : chaque orbite vue seule, les vitesses
moins les offsets, la tendance et les autres orbites, en fonction de la
phase (le temps modulo la période, 0 à la conjonction, quand la planète
passe devant l'étoile si elle transite). Un cercle marque une nuit plus
probablement aberrante que non."""),
    'residuals': dict(
        en="""12. What is left: the orbits taken out of the velocities, the
FIP again. A signal the others hid may show; or the peaks of the
activity of the star.""",
        fr="""12. Ce qui reste : les orbites retirées des vitesses, le FIP de
nouveau. Un signal que les autres cachaient peut apparaître ; ou les pics
de l'activité de l'étoile."""),
    'main': dict(
        en="""The analysis, step after step. Each step prints what it found.""",
        fr="""L'analyse, étape après étape. Chaque étape affiche ce qu'elle a
trouvé."""),
    'results': dict(
        en="""The numbers, for later: results.json.""",
        fr="""Les nombres, pour plus tard : results.json."""),
    'report': dict(
        en="""koloa's full analysis, with its PDF report (the GP of the
activity, the activity indicators, the duck test of each signal, the
detection map...), from the same data: uncomment these lines.""",
        fr="""L'analyse complète de koloa, avec son rapport PDF (le GP de
l'activité, les indicateurs d'activité, le test du canard de chaque
signal, la carte de détection...), à partir des mêmes données :
décommentez ces lignes."""),
    'guard': dict(
        en="""Run main() only when the script is run, not when it is imported
(by the processes of the FIP, or by another script).""",
        fr="""Lancer main() seulement quand le script est lancé, pas quand il
est importé (par les processus du FIP, ou par un autre script)."""),
}

#: the titles of the sections of the script, in French
HEADERS_FR = {
    '0. The tools': '0. Les outils', '1. The settings': '1. Les réglages',
    '2. The star': "2. L'étoile", '3. The velocities': '3. Les vitesses',
    '4. A first look': "4. Un premier coup d'oeil",
    '5. The noise': '5. Le bruit', '6. The FIP': '6. Le FIP',
    '7. The peaks': '7. Les pics', '8. The orbits': '8. Les orbites',
    '9. The minimum masses': '9. Les masses minimales',
    '10. The acceleration': "10. L'accélération",
    '11. The folds': '11. Les repliements',
    '12. What is left': '12. Ce qui reste', 'The analysis': "L'analyse"}

#: the script, its comments marked '#> key' (replaced by COMMENTS[key]) and
#: its settings '{{NAME}}'
SCRIPT = '''#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
{{TITLE}}

{{MADE}}
"""
#> intro

# =============================================================================
# 0. The tools
# =============================================================================
#> imports
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from koloa import kepler, published, secular, stars
from koloa.aliases import same_family
from koloa.bandfip import as_result, banded_fip
from koloa.data import merge
from koloa.detailed import distinct, read_files
from koloa.fip import inflate_to_fit, oafip
from koloa.fit import RVModel
from koloa.gather import load as load_archives
from koloa.utils import blas_threads

# =============================================================================
# 1. The settings
# =============================================================================
#> settings
HERE = Path(__file__).resolve().parent
STAR = {{STAR}}
FIGURES = HERE / 'figures'
SEED = 1

#> settings_files
FILES = {{FILES}}

#> settings_archives
ARCHIVES = HERE / 'data' / 'archives' / {{ARCHIVE_FOLDER}}
USE_DACE = {{USE_DACE}}
USE_CARMENES = {{USE_CARMENES}}
USE_VIZIER = {{USE_VIZIER}}

#> settings_exclude
EXCLUDE = {{EXCLUDE}}

#> settings_trend
TREND = {{TREND}}

#> settings_fip
KMAX = {{KMAX}}
NSWEEP = {{NSWEEP}}
NBURN = {{NBURN}}
PMIN = {{PMIN}}
PMAX = {{PMAX}}

#> settings_threshold
THRESHOLD = 0.01

#> settings_gp
FIP_GP = {{FIP_GP}}
ROTATION = {{ROTATION}}

#> settings_periods
PERIODS = {{PERIODS}}

#> settings_nightly
NIGHTLY = {{NIGHTLY}}

#> settings_mass
MSTAR = {{MSTAR}}
MSTAR_ERR = {{MSTAR_ERR}}


# =============================================================================
# 2. The star
# =============================================================================
def read_star():
    #> read_star
    try:
        import yaml
        with open(HERE / 'star.yaml') as handle:
            return yaml.safe_load(handle)
    except ImportError:
        with open(HERE / 'star.json') as handle:
            return json.load(handle)


# =============================================================================
# 3. The velocities
# =============================================================================
def read_velocities():
    #> read_files
    series = []
    if FILES:
        series = read_files([str(HERE / path) for path, _ in FILES],
                            [label for _, label in FILES])
        for (path, _), part in zip(FILES, series):
            print(f'{path}: {part.n} velocities, ' + ', '.join(
                f'{inst} {np.sum(part.inst == inst)}'
                for inst in part.instruments))
    files = (series[0] if len(series) == 1 else merge(series)) \\
        if series else None

    #> read_archives
    if (USE_DACE or USE_CARMENES) and (ARCHIVES / 'rv' / 'all_rv.csv').exists():
        gathered = load_archives(str(ARCHIVES), photometry=False)['rv']
        carmenes = np.char.startswith(gathered.inst.astype(str), 'CARM')
        for use, sel, tag in ((USE_DACE, ~carmenes, 'DACE'),
                              (USE_CARMENES, carmenes, 'DR1')):
            if not use or not np.any(sel):
                continue
            part = gathered.select(sel)
            #> distinct
            if files is not None:
                part = distinct(part, files.instruments, files.time, tag)
            if part is not None:
                series.append(part)
                print(f'{tag}: {part.n} velocities, ' + ', '.join(
                    f'{inst} {np.sum(part.inst == inst)}'
                    for inst in part.instruments))

    #> read_vizier
    folder = ARCHIVES / 'rv' / 'published'
    if USE_VIZIER and folder.exists():
        pub = published.load(str(folder))
        if pub is not None:
            pub = published.enough(published.new_spectra(pub, series))
        if pub is not None:
            series.append(pub)
            print(f'VizieR: {pub.n} velocities, ' + ', '.join(
                f'{inst} {np.sum(pub.inst == inst)}'
                for inst in pub.instruments))

    if not series:
        raise SystemExit('No velocity: give a file in FILES, or set '
                         'USE_DACE, USE_CARMENES or USE_VIZIER.')
    #> merge
    data = series[0] if len(series) == 1 else merge(series, name=STAR)
    #> exclude
    if EXCLUDE:
        drop = [name.upper() for name in EXCLUDE]
        data = data.select(~np.isin(np.char.upper(data.inst.astype(str)),
                                    drop))
    print(f'Together: {data.n} velocities of {len(data.instruments)} '
          f'instruments over {data.baseline:.0f} days')
    return data


# =============================================================================
# 4. A first look
# =============================================================================
def plot_series(data, path):
    #> plot_series
    fig, ax = plt.subplots(figsize=(10, 4))
    for inst in data.instruments:
        sel = data.inst == inst
        ax.errorbar(data.time[sel], data.rv[sel] - np.median(data.rv[sel]),
                    data.err[sel], fmt='o', ms=3, lw=0.6,
                    label=f'{inst} ({int(sel.sum())})')
    ax.set_xlabel('BJD - 2400000')
    ax.set_ylabel('RV - median [m/s]')
    ax.set_title(f'{STAR}: the velocities')
    ax.legend(fontsize=7, ncol=3)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# =============================================================================
# 5. The noise
# =============================================================================
def fit_noise(nights, planets=()):
    #> fit_noise
    seq_jitter = ('instrument' if len(nights.instruments) > 1
                  and nights.nseq < nights.n else None)
    orbits = [dict(period=per, period_range=(0.98 * per, 1.02 * per))
              for per in planets]
    model = RVModel(nights, orbits, likelihood='mixture', unit='both',
                    trend=TREND, seq_jitter=seq_jitter)
    return model.fit(nstart=2, quiet=True)


# =============================================================================
# 6. The FIP
# =============================================================================
def run_fip(nights, fit, label, seed):
    #> run_fip
    #> inflate
    inflated, info = inflate_to_fit(fit)
    print(f'{label}: errors inflated by ' + ', '.join(
        f'{inst} {val:.2f}' for inst, val in info['inflation'].items())
        + ' m/s')
    #> oafip
    with blas_threads(1):
        if FIP_GP == 'banded':
            #> banded
            band = banded_fip(inflated, pmin=PMIN, pmax=PMAX, kmax=KMAX,
                              nsweep=NSWEEP, nburn=NBURN, nchains=2,
                              seed=seed, label=label, trend=TREND)
            return as_result(band)
        gp = (dict(kind='sho', period=ROTATION)
              if FIP_GP == 'sho' and ROTATION else None)
        return oafip(inflated, kmax=KMAX, outliers='both', nsweep=NSWEEP,
                     nburn=NBURN, nchains=2, seed=seed, progress=True, gp=gp,
                     trend=TREND, pmin=PMIN, pmax=PMAX, label=label,
                     nightly=False)


# =============================================================================
# 7. The peaks
# =============================================================================
def peaks_of(res, width):
    #> peaks
    found = []
    for peak in sorted(res.peaks, key=lambda peak: peak['fip']):
        period = float(peak['period'])
        if any(same_family(period, old['period'], width) for old in found):
            continue
        found.append(dict(period=period,
                          alone=float(res.fip_containing(period, width)),
                          family=float(res.family_containing(period, width))))
    return sorted(found, key=lambda peak: (peak['family'], peak['alone']))


def plot_fip(res, peaks, known, path, title):
    #> plot_fip
    period = 1.0 / np.asarray(res.freq)
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(period, -np.log10(np.clip(res.family, 1e-300, 1.0)), color='k',
            lw=1.0, label='the period or any of its aliases')
    ax.plot(period, -np.log10(np.clip(res.fip, 1e-300, 1.0)),
            color='tab:orange', lw=0.8, label='the period alone')
    ax.axhline(2, color='0.4', ls=':', lw=0.8, label='FIP = 1 %')
    for planet in known:
        ax.axvline(planet['P'], color='tab:red', ls='--', lw=0.8)
    for rank, peak in enumerate(peaks[:5]):
        ax.annotate(f'#{rank + 1}', (peak['period'],
                                     -np.log10(max(peak['family'], 1e-300))),
                    textcoords='offset points', xytext=(0, 4), ha='center',
                    fontsize=8)
    ax.set_xscale('log')
    ax.set_xlabel('period [days]')
    ax.set_ylabel('-log10 FIP')
    ax.set_title(title)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# =============================================================================
# 8. The orbits
# =============================================================================
def fit_orbits(nights, periods, width):
    #> fit_orbits
    orbits = []
    for per in periods:
        low = max(1.0 / (1.0 / per + 0.5 * width), 0.9 * per)
        high = min(1.0 / max(1.0 / per - 0.5 * width, 1e-9), 1.1 * per)
        orbits.append(dict(period=per, eccentric=True,
                           period_range=(low, high)))
    seq_jitter = ('instrument' if len(nights.instruments) > 1
                  and nights.nseq < nights.n else None)
    model = RVModel(nights, orbits, likelihood='mixture', unit='both',
                    trend=TREND, seq_jitter=seq_jitter)
    fit = model.fit(nstart=4, quiet=True)
    #> errors
    out = []
    for orbit in fit.orbits():
        out.append(dict(P=orbit['P'][0], K=orbit['K'][0], e=orbit['e'][0],
                        omega=orbit['omega'][0], tp=orbit['tp'][0],
                        tc=orbit['tc'][0],
                        P_err=0.5 * (orbit['P'][1] + orbit['P'][2]),
                        K_err=0.5 * (orbit['K'][1] + orbit['K'][2]),
                        e_err=0.5 * (orbit['e'][1] + orbit['e'][2]),
                        tc_err=0.5 * (orbit['tc'][1] + orbit['tc'][2])))
    return fit, out


# =============================================================================
# 9. The minimum masses
# =============================================================================
def masses(orbit):
    #> masses
    if not MSTAR:
        return None
    error = lambda val: float(val) if np.isfinite(val) else 0.0
    return stars.minimum_mass(orbit['K'], orbit['P'], orbit['e'], MSTAR,
                              error(orbit['K_err']), error(orbit['P_err']),
                              error(orbit['e_err']), MSTAR_ERR or 0.0)


# =============================================================================
# 10. The acceleration
# =============================================================================
def acceleration_of(fit):
    #> acceleration
    if TREND < 1:
        return None
    return secular.acceleration(fit)


# =============================================================================
# 11. The folds
# =============================================================================
def plot_folds(nights, fit, orbits, path):
    #> folds
    model = fit.model
    fig, axes = plt.subplots(1, len(orbits), figsize=(4.5 * len(orbits), 3.5),
                             squeeze=False)
    systematics = model.systematics(fit.theta)
    for ip, ax in enumerate(axes[0]):
        others = np.zeros(nights.n)
        for jp in range(len(orbits)):
            if jp != ip:
                others += model.planet_rv(fit.theta, jp)
        shown = nights.rv - systematics - others
        per, tperi, ecc, omega, amp = model.orbit(fit.theta, ip)
        tconj = kepler.tp_to_tc(tperi, per, ecc, omega)
        phase = ((nights.time - tconj) / per) % 1.0
        for inst in nights.instruments:
            sel = nights.inst == inst
            ax.errorbar(phase[sel], shown[sel], nights.err[sel], fmt='o',
                        ms=3, lw=0.5, label=str(inst))
        doubtful = fit.reliability < 0.5
        ax.plot(phase[doubtful], shown[doubtful], 'o', mfc='none', mec='k',
                ms=8)
        grid = np.linspace(0, 1, 400)
        ax.plot(grid, kepler.rv_keplerian(tconj + grid * per, per, tperi,
                                          ecc, omega, amp), 'k', lw=1.2)
        ax.set_title(f'P = {per:.4f} d, K = {amp:.2f} m/s, e = {ecc:.2f}',
                     fontsize=9)
        ax.set_xlabel('phase (0 = conjunction)')
        ax.set_ylabel('RV [m/s]')
    axes[0][0].legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# =============================================================================
# 12. What is left
# =============================================================================
def residuals_of(nights, fit, norbit):
    #> residuals
    rest = nights.select(np.ones(nights.n, dtype=bool))
    for ip in range(norbit):
        rest.rv = rest.rv - fit.model.planet_rv(fit.theta, ip)
    return rest


# =============================================================================
# The analysis
# =============================================================================
def main():
    #> main
    FIGURES.mkdir(exist_ok=True)
    star = read_star()
    known = [planet for planet in star.get('planets') or []
             if planet.get('P')]
    print(f'{STAR}: {star.get("main")}, {star.get("sptype") or "?"}, '
          f'{len(known)} known planet(s)')
    data = read_velocities()
    plot_series(data, FIGURES / 'series.pdf')

    #> nightly
    nights = data.nightly() if NIGHTLY else data
    width = 1.0 / nights.baseline
    print(f'{nights.n} ' + ('nights' if NIGHTLY else 'exposures'))

    # the first pass
    noise = fit_noise(nights)
    first = run_fip(nights, noise, 'FIP, first pass', SEED)
    peaks = peaks_of(first, width)

    #> second_pass
    periods = [peak['period'] for peak in peaks
               if peak['family'] < THRESHOLD]
    periods += [planet['P'] for planet in known
                if 1.0 < planet['P'] < nights.baseline] + list(PERIODS)
    chosen = []
    for per in sorted(periods):
        if all(abs(per / old - 1) > 0.02 for old in chosen):
            chosen.append(per)
    fip = first
    if chosen:
        noise = fit_noise(nights, chosen)
        fip = run_fip(nights, noise, 'FIP, second pass', SEED + 1)
        peaks = peaks_of(fip, width)
    plot_fip(fip, peaks, known, FIGURES / 'fip.pdf', f'{STAR}: the FIP')
    for rank, peak in enumerate(peaks[:5]):
        print(f'peak #{rank + 1}: P = {peak["period"]:.4f} d, FIP (period '
              f'or alias) {peak["family"]:.1e}, alone {peak["alone"]:.1e}')

    signals = [peak['period'] for peak in peaks
               if peak['family'] < THRESHOLD]
    signals += [per for per in PERIODS
                if all(abs(per / old - 1) > 0.02 for old in signals)]
    results = dict(star=STAR, points=int(nights.n),
                   instruments=list(map(str, nights.instruments)),
                   peaks=peaks, signals=[], acceleration=None,
                   jerk=None)
    fit = noise
    if signals:
        fit, orbits = fit_orbits(nights, signals, width)
        for orbit in orbits:
            mass = masses(orbit)
            print(f'orbit: P = {orbit["P"]:.4f} +- {orbit["P_err"]:.4f} d, '
                  f'K = {orbit["K"]:.2f} +- {orbit["K_err"]:.2f} m/s, '
                  f'e = {orbit["e"]:.2f} +- {orbit["e_err"]:.2f}'
                  + (f', m sin i = {mass["earth"][0]:.1f} Earth = '
                     f'{mass["neptune"][0]:.2f} Neptune = '
                     f'{mass["jupiter"][0]:.3f} Jupiter masses'
                     if mass else ''))
            results['signals'].append(dict(orbit, masses=mass))
        plot_folds(nights, fit, orbits, FIGURES / 'folds.pdf')
        rest = residuals_of(nights, fit, len(orbits))
        left = run_fip(rest, fit_noise(rest), 'FIP of the residuals',
                       SEED + 2)
        plot_fip(left, peaks_of(left, width), known,
                 FIGURES / 'fip_residuals.pdf',
                 f'{STAR}: the FIP of the residuals')
    else:
        print(f'no signal below a FIP of {THRESHOLD:.0%}')

    acc = acceleration_of(fit)
    if acc:
        value, low, high = acc['accel']
        print(f'acceleration of the star: dv/dt = {value:+.3f} '
              f'(-{low:.3f} +{high:.3f}) m/s/yr')
        results['acceleration'] = dict(value=value, minus=low, plus=high)
        if 'jerk' in acc:
            value, low, high = acc['jerk']
            print(f'its change: d2v/dt2 = {value:+.3f} (-{low:.3f} '
                  f'+{high:.3f}) m/s/yr^2')
            results['jerk'] = dict(value=value, minus=low, plus=high)

    #> results
    with open(HERE / 'results.json', 'w') as handle:
        json.dump(results, handle, indent=1, default=float)
    print(f'figures in {FIGURES}, numbers in results.json')

    #> report
    # from koloa.detailed import detailed_analysis
    # detailed_analysis([str(HERE / path) for path, _ in FILES] or None,
    #                   instruments=[label for _, label in FILES] or None,
    #                   target=STAR, outdir=str(HERE / 'report'),
    #                   literature=[data], dace=False, carmenes=False,
    #                   kmax=KMAX, nsweep=NSWEEP, nburn=NBURN,
    #                   trend=TREND >= 1, curvature=TREND >= 2,
    #                   rotation=ROTATION)


#> guard
if __name__ == '__main__':
    main()
'''

#: the line of the README on the archives, by language
ARCHIVE_LINE = dict(
    en="""- `data/archives/{folder}/`: the archives gathered by koloa: `rv/all_rv.csv`
  (DACE and CARMENES DR1), `rv/dace/`, `rv/carmenes/`, `rv/published/`
  (VizieR, with `published.json` naming each source), `manifest.json`
  (what each archive gave), `target.json`, `phot/tess.csv` (TESS).
""",
    fr="""- `data/archives/{folder}/` : les archives récupérées par koloa :
  `rv/all_rv.csv` (DACE et CARMENES DR1), `rv/dace/`, `rv/carmenes/`,
  `rv/published/` (VizieR, avec `published.json` qui nomme chaque source),
  `manifest.json` (ce que chaque archive a donné), `target.json`,
  `phot/tess.csv` (TESS).
""")

#: the README of the kit, by language
README = dict(
    en="""# koloa: an analysis of {star}

{made}

## What is here

- `analysis.py`: the analysis, step by step, each step explained in its
  comments. Run it with `python analysis.py` (koloa installed:
  `pip install git+https://github.com/eartigau/koloa.git`). It writes its
  figures to `figures/` (PDF) and its numbers to `results.json`.
- `star.yaml` (and `star.json`, the same): what SIMBAD, the NASA Exoplanet
  Archive and TESS say of the star.
{files_line}{archive_line}
## The data

{sources}

{warning}

## Changing it

The settings are at the top of `analysis.py` (section 1): the archives
used, the instruments left out, the trend, the FIP. Each step is a
function: change one, add one, call it from `main()`.
""",
    fr="""# koloa : une analyse de {star}

{made}

## Ce qu'il y a ici

- `analysis.py` : l'analyse, étape par étape, chaque étape expliquée dans
  ses commentaires. Lancez-la avec `python analysis.py` (koloa installé :
  `pip install git+https://github.com/eartigau/koloa.git`). Elle écrit ses
  figures dans `figures/` (PDF) et ses nombres dans `results.json`.
- `star.yaml` (et `star.json`, le même) : ce que SIMBAD, la NASA Exoplanet
  Archive et TESS disent de l'étoile.
{files_line}{archive_line}
## Les données

{sources}

{warning}

## La modifier

Les réglages sont en haut de `analysis.py` (section 1) : les archives
utilisées, les instruments écartés, la tendance, le FIP. Chaque étape est
une fonction : changez-en une, ajoutez-en une, appelez-la depuis `main()`.
""")


# =============================================================================
# Define functions
# =============================================================================
def _comment(text: str, indent: str) -> List[str]:
    """a block of text as comment lines, wrapped, at an indentation (the
    space before a French colon or semicolon kept with its word)"""
    out = []
    text = re.sub(r' ([:;?!])', '\u00a0\\1', text)
    for para in text.split('\n\n'):
        lines = para.split('\n')
        # a list (lines that start with a number or a dash) keeps its lines
        if any(re.match(r'\s*(\d+\.|-)\s', line) for line in lines[1:]):
            out += [f'{indent}# {line}'.rstrip() for line in lines]
        else:
            out += [f'{indent}# {line}'.rstrip() for line in textwrap.wrap(
                ' '.join(lines), width=76 - len(indent),
                break_on_hyphens=False)]
        out.append(f'{indent}#')
    return [line.replace('\u00a0', ' ') for line in out[:-1]]


def render(settings: Dict[str, Any], lang: str = 'en') -> str:
    """
    The script, its settings filled in and its comments in a language

    :param settings: dict, the values of the '{{NAME}}' of SCRIPT (as
                     Python source)
    :param lang: str, en or fr

    :return: str, the script
    """
    lang = 'fr' if lang == 'fr' else 'en'
    out = []
    for line in SCRIPT.split('\n'):
        mark = re.match(r'^(\s*)#> (\w+)$', line)
        if mark:
            out += _comment(COMMENTS[mark.group(2)][lang], mark.group(1))
            continue
        if lang == 'fr' and line[2:] in HEADERS_FR:
            line = '# ' + HEADERS_FR[line[2:]]
        out.append(line)
    text = '\n'.join(out)
    for key, val in settings.items():
        text = text.replace('{{' + key + '}}', str(val))
    left = re.findall(r'\{\{[A-Z_]+\}\}', text)
    if left:
        raise ValueError(f'settings not given: {left}')
    return text


def _yaml(value: Any, indent: int = 0) -> str:
    """a value as YAML (dicts, lists, numbers, strings: what star.yaml
    holds), without PyYAML"""
    pad = ' ' * indent

    def scalar(val):
        if val is None:
            return 'null'
        if isinstance(val, bool):
            return 'true' if val else 'false'
        if isinstance(val, (int, float)):
            return repr(float(val)) if isinstance(val, float) else str(val)
        return json.dumps(str(val), ensure_ascii=False)
    if isinstance(value, dict):
        if not value:
            return pad + '{}\n'
        out = ''
        for key, val in value.items():
            if isinstance(val, (dict, list)) and val:
                out += f'{pad}{key}:\n' + _yaml(val, indent + 2)
            else:
                out += f'{pad}{key}: ' + (_yaml(val).strip() if isinstance(
                    val, (dict, list)) else scalar(val)) + '\n'
        return out
    if isinstance(value, list):
        if not value:
            return pad + '[]\n'
        out = ''
        for val in value:
            if isinstance(val, dict) and val:
                body = _yaml(val, indent + 2)
                out += pad + '- ' + body[indent + 2:]
            elif isinstance(val, list) and val:
                out += pad + '-\n' + _yaml(val, indent + 2)
            else:
                out += pad + '- ' + scalar(val) + '\n'
        return out
    return pad + scalar(value) + '\n'


def star_info(target: str, root: str = '') -> Dict[str, Any]:
    """
    What the kit says of the star (star.yaml): SIMBAD's identifiers,
    position, spectral type, mass, periods of variability, CARMENES, the
    known planets and the TOIs (koloa.gui.resolve_star, the disk first)

    :return: dict
    """
    if not target.strip():
        return dict(name='', main='', note='no SIMBAD name was given')
    from koloa.gui import known_periods, resolve_star
    ident = resolve_star(target, root)
    star = ident.get('star') or {}
    keep = ('name', 'main', 'tic', 'gaia_dr3', 'hip', 'hd', 'gj', 'ra', 'dec',
            'plx', 'aliases')
    out = {key: ident.get(key) for key in keep if key in ident}
    out.update(sptype=star.get('sptype') or ident.get('sptype'),
               mass=star.get('mass'), mass_err=star.get('mass_err'),
               mass_source=star.get('source'))
    # the periods of variability SIMBAD lists (a rotation, say), each once
    out['variability'] = []
    for item in ident.get('variability') or []:
        item = {key: val for key, val in item.items() if val is not None}
        if item not in out['variability']:
            out['variability'].append(item)
    if ident.get('archive_rotation'):
        out['archive_rotation'] = ident['archive_rotation']
    if ident.get('carmenes'):
        out['carmenes'] = ident['carmenes']
    out['planets'] = known_periods(target)
    out['tois'] = ident.get('tois') or []
    out['archive_host'] = ident.get('archive_host')
    return out


def build(opts: Dict[str, Any], lang: str = 'en') -> bytes:
    """
    The kit of the page as it is: its files, the archives of its star, its
    star, the script with its settings, the README; a .tar.gz

    :param opts: dict, the fields of the page, as for the report (target,
                 files [{path, label}], root, dace, carmenes, vizier,
                 exclude, trend, curvature, kmax, nsweep, nburn, pmin, pmax,
                 periods, fip_gp, rotation, exposures), and mstar,
                 mstar_err (the star's card)
    :param lang: str, the language of the comments (en or fr)

    :return: bytes, the .tar.gz
    """
    from koloa.gather import folder_name
    from koloa.gui import DEFAULTS, _files, _number, trend_order
    target = str(opts.get('target') or '').strip()
    paths, labels = _files(opts)
    folder = folder_name(target) if target else 'series'
    top = f'koloa_{folder}_analysis'
    root = opts.get('root') or 'archives'
    archive = os.path.join(root, folder) if target else ''
    has_archive = bool(archive) and os.path.isdir(os.path.join(archive, 'rv'))
    if not paths and not has_archive:
        raise ValueError('nothing to put in the kit: no file, and no '
                         'archive gathered for the star')
    # the archives the page ticked; with no file, every one gathered
    manifest = {}
    if has_archive and os.path.exists(os.path.join(archive,
                                                   'manifest.json')):
        with open(os.path.join(archive, 'manifest.json')) as handle:
            manifest = json.load(handle).get('archives', {})
    asked = {key: bool(opts.get(key)) for key in
             ('dace', 'carmenes', 'vizier')}
    if not paths and not any(asked.values()):
        asked = dict(dace='dace' in manifest, carmenes='carmenes' in manifest,
                     vizier='published' in manifest)
    star = star_info(target, root)
    # the settings of the page
    number = lambda key, kind=float: _number(opts.get(key), kind)
    fipgp = opts.get('fip_gp', 'banded')
    fipgp = {True: 'banded', False: 'none'}.get(fipgp, fipgp) or 'banded'
    rotation = number('rotation')
    if fipgp == 'sho' and not rotation:
        fipgp = 'banded'
    mstar = _number(opts.get('mstar')) or star.get('mass')
    mstar_err = _number(opts.get('mstar_err')) or star.get('mass_err') or 0.0
    split = lambda key: [val for val in str(opts.get(key) or '')
                         .replace(',', ' ').split() if val]
    when = time.strftime('%Y-%m-%d %H:%M')
    made = (f'Fait par koloa {_version()} le {when}, depuis sa page.'
            if lang == 'fr' else
            f'Made by koloa {_version()} on {when}, from its page.')
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode='w:gz') as tar:
        def add_text(name, text):
            data = text.encode('utf-8')
            info = tarfile.TarInfo(f'{top}/{name}')
            info.size, info.mtime, info.mode = len(data), time.time(), 0o644
            tar.addfile(info, io.BytesIO(data))
        # the files of velocities
        files = []
        for rank, (path, label) in enumerate(zip(paths, labels)):
            name = f'data/files/{rank}_{os.path.basename(path)}'
            tar.add(os.path.expanduser(path), arcname=f'{top}/{name}')
            files.append((name, label or None))
        # the archives of the star: its velocities and its light curve (not
        #   the raw files of TESS)
        if has_archive:
            base = f'{top}/data/archives/{folder}'
            tar.add(os.path.join(archive, 'rv'), arcname=f'{base}/rv')
            for name in ('manifest.json', 'target.json',
                         os.path.join('phot', 'tess.csv')):
                if os.path.exists(os.path.join(archive, name)):
                    tar.add(os.path.join(archive, name),
                            arcname=f'{base}/{name}')
        # the star
        add_text('star.yaml', _yaml(star))
        add_text('star.json', json.dumps(star, indent=1, default=str))
        # the script
        fr = lang == 'fr'
        settings = dict(
            TITLE=(f"koloa : l'analyse de {target or 'une série'}" if fr
                   else f"koloa: the analysis of {target or 'a series'}"),
            MADE=made, STAR=repr(target or 'the series'),
            FILES='[' + ''.join(f'\n    ({name!r}, {label!r}),'
                                for name, label in files)
            + ('\n]' if files else ']'),
            ARCHIVE_FOLDER=repr(folder), USE_DACE=asked['dace'],
            USE_CARMENES=asked['carmenes'], USE_VIZIER=asked['vizier'],
            EXCLUDE=repr(split('exclude')), TREND=trend_order(opts),
            KMAX=number('kmax', int) or DEFAULTS['kmax'],
            NSWEEP=number('nsweep', int) or DEFAULTS['nsweep'],
            NBURN=number('nburn', int) or DEFAULTS['nburn'],
            PMIN=number('pmin') or DEFAULTS['pmin'],
            PMAX=repr(number('pmax')), FIP_GP=repr(fipgp),
            ROTATION=repr(rotation),
            PERIODS=repr([float(val) for val in split('periods')]),
            NIGHTLY=not opts.get('exposures'),
            MSTAR=repr(round(float(mstar), 4) if mstar else None),
            MSTAR_ERR=repr(round(float(mstar_err), 4) if mstar else None))
        add_text('analysis.py', render(settings, lang))
        add_text('README.md', _readme(target or 'a series', folder,
                                      manifest, files, made, lang))
    return buf.getvalue()


def _version() -> str:
    """koloa's version"""
    import koloa
    return getattr(koloa, '__version__', '')


def _readme(star: str, folder: str, manifest: Dict[str, Any], files: List,
            made: str, lang: str) -> str:
    """the README of a kit: what is in it, where its data come from, which
    may not be public"""
    fr = lang == 'fr'
    lines = []
    for name, label in files:
        lines.append(f'- `{name}`' + (f' ({label})' if label else '')
                     + (' : votre fichier' if fr else ': your file'))
    names = dict(dace='DACE', carmenes='CARMENES DR1 (Ribas et al. 2023)',
                 published='VizieR', tess='TESS (MAST)')
    for key, entry in manifest.items():
        if entry.get('status') != 'ok':
            continue
        what = ', '.join(f'{inst} {num}' for inst, num in
                         (entry.get('instruments') or {}).items())
        lines.append(f'- {names.get(key, key)}: '
                     + (what or f'{len(entry.get("sectors") or [])} sectors'))
        for src in entry.get('sources') or []:
            lines.append(f'  - {src["reference"]} ({src["catalogue"]}): '
                         f'{src["n"]}')
    private = bool(files) or any(
        'not public' in str(entry.get('key', ''))
        for entry in manifest.values())
    warning = ''
    if private:
        warning = ('**Attention** : ce paquet peut contenir des données non '
                   'publiques (vos fichiers ; DACE avec une clé). Gardez-le '
                   'dans votre équipe.' if fr else
                   '**Note**: this kit may hold data that are not public '
                   '(your files; DACE with a key). Keep it within your '
                   'team.')
    files_line = ''
    if files:
        files_line = ('- `data/files/` : les fichiers de vitesses.\n' if fr
                      else '- `data/files/`: the files of velocities.\n')
    return README['fr' if fr else 'en'].format(
        star=star, made=made, folder=folder, files_line=files_line,
        archive_line=(ARCHIVE_LINE['fr' if fr else 'en'].format(folder=folder)
                      if manifest else ''),
        sources='\n'.join(lines) or '-', warning=warning)


# =============================================================================
# End of code
# =============================================================================
