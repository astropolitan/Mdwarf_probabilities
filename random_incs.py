# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 13:06:40 2016

@author: Madeleine

RANDOM_INCS.py

 This function generates random inclinations based on a probability
 distribution with set values for period and radius.

Files used:
 * generated_data.txt - w
 * all_data.txt - w
 * gen_stats.txt - a

Modification history:
 4/4/16   - revised equation for data1, added prob distrib for radius
 4/6/16   - added truth/measured/assumed elements
 4/8/16   - added code to see if my method works, o_and_i(), etc.
 4/11/16  - added bias()
 4/13/16  - write to file
 4/20/16  - fixed bug that generated negative vsini's
 4/22/16  - added SINI CUT
 6/13/16  - began revising code to fit PEP8 standards, will continue while
             reviewing code...
 6/25/16  - added clarifying comments, rearranged order of all_data data
             sublists
 6/26/16  - removed writing to 'all_data.txt', revised to fit new bias.py,
             changed how it writes to 'generated_data.txt'
 8/11/16  - added generate_pop(), gen_inc(), gen_period(), gen_radius();
             added lots of kwargs to generate_pop() for testing purposes;
             added transcribe(); removed o_and_i()
 8/15/16  - removed transcribe() and moved writing commands into
             generate_pop(), removed orbit()
 9/19/16  - added dictionary capability to generate_pop() by importing from
             STAR_DICT.py
 9/22/16  - fixed 'division by zero' error in line 287 for calculating
             percentages
 10/28/16 - temporarily(?) commenting out return statement
 11/14/16 - added 'top20' kwarg to generate_pop()
 11/21/16 - added 'md' kwargs to multiple functions in order to work with
             M-dwarf parameters
 11/22/16 - added limits to period generator
 11/28/16 - have it import * from FIND_SINI.py, started gen_planet()
 12/1/16  - continued gen_planet(), removed planet lines from generate_pop()
 12/9/16  - added planet_eval() for use with M-dwarf parameters, changed
             planet generation in generate_pop(md=1)

* LAST REVIEWED: 7/27/16
"""


import numpy as np
import math
import random
import datetime

from find_inc import find_inc
from find_sini import *
from bias import *
from star_dict import add_entry


def f(x):
    return math.sin(x)


def gen_incs(n):
    """A function to generate inclination based on a probability
    distribution.

    Parameters
    ----------
    n : number
        The size of the desired population

    Returns
    -------
    incs : list
        A list containing n generated inclinations

    Examples
    --------
    >>> gen_incs(5)
    array([ 0.89550759,  1.54404243,  0.87317641,  1.57114564,  1.71960986])
    >>> gen_incs(5)
    array([ 1.05294527,  1.48407997,  0.05819425,  2.10669321,  0.20940054])
    """
    data = np.linspace(0, 3.14, num=100)         # List of 100 nums btwn 0 & pi
    x = data
    y = [i/len(data) for i in range(len(data))]  # cumulative num / total num

    draw = np.random.rand(n)
    incs = np.interp(draw, y, x)                 # Inclinations are in radians.

    # added 12/12
    # source: http://code.activestate.com/recipes/577264-random-numbers-with-
    #           arbitrary-probability-distribu/
#    incs = []
#
#    xmin = 0.0
#    xmax = math.pi
#
#    numSteps = 10000  # bigger the better but slower!
#    ymin = f(xmin)
#    ymax = ymin
#    for i in range(numSteps):
#        x = xmin + (xmax - xmin) * float(i) / numSteps
#        y = f(x)
#        if y < ymin:
#            ymin = y
#        if y > ymax:
#            ymax = y
#
#    for i in range(n):
#        while True:
#            # generate a random number between 0 to 1
#            xr = random.random()
#            yr = random.random()
#            x = xmin + (xmax - xmin) * xr
#            y = ymin + (ymax - ymin) * yr
#            if y <= f(x):
#                incs += [xr]
#                break

    return incs


def gen_period(md=0):
    """A function that will generate a rotational period (seconds).

    Parameters
    ----------
    *md : 0 or 1
        A kwarg to indicate whether or not to use M-dwarf parameters

    Returns
    -------
    period : number
        A randomly-generated period in seconds

    Examples
    --------
    >>> gen_period()
    9590.708161387285
    >>> gen_period(md=1)
    61709.92061434599
    """
    if md == 0:
        low = 2 * 3600    # 2 hours -> seconds
        hi = 8 * 3600     # 8 hours -> seconds
    if md == 1:
        low = 2.4 * 3600       # 0.1 days -> seconds
        hi = 24 * 3600         # 1.0 day -> seconds
    period = random.uniform(low, hi)
    return period


def gen_radius(md=0):
    """Randomly generates a radius based on a Gaussian curve.

    Parameters
    ----------
    *md : 0 or 1
        A kwarg to indicate whether or not to use M-dwarf parameters

    Returns
    -------
    r : number
        One randomly generated value for radius (kilometers)

    Examples
    --------
    >>> gen_radius()
    75321.38125820422
    >>> gen_radius(md=1)
    207582.33639403537
    """
    R = 71492  # km = 1 Jupiter radius (NASA)
    if md == 1:
        Rsun = 695700  # km (source: NASA fact sheet)
        R = 0.3 * Rsun
    stdv = R * 0.1 / 3
    r = random.gauss(R, stdv)
    return r


def gen_planet(n_pl, n, ar, incs, freq=1):
    """A function to generate planets around stars.

    Parameters
    ----------
    n_pl : number
        The number of desired planets
    n : number
        The number of stars in the population
    ar : list
        A list of length n_pl with the a/R* values for the planets
    incs : list
        A list of inclinations of stars in the population
    *freq : number
        A fraction representing the intrinsic frequency of stars having an
        orbiting exoplanet (1.0 = 100% by default)

    Returns
    -------
    pl_data : list
        A list containing data about whether or not a star has orbiting planets
        and/or whether the exoplanet's transit is visible
    transits : list
        A list of indices of stars with observable transits
    """
    # Randomly assign a planet to an object based on intrinsic frequency.
    num_planets = int(freq * n)
    transits = []                            # Indices of dwarfs w/ transits
    pl_data = [[[], []] for x in range(n)]   # Full list of planetary data
    for i in range(n_pl):  # for each exoplanet...
        planets = ([1] * num_planets) + ([0] * (n - num_planets))
        random.shuffle(planets)
        for j in range(n):  # for each star in the population...
            pl_data[j][0] += [planets[j]]       # Assign a planet to each obj
            if planets[j] == 1:  # if this star has a planet...
                arcosi = abs(ar[i]*math.cos(incs[j]))
                if arcosi < 1:  # if (a/R)*cosi < 1:
                    obs = [1]
                    transits += [j]
                else:
                    obs = [0]
            else:
                obs = [0]
            pl_data[j][1] += obs
    return pl_data, transits


def planet_eval(ar, incs):
    """A function that evaluates exoplanets' transitability(?) when the
    planet's existence doesn't need to be determined. (An alternate to
    gen_planet() for use with M-dwarf data.)

    Parameters
    ----------
    ar : list
        A list of lists with the a/R* values for the planets of each star
    incs : list
        A list of inclinations of stars in the population

    Returns
    -------
    obs : list
        A list of lists representing whether or not the exoplanet's transit
        is observable
    transits : list
        A list of indices of stars with observable transits
    """
    transits = []
    obs = [[] for x in range(len(ar))]
    for i in range(len(ar)):
        for j in range(len(ar[i])):
            arcosi = abs(ar[i][j]*math.cos(incs[i]))
            if arcosi < 1:  # if (a/R)*cosi < 1:
                obs[i] += [1]
                transits += [i]
            else:
                obs[i] += [0]
    return obs, transits


def generate_pop(n=100, cut=0.95, ar=[20], freq=1, vsini_e=0.1, period_e=0.05,
                 radius_e=0.1, top20=1, md=0, n_pl=1):
    """A function that will generate a population of n stars with the
    specified parameters.

    Parameters
    ----------
    *n : number
        The number of stars in the population (100 by default)
    *cut : number
        The sini cut for this population (0.95 by default)
    *ar : list
        Semimajor axis / stellar radius ([20] by default)
    *freq : number
        A fraction representing the intrinsic frequency of stars having an
        orbiting exoplanet (1.0 = 100% by default)
    *vsini_e : number
        Rate of error on the vsini values (0.1 by default, from literature)
    *period_e : number
        The error rate assumed for the rotational period (0.05 by default,
        averaged from literature)
    *radius_e : number
        The error rate assumed for the radii (10% by default)
    *top20 : 0 or 1
        Optional keyword that indicates whether or not to bias to the top 20%
        (for bias.bias())
    *md : 0 or 1
        A kwarg to indicate whether or not to use M-dwarf parameters
    *n_pl : number
        The number of planets to be generated per star

    Returns
    -------
    all_data : list
        A list of the generated data
    transits : list
        A list of the indices of dwarfs with transiting planets
    """
#==============================================================================
#   GENERATING TRUTHS
#   -----------------
#  First, we generate a sample of stellar data. These "truths" represent a
#  dwarf's actual measurements, but this is not what we will measure.
#
#   true inc = generated
#   true period = set at 4 hours = 14400 seconds
#   true radius = randomly generated from Gaussian curve
#==============================================================================
    all_data = []

    incs = gen_incs(n)
    for i in incs:
        r = gen_radius(md)
        P = gen_period(md)
        vsin_i = (2*math.pi*r/P)*math.sin(i)
        sin_i = find_sini(vsin_i, P, r)
        all_data += [[i, vsin_i, P, round(r, 4), sin_i]]

    # We have just created a list 'all_data' containing "true" elements for
    # each generated dwarf:
    key = ['inclination (rads)', 'vsini (km/s)', 'period (seconds)',
           'radius (km)', 'sini (rads)']

    if md == 1:
        planet_data = planet_eval(ar, incs)
        observed = planet_data[0]
        transits = planet_data[1]
        for i in range(n):
            all_data[i] += [ar[i], observed]
        key += ['exoplanet a/R*(s)', 'transit seen?']
    else:
        planet_data = gen_planet(n_pl, n, ar, incs, freq)
        planets = planet_data[0]
        transits = planet_data[1]
        for i in range(n):
            all_data[i] += planets[i]
        key += ['exoplanet(s)?', 'transit seen?']

#==============================================================================
#   MEASUREMENTS/ASSUMPTIONS
#   ------------------------
#  The following code generates and calculates the measurements we would
#  obtain from observations and making assumptions.
#
#   measured vsini = 2*pi*R * sini_t + error
#                   --------
#                   period_t
#   measured period = period_t + error
#   assumed radius = 1 Jupiter radius
#
#   assumed vsini error = input percentage * vsini
#   assumed period error = input percentage * period
#   assumed radius error = input percent * radius (10% by default)
#==============================================================================
    results = 0                                 # successful transits calc'd
    for i in range(n):
        # "true" elements:
        it = all_data[i][0]                     # true inc
        vt = all_data[i][1]                     # true vsini
        pt = all_data[i][2]                     # true period
        rt = all_data[i][3]                     # true radius
        sini_t = all_data[i][4]                 # true sini

        # "measured" or assumed elements:
        ra = 71492                              # assumed radius in km
        if md == 1:
            ra = 0.3 * 695700                   # assumed radius for M-dwarfs
        ra_e = radius_e * ra                    # assumed radius error
        p_e = period_e * pt                     # period error in secs
        pm = pt + p_e                           # measured period
        vm_e = vsini_e * vt                     # measured vsini error
        vm = vt + vm_e                          # measured vsini

        # Calculate measured sini.
        sini_m = find_sini(vm, pm, ra)
        sini_u = sini_unc(vm, pm, ra, vm_e, p_e, ra_e)
        diff = abs(sini_m - sini_t)  # added abs() 10/24

        # Calculate measured inclination (in radians).
        im = find_inc(vm, pm, ra)

        if diff <= sini_u:      # This means the sini we calculated with the
            all_data[i] += [1]  # measured values is off of the true value
            results += 1        # within the calculated error.
        else:
            all_data[i] += [0]
        add = [im, vm, pm, ra, sini_m, sini_u]
        all_data[i] += add

    # *** 9/19 CHECK THIS!!  10/21 left off here
    key += ['successfully calculated inc?', 'measured inc', 'measured vsini',
            'measured period', 'measured radius', 'measured sini',
            'sini uncertainty']

    m_sinis = [all_data[x][12] for x in range(n)]
    m_sinius = [all_data[x][13] for x in range(n)]
    bias1 = bias(sini=m_sinis, sini_u=m_sinius,
                 cut=cut, transits=transits, top20=top20)
    for i in range(n):
        if i in bias1[0]:            # If index is in the top 20 selected:
            all_data[i] += [1]
        else:
            all_data[i] += [0]
        if i in bias1[1]:            # If a transit can be spotted:
            all_data[i] += [1]
        else:
            all_data[i] += [0]

#    star_catalog = {'00KEY': key}
#    for i in range(len(all_data)):
#        star_catalog = add_entry(star_catalog, i, all_data[i])

    date = 'Generated on: ' + str(datetime.datetime.now()) + '\n'

    file = open('gen_pop.txt', 'w')
    file.write(date)
    for i in range(n):
        file.write('{:06.4f}   {:06.4f} \n'.format(float(m_sinis[i]),
                                                   float(m_sinius[i])))
    file.write(str(transits))
    file.close()

    g = open('all_data.txt', 'w')
    g.write(date + str(all_data) + ',' + str(transits))
    g.close()

    return

#    f = open('generated_data.txt', 'w')
#    gen_data = [key] + [all_data]
#    f.write(date + str(gen_data))
#    f.close()

#    h = open('gen_stats.txt', 'a')
#    perc = results / n * 100
#    perc = round(perc, 4)
#    if len(transits) != 0:
#        perct = ct / len(transits) * 100
#        perct = round(perct, 4)
#    else:
#        perct = 0.0
#    h.write(str(results) + ',        ' + str(n)+',      '+str(perc)+',      '+
#            str(ct)+',     '+str(len(transits))+',      '+str(perct)+','+'\n')
#    h.close()

#    return all_data, transits, star_catalog

#sample = generate_pop()
#stardata = sample[0]
#transits = sample[1]
#cat = sample[2]
