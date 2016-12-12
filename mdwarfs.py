# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 11:02:08 2016

@author: Madeleine

MDWARFS.py

 A function to work with M-dwarf parameters.

Modification history:
 12/8/16  - added Aurora's planet probabilites, added choose_planet() for
             generating a/R* values based on probs
 12/9/16  - fixed probabilities, changed choose_planet()
"""

import math
import random

from random_incs import *


def ar_gen(P=10):
    """This function generates a random a/R* value using parameters for
    M-dwarf stars.

    Parameters
    ----------
    P : number
        Period of the exoplanet's orbit (in days)

    Returns
    -------
    ar : number
        The ratio of semimajor axis (a) to stellar radius (R*), aka the
        exoplanet's orbital distance (in km)
    """
    P *= 3600 * 24  # converts days to seconds
    G = 6.67e-11  # Nm^2/kg^2
    Msun = 2e30  # kg
    M = 0.2 * Msun
    a = (P**2 * G * M / (4*math.pi**2))**(1/3)  # meters
    a /= 1000                                   # convert to kilometers
    #R = gen_radius(md=1)  # in km
    Rsun = 695700  # km (source: NASA fact sheet)
    R = 0.3 * Rsun
    ar = a/R
    return ar


probs = [0.000, 0.008, 0.18, 0.18, 0.36, 0.51, 0.32, 0.21, 0.42, 0.080,
         0.000, 0.006, 0.17, 0.42, 1.1, 1.4, 0.81, 1.6, 1.7, 0.16,
         0.000, 0.004, 0.23, 0.96, 2.7, 3.8, 4.6, 5.8, 4.2, 1.1,
         0.002, 0.009, 0.42, 1.8, 6.4, 9.3, 10, 12, 9.6, 4.5,
         0.061, 0.27, 1.2, 2.5, 6.7, 13, 14, 12, 8.3, 10,
         0.46, 1.4, 3.5, 5.7, 10, 13, 16, 6.4, 10, 19,
         0.40, 1.5, 4.4, 5.5, 10, 12, 11]
periods = [0.7, 1, 2, 4, 7, 12, 20, 40, 80, 100,
           0.7, 1, 2, 4, 7, 12, 20, 40, 80, 100,
           0.7, 1, 2, 4, 7, 12, 20, 40, 80, 100,
           0.7, 1, 2, 4, 7, 12, 20, 40, 80, 100,
           0.7, 1, 2, 4, 7, 12, 20, 40, 80, 100,
           0.7, 1, 2, 4, 7, 12, 20, 40, 80, 100,
           0.7, 1, 2, 4, 7, 12, 20]
orb_dist = [ar_gen(p) for p in periods]
dist_rounded = [int(round(x)) for x in orb_dist]


def choose_planet():
    """A function to return one star's planets' a/R* values by choosing from a
    probability distribution.

    Parameters
    ----------
    None

    Returns
    -------
    ar_list : list
        A list of value(s) of a/R* based on planet probability

    Examples
    --------
    >>> choose_planet()
    [43, 30]
    >>> choose_planet()
    [96, 13, 43, 30]
    """
    ar_list = []
    for x in range(len(probs)):
        num = random.uniform(0, 100)
        if num <= probs[x]:
            ar_list += [dist_rounded[x]]
    return ar_list
