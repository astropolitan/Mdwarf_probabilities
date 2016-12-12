# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 11:46:35 2016

@author: Madeleine

FIND_SINI.py

 Function to calculate sin(inclination) of an object given its vsini in km/s,
 period in secs, and radius in km. Also can calculate sini absolute
 uncertainty.

     * find_sini(v, p, r, round_= ) calculates the sin(inc) of an object when
         given its vsini, period, and radius
     * sini_unc(v,p,r,x,y,z,deg= ) calculates the absolute uncertainty of
         sin(i) when given vsini, period, radius, and their respective
         uncertainities. Can also convert to degrees.

Modification history:
 6/6/16   - added optional keword 'rnd' to find_sini()
 6/10/16  - changed rnd to be either 0 or 1
 6/13/16  - revised docstrings, changed code to PEP8 standards
 6/14/16  - changed keyword 'rnd' to 'round_'
 11/28/16 - added sini_unc() from UNCERTAINTY.py

* LAST REVIEWED: 6/13/16
"""

import math


def find_sini(v, p, r, round_=0):
    """Calculates the sin(i) to be used to find inclination uncertainty
        sini = v * p
              --------
              (2*pi*r)

    Parameters
    ----------
    v : number
        Vsini (in kilometers/second)
    p : number
        Period (in seconds)
    r : number
        Radius (in kilometers)
    *round_ : 0 or 1
        Optional keyword to round sini's greater than 1.0 when set to 1

    Returns
    -------
    result : number
        The sin(inclination) of the object

    Example
    -------
    >>> find_sini(35,14000,71492)
    1.0908342487974698
    >>> find_sini(35,14000,71492,1)
    1.0
    """
    arg_top = v * p
    arg_bot = 2 * math.pi * r
    result = arg_top / arg_bot
    if round_ == 1:
        if result >= 1:
            return 1.0
    return result


def sini_unc(v, p, r, x, y, z, deg=0):
    """
    A function for calculating the **ABSOLUTE** uncertainty of a calculated
    sin(i) using:
            unc/sini = sqrt( (x/v)^2 + (y/p)^2 + (z/r)^2 )
        so:
            unc = sini * sqrt( (x/v)^2 + (y/p)^2 + (z/r)^2 )

    Parameters
    ----------
    v : number
        Vsini
    p : number
        Period
    r : number
        Radius
    x : number
        Uncertainty of vsini
    y : number
        Uncertainty of the period
    z : number
        Uncertainty of the radius
    deg : 0 or 1
        Keyword that converts result to degrees when set to 1

    Returns
    -------
    res : number
        Total absolute uncertainty of the sin(i)

    Example
    -------
    >>> sini_unc(23,14000,71492,3,900,7150)
    0.12651281002313186
    >>> sini_unc(23,14000,71492,3,900,7150,1)
    7.2487
    """
    arg = ((x/v)**2) + ((y/p)**2) + ((z/r)**2)
    sq = math.sqrt(arg)
    sini = float(find_sini(v, p, r))
    result = sini*sq
    if deg == 1:       # convert to degrees
        deg_result = math.degrees(result)
        result = round(deg_result, 4)
    return result
