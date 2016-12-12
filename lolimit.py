# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 11:15:07 2016

@author: Madeleine

LOLIMIT.py

 A helper function to test the lower limit of uncertainty for a sini to see
 whether it will yield a result.

     * lolimit(sini,siniu,inc= ) calculates the lower limit of a sini according
         to its calc'd uncertainty

Mod history:
 6/10/16  - added optional keyword 'inc' to lolimit() that calculates the
             inclination with the new low sini if set to 1.
 6/14/16  - revised code to comply with PEP8 standards
 9/1/16   - added function test_limits() from TEST_LIMITS.py
 11/28/16 - have it import * from FIND_SINI.py

* LAST REVIEWED: 10/24/16
"""

import math

from find_sini import *

# 6/10 suggestion: create way to return the low sini


def lolimit(sini, siniu, inc=0, hi=0):
    """Calculates the lower limit of sini; if the new sini is less than 1.0,
    we can calculate an inclination by setting keyword 'inc' to 1.

    Parameters
    ----------
    sini : number
        The sin(inc)
    siniu : number
        The uncertainty of the sin(i)
    *inc : 0 or 1
        An optional keyword that will calculate the inclination when set to 1
    *hi : 0 or 1
        An optional keyword that calculates the upper limit

    Returns
    -------
    *incl : number
        The inclination of the lower limit sini
    low : number
        The lower limit of sini
    sini : number
        Original sini

    Example
    -------
    >>> lolimit(1.0908342487974698,0.12651281002313186)
    0.964321438774338
    >>> lolimit(1.0908342487974698,0.12651281002313186,1)
    1.3028681155989101
    """
    low = sini - siniu
    if low <= 1.0:
        if inc == 1:
            incl = math.asin(low)
            return incl
        else:
            return low
    if hi == 1:
        high = sini + siniu
        # What was I doing here?...
    else:          # Unable to calculate inclination:
        return sini     # Used for calculations.


def test_limits(v, p, r, vu, pu, ru):
    """A function to test the limits of sini to see if it will yield a better
    number.

    Parameters
    ----------
    v : number
        Vsini
    p : number
        Period
    r : number
        Radius
    vu : number
        Uncertainty of vsini
    pu : number
        Uncertainty of the period
    ru : number
        Uncertainty of the radius

    Returns
    -------
    result :
        ...

    Examples
    --------
    >>>
    """
    sini = find_sini(v, p, r)
    unc = sini_unc(v, p, r, vu, pu, ru)

    result = lolimit(sini, unc)
    return result
