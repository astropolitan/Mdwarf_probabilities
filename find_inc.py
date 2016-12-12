# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 12:54:27 2016

@author: Madeleine

FIND_INC.py

 Basic function for finding the inclination of a star given:
     v = vsini
     p = period
     r = radius

     * find_inc(v,p,r,deg= ) calculates the inclination of an object in either
         radians or degrees given its vsini, period, and rotation

Modification history:
 6/6/16  - added optional keyword 'deg' to find_inc(),
           added find_sini.py as a helper function
 6/10/16 - changed keyword 'deg' to 0 or 1
 6/13/16 - updated docstrings, revised code for PEP8 standards

* LAST REVIEWED: 6/10/16
"""


import math
from find_sini import find_sini


def find_inc(v, p, r, deg=0):
    """When given inputs of vsini, period, and radius, turns them into
    floats and calculates the inclination of that star (in radians) using:
        i = arcsin(v * p)
                   ------
                  (2*pi*r)

    Parameters
    ----------
    v : number
        Vsini (in kilometers/second)
    p : number
        Period (in seconds)
    r : number
        Radius (in kilometers)
    *deg : 0 or 1
        Optional keyword to convert inclinations to degrees when set to 1

    Returns
    -------
    result : number
        The inclination of the object

    Example
    -------
    >>> find_inc(23,14000,71492)
    0.79925082577082
    >>> find_inc(23,14000,71492,1)
    45.7937
    """
    v = float(v)
    p = float(p)
    r = float(r)
    sini = find_sini(v, p, r, 1)    # Using 'round_' keyword to avoid sinis > 1
    result = math.asin(sini)
    if deg == 1:
        inc = math.degrees(result)
        result = round(inc, 4)
    return result
