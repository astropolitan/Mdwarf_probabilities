# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 12:51:22 2016

@author: Madeleine

BIAS.py

 A helper function for random_incs.py that takes as input a list of data, a
 list of indices, and a sini cut, and outputs the number of transits that
 would be detected by my method.

     * bias(sini, sini_u, transits, cut= ) uses the value(s) of sini to bias
         a population of dwarfs or a single dwarf.
     * help_bias(sini, siniu, cut= , i= ) is a helper function to assit bias()
         in determining whether an object would have been chosen or not.

Modification history:
 6/26/16  - made bias() less monolithic by adding help_bias(), fixed code to
             adhere to PEP8 standards, revised docstrings
 10/24/16 - fixed help_bias(), fixed bias() to accomodate for different-sized
             stellar populations
 11/4/16  - revised bias()
 11/14/16 - added 'top20' keyword for bias()

* LAST REVIEWED: 10/24/16
"""

from lolimit import lolimit


def bias(sini, sini_u, cut=0.95, transits=None, top20=1):
    """Biases the population of BDs to the top candidates based *only* on
    their measured values of sini and the sini cut-off. **Chooses objects
    that will be observed, not which objects will show a transiting
    exoplanet.**

    Parameters
    ----------
    sini : number or list
        The sini(s) of the object
    sini_u : number or list
        The uncertainty/ies of the sini value(s)
    *cut : number
        The sini cut-off to be used for this bias
    *transits : list
        A list of the indices of objects with transiting exoplanets
    *top20 : 0 or 1
        Optional keyword that indicates whether or not to bias to the top 20%

    Returns
    -------
    **high_inds : list
        The indices *[of the top 20 biased] BDs (when sini is a list)
    **transited : list
        The indices of BDs *[from the top 20] with a visibly transiting planet
    **ct : number
        The total number of BDs inclined close to 90 degrees
    res : list
        A list of results for a single object's analysis

    **Returns only if original inputs are lists of data.
    """
    # For a list of sinis:
    if type(sini) == list:
        indices = []     # indices of the BDs with sini's within the parameters
        lows = []        # indices of the BDs that used lower-limit sini
        highs = []       # This will be a list of top sini values and indices.
        for i in range(len(sini)):
            res = help_bias(sini[i], sini_u[i], cut)
            if res[0] == 1:         # The lower limit was tried.
                lows += [i]
                sini[i] = res[2]    # Reset the sini to its lower limit.
            if res[1] == 1:         # This object would have been chosen.
                highs += [[sini[i], i]]
                indices += [i]
        ct = len(indices)
        if top20 == 1:
            print('Taking top 20%')     # for testing purposes
            lgth = len(sini)
            top = int(0.2 * lgth)
            highs = sorted(highs, reverse=True)   # makes list of highest sinis
            highs = highs[:top]     # The list is reduced to the top 20% sinis.
#           used_low = [str(x) + '*' for x in indices if x in lows]
        high_inds = [x[-1] for x in highs]
        transited = [x for x in high_inds if x in transits]
        return high_inds, transited, ct

    # To evaluate just one object's data:
    else:
        res = help_bias(sini, sini_u, cut)
        return res


def help_bias(sini, sini_u, cut=0.95):
    """A helper function for bias() that evaluates a sini to see whether the
    lower limit should be used and whether or not the object would have been
    biased.

    Parameters
    ----------
    sini : number
        The sini of the object
    sini_u : number
        The uncertainty of the sini value
    cut : number
        The value of the sini cut-off

    Returns
    -------
    low : 0 or 1
        Indicates whether (0) or not (1) the lower limit test was used
    biased : 0 or 1
        Indicates whether we would have chosen this object in our survey
    sini : number
        The object's sini, a different value than that which was input if
        low = 1

    Example
    -------
    >>> help_bias(0.8338721104131719, 0.11948557079946316)
    (0, 0, 0.8338721104131719)
    >>> help_bias(1.199411835930064, 0.1718637738930092)
    (1, 0, 1.199411835930064)
    """
    if sini > 1.0:              # If sini > 1.0, try the lower limit.
        nsini = float(lolimit(sini, sini_u))
        if sini != nsini:
            low = 1
        else:
            low = 0
        sini = nsini
    else:
        low = 0
    if cut <= sini <= 1.0:      # If sini is above the sini cut-off:
        biased = 1
    else:
        biased = 0
    return low, biased, sini
