# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 13:06:07 2016

@author: Madeleine

SINI_CUT.py

 A function to plot the ratio of [transits : # of BDs] as a function of
 sini cut to evaluate at which sini cut we can get the best bias of BDs with
 transits

Files used:
 * all_data.txt - r
 * sini_cuts.txt - a

Modification history:
 4/25/16  - incorporating fractional increase, reading in all_data from file
             generated by RANDOM_INCS.py
 4/27/16  - generates a plot and prints the ideal sini cut. **still appending
             elems of all_data to be 117 elements long!?**
 6/14/16  - began revising code to comply with PEP8 guidelines
 6/30/16  - revised to fit the new bias() parameters, added an else stmt to the
             frac_inc calculations
 7/27/16  - fixed first for-loop that converts data from text file into lists
             of floats, added list 'details'
 8/1/16   - added gen_cuts() to make it usable in DWARF_DATAGEN.py
 8/15/16  - revised for new version of RANDOM_INCS.py
 10/27/16 - revised for agreement with DWARF_DATAGEN.py
 11/3/16  - added eval_cut2()
 11/11/16 - added fractional increase to eval_cut2()
 11/14/16 - added 'top20' kwarg to eval_cut2()
 11/18/16 - added if/else statements to navigate top20 kwarg, made output
             'maxtr' into a list from a number
 12/2/16  - clarified variables in eval_cut2()

* LAST REVIEWED: 6/30/16
"""

import numpy as np
import matplotlib.pyplot as plt

from bias import bias


def gen_cuts(num=100):
    """Generates a list of sini cuts.

    Parameters
    ----------
    *num : number
        The number of cuts desired (100 by default)

    Returns
    -------
    cuts : list
        A list of floats between 0 and 1; the length is equal to num
    """
    cuts = np.linspace(0, 1, num=num)
    return cuts

cuts = gen_cuts()


def eval_cut2(cuts, sinis, sinius, transits, frinc=0, top20=0):
    """Less-monolithic version of eval_cut().

    Parameters
    ----------
    cuts : list
        A list of possible sini cuts
    sinis : list
        List of sini values for a population
    sinius : list
        List of sini uncertainty values for a population
    transits : list
        A list of indices that represent the stars with a transiting planet
    *frinc : 0 or 1
        Optional keyword that toggles fractional increase analysis on or off
    *top20 : 0 or 1
        Optional keyword that indicates whether or not to bias to the top 20%
        (for bias.bias())

    Returns
    -------
    ratios : list
        A list of the ratios using different sini cuts on one population
    ideal : number
        The ideal sini cut-off value
    maxtr : list
        List containing highest number of transits spotted and total number
        of objects observed in that run
    fr_incs : list
        ...
    ideal_fr : number
        ...
    """
    ideal = 0
    ideal_fr = 0
    ratios = []
    fr_incs = []
    top_cut = 0                 # Peak sini cut-off for highest planet ratio.
    top_fr = 0                  # Peak sini cut for best fractional increase.
    maxtr = [0, 0]              # Highest # of transits spotted, total observed
    for c in cuts:
        results = bias(sinis, sinius, c, transits, top20=top20)
        top20 = results[0]      # Indices of top 20% BDs that would be biased.
        toptr = results[1]      # Indices with detectable transits.
        allobs = results[2]     # Total num of BDs observed with sini cut = c
        tot_rat = float(len(transits)) / float(len(sinis))
        # 'tot_rat' is now the ratio of detected transits : entire population
        if top20 == 1:
            if len(toptr) == 0:
                tr_rat = 0
            else:  # ratio of transits detected : # stars observed
                tr_rat = float(len(toptr)) / float(len(top20))
        else:   # when we are looking at the entire population of stars
                # ratio of detectable transits : entire population size
            if len(top20) == 0:
                tr_rat = 0
            else:
                tr_rat = float(len(toptr)) / allobs
        ratios += [tr_rat]
        if tr_rat > top_cut:
            ideal = c
            top_cut = tr_rat
        if len(toptr) > maxtr[0]:
            maxtr = [len(toptr), allobs]

        # This part deals with FRACTIONAL INCREASE:
        #   by how much does the hit rate increase with a step-up in cut-off?
        if frinc == 1:
            if tr_rat == 0 or tot_rat == 0:
                frac_inc = 0
            else:
                frac_inc = tr_rat - tot_rat
            fr_incs += [frac_inc]
            if frac_inc > top_fr:
                ideal_fr = c
                top_fr = frac_inc
            return ratios, ideal, maxtr, fr_incs, ideal_fr

    return ratios, ideal, maxtr


def eval_cut(cuts, subs=0):
    """A function that evaluates the generated population and plots the
    results.

    Parameters
    ----------
    cuts : list
        A list of possible sini cuts
    *subs : 0 or 1
        1 if there will be subplots (mostly for use with DWARF_DATAGEN.py)

    Returns
    -------
    ratios : list
        A list of the ratios using different sini cuts on one population
    ideal : number
        The ideal sini cut-off value
    """
    ideal = 0
    data = []
    with open('all_data.txt') as file:
        for line in file:
            data.append(line.strip().split('],'))

    data = data[1:]                         # Cut off date line.
    data = data[0]
    data[-2] = data[-2][:-1]                # Get rid of last list bracket
    all_data = []
    for i in data[:-1]:
        elem = []
        elem = [float(x) for x in i[2:].split(',')]
        all_data.append(elem)
    tr_list = data[-1][1:-1].split(',')
    try:
        transits = [int(x) for x in tr_list]
    except:
        transits = []

    sini_list = [all_data[x][12] for x in range(len(all_data))]
    siniu_list = [all_data[x][13] for x in range(len(all_data))]

    ratios = []
    details = []
    top_cut = 0                             # Peak sini cut-off starts at zero.
    for i in cuts:
        results = bias(sini_list, siniu_list, i, transits)
        top20 = results[0]      # Indices of top 20% BDs that would be biased.
        toptr = results[1]      # Indices of top 20% BDs w detectable transits.
        alltr = results[2]      # Total num of BDs inclined close to 90 degs.
        both = [x for x in toptr if x in top20]
        tot_rat = float(len(transits)) / float(len(all_data))
        # 'tot_rat' is now the ratio of detected transits : all transits
        if len(both) == 0:
            tr_rat = 0
        else:
            tr_rat = float(len(both)) / float(len(top20))
        # Check: what is the purpose of frac_inc ??
        if tr_rat == 0 or tot_rat == 0:
            frac_inc = 0
        else:
            frac_inc = tr_rat - tot_rat
        ratios += [tr_rat]
        #ratios += [frac_inc]
        if frac_inc > top_cut:
            ideal = i
            top_cut = frac_inc
        details += [[i, top20, toptr, alltr, both]]

    if subs == 0:
        plt.plot(cuts, ratios)
        print(ideal)

        f = open('sini_cuts.txt', 'a')
        f.write(str(ideal)+'\n')
        f.close()

    return ratios, ideal
