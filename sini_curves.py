# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 10:57:46 2016

@author: Madeleine

SINI_CURVES.py

 This function will average hit rate data for randomly-generated ultracool
 dwarf populations with transiting exoplanets at a certain distance (a/R*).
 The averaged data will be used to construct a curve to show the ideal sini
 cut-off at which we can observe the highest amount of transiting exoplanets.

Mod history:
 11/14/16 - added savefig for plots
 11/18/16 - fixed averaging of hit-rates, added 'ideal' annotation
 11/21/16 - added 'md' kwarg to plot() so we can use M-dwarf parameters
 11/22/16 - altered title of plot based on md kwarg
 12/2/16  - fixed some md parameters in plot(), added top20 kwarg
 12/8/16  - import from MDWARFS.py, added a/R* generator in plot() for md=1
 12/9/16  - changed exoplanet generator for md=1, changed title for M-dwarf
             plots
"""

import ast
import os
import numpy as np

from random_incs import *
from sini_cut import *
from mdwarfs import *


def plot(n=100, pop=100, ar=[20], numcuts=100, md=0, n_pl=1, top20=0):
    """This function produces n-number of plots of Hit Rate v. Sini Cut-off
    to show how the probability of finding a transiting exoplanet changes
    with a varying sini cut-off. The plots are saved to a folder.

    Parameters
    ----------
    *n : number
        The number of trials desired
    *pop : number
        The size of the populations that will be generated
    *ar : number
        Distance of exoplanet orbits (default [20], random when md=1)
    *numcuts :number
        The desired number of sini-cuts to generate
    *md : 0 or 1
        A kwarg to indicate whether or not to use M-dwarf parameters
    *n_pl : number
        The number of planets to be generated per star (default at 1, random
        for M-dwarfs when md=1)
    *top20 : 0 or 1
        Optional keyword that indicates whether or not to bias to the top 20%
        (for bias.bias())

    Returns
    -------
    None
    """
    hitrates = []
    ideal = []
    maxavgs = []
    cuts = gen_cuts(numcuts)
    if md == 1:
        ar = [choose_planet() for star in range(pop)]
        # ar ^^^ will be a list of lists of a/R*s (len(ar) == pop)
    for x in range(n):
        generate_pop(n=pop, ar=ar, top20=top20, md=md, n_pl=n_pl)

        # Read data from text document and then evaluate
        data = []
        sini_list = []
        siniu_list = []
        with open('gen_pop.txt') as file:
            for line in file:
                data.append(line.strip('\n'))
        transits = ast.literal_eval(data[-1])
        data = data[1:-1]
        for i in range(len(data)):
            elem = data[i].split('   ')
            sini = float(elem[0])
            siniu = float(elem[1].strip())
            sini_list += [sini]
            siniu_list += [siniu]

        evals = eval_cut2(cuts, sini_list, siniu_list, transits, top20=top20)
        hitrates += [evals[0]]
        avg_hr = [sum(e)/(x+1) for e in zip(*hitrates)]
        index = avg_hr.index(max(avg_hr))
        ideal += [cuts[index]]
        avg_ideal = np.mean(ideal)

        maxtr = evals[2]
        maxavgs += [maxtr[0]]
        avgid = np.mean(maxavgs)        # average num of transits at ideal cut
        high = 100 * max(avg_hr)
        # Average of transits spotted : stars observed at ideal sini cut-off
        obs = avgid/max(avg_hr)    # number of stars observed at ideal sini cut
        frac = 'Highest average HR = {:3.1f}% = {:3.1f}/{:3.1f}'.format(high,
                                                                        avgid,
                                                                        obs)

        plt.figure()
        plt.plot(cuts, avg_hr)
        if md == 1:
            s_type = 'M-Dwarfs'
            title = 'Hit Rate v. sin(i) Cut-Off for a population of {0} {1}\n\
Average of {2} trials'.format(pop, s_type, x+1)
        else:
            s_type = 'Ultracool Dwarfs'
            title = 'Hit Rate v. sin(i) Cut-Off for a population of {0} {1}\n\
{2} planet(s) per star, a/R* = {3}, Average of {4} trials'.format(pop, s_type,
                                                                  n_pl, ar,
                                                                  x+1)
        plt.title(title)
        plt.xlabel('sini cut-off')
        plt.xticks(np.arange(0, 1.1, 0.1))
        plt.ylabel('Hit Rate\n [# transits seen : # stars observed]')
        if max(avg_hr) < 0.40:
            plt.ylim(0, 0.40)
        else:
            plt.ylim(0, (max(avg_hr)+0.1))
        label = 'ideal sini cut-off: {:1.5f} \n'.format(avg_ideal)
        label += frac
        ax = plt.gca()
        plt.text(.5, .9,
                 label,
                 horizontalalignment='center',
                 transform=ax.transAxes)

        cwd = os.getcwd()
        if md == 0:
            path = cwd + '\\sini_curves'
        if md == 1:
            path = cwd + '\\mdwarf_curves'
        name = 'avg' + str(x+1) + '.png'
        plt.savefig(os.path.join(path, name))
        plt.close()
        status = str(x+1) + ' plots complete.'
        print(status)

    return
