# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 12:41:05 2016

@author: Madeleine

DWARF_DATAGEN.py

 A function that generates random dwarf data points for evaluation and to use
 in tests of my biasing method. Changes all possible variables to see how they
 affect the "hit rate", or the number of transits spotted.

Modification history:
 8/1/16   - changed course of this program after talking to Phil on 7/27
 8/11/16  - incorporated generate_pop() from RANDOM_INCS.py, added variables
 8/15/16  - added a function specifically to work with varying uncertainties
             separately from other variables (needs lots of work...)
 9/19/16  - added star_catalog because of updates to generate_pop()
 9/22/16  - added plots to datagen()'s for-loop nest
 9/24/16  - figured out how to do subplots!
 9/29/16  - figured out how to overplot data, SWITCHED f & d for-loops!!!,
             started adjusting colors
 9/30/16  - fixed colors, added jitter so data points don't overlap
 10/3/16  - plotting 10 plots for 10 different sini cuts, added max ratio
             parameter, made all plots use same color/marker combo
 10/6/16  - added linear trendlines to plots but commented them out
 10/7/16  - 10 sublots on one figure in sini-cut loop
 10/10/16 - fixed titles on subplots, added legend, made x-labels only on
             bottom plots
 10/11/16 - fixed ratio calculation in innermost for-loop!
 10/20/16 - switched a/R* loop with cuts loop to make 'sini cut-off' the
             innermost for-loop, removed intrinsic frequency (doesn't seem to
             be helping)
 10/27/16 - removed sini cut-off for-loop, replaced with eval_cut() from
             SINI_CUT.py
 10/28/16 - continued cleaning up subplots
 10/31/16 - added maxy limit
 11/3/16  - incorporated eval_cut2()
 11/11/16 - added fractional increase plots to each subplot

* LAST REVIEWED: 10/27/16
"""

import numpy as np
import matplotlib.pyplot as plt
import ast

from random_incs import *
from sini_cut import *


#==============================================================================
#   VARIABLES
#   ---------
#   a/R*: 0 -> 100 in increments of 10
#   sini cut: 0 -> 1 in varying increments
#   number of stars in population: 50 -> 1000(?) in increments of 50
#   number of stars with exoplanets (intrinsic frequency): 0 -> 100% in
#           increments of 10%
#==============================================================================

ar = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
pop = [x * 10 for x in ar]
int_freq = [x/100 for x in ar]

errors = np.linspace(0, 1, num=21)
vsini_e = errors[:]
pd_e = errors[:]
rad_e = errors[:]


def datagen(numcuts=100, unc=0):
    """A function to generate and evaluate multiple populations of stars
    while altering certain variables.

    Parameters
    ----------
    *numcuts : number
        The number of sini cut-off values desired
    *unc : 0 or 1
        Optional keyword to decide whether to vary uncertainty values or not
    """
    cuts = gen_cuts(numcuts)
    if unc == 1:
        test_unc_prec()
    for s in range(len(pop)):           # Different population sizes.
        plt.figure(s+1, figsize=(25, 12))

        for d in range(len(ar)):        # Test different a/R* values.
            generate_pop(n=pop[s], ar=ar[d], top20=0)

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

            evals = eval_cut2(cuts, sini_list, siniu_list, transits,
                              frinc=1, top20=0)     # change top20 as necessary
            hitrate = evals[0]
            ideal = evals[1]
            num_trans = evals[2]
            frincs = evals[3]
            ideal_frac = evals[4]

            if len(transits) != 0:
                tr_ratio = float(num_trans) / len(transits)
                tr_perc = round((100 * tr_ratio), 2)
            else:
                tr_perc = 0

            plt.subplot(5, 2, d+1)
            print(cuts, hitrate)
            plt.plot(cuts, hitrate,
                     c='g',
                     label='Hit Rate\n v. Sini Cut-off')
            plt.plot(cuts, frincs,
                     label='Fractional Increase\n v. Sini Cut-off')
            if numcuts > 100:
                size = 10
                marker = '.'
            else:
                size = 16
                marker = '^'
            plt.scatter(cuts, hitrate,
                        s=size,       # marker size
                        c='magenta',
                        marker=marker)

            if d == 1:
                ax.legend(
                           bbox_to_anchor=(1.05, 1),
                           loc=2, borderaxespad=0.
                           )
            maxhit = round(max(hitrate) * 100, 2)
            if maxhit != 0:
                title = 'a/R* = {0} \n Max hit rate: {1}% of stars observed \
had a transiting exoplanet at sini cut of {2}\n {3}/{4} observable transits \
were detected ({5}%)'.format(
                            ar[d], maxhit, round(ideal, 5),
                            num_trans, len(transits), tr_perc
                            )
            else:
                title = 'a/R* = {0} \n No transits seen in this \
population.'.format(ar[d])
            ax = plt.gca()
            plt.text(.5, .65,
                     title,
                     horizontalalignment='center',
                     transform=ax.transAxes)
            plt.xticks(np.arange(0, 1.1, 0.1))
            if d == 8 or d == 9:
                plt.xlabel('sin(i) cut')
            # plt.ylabel('Hit Rate')
            plt.xlim([0, 1])
            maxy = 0.5
            if max(hitrate) > 0.45:
                maxy = max(hitrate) + 0.1
            plt.ylim([-0.02, maxy])

        subtitle = 'Population size: {0} stars \n\
Hit Rate = # of transits seen / # stars observed'.format(pop[s])
        plt.suptitle(subtitle, fontsize=15.0)

        print('Plot', s+1, 'complete.')

    #result = 'some sort of analysis'
        return


#==============================================================================
#   error in vsini :
#   error in rotational period :
#   error in radius :
#==============================================================================


def test_unc_prec():
    """A function that tests variations of vsini, period, and radius
    uncertainties. This affects the "assumed" uncertainties that are taken
    into account during star population generation in RANDOM_INCS.py's
    generate_pop().
    """
    # Precision of vsini error.
    for v in vsini_e:

        # Precision of error of rotational period.
        for p in pd_e:

            # Precision of error in radius.
            for r in rad_e:
                # generate populations of stars!
                new_pop = generate_pop(n, c, d, f, v, p, r)
                print(new_pop)
                break
