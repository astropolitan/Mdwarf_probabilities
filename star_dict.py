# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 14:33:05 2016

@author: Madeleine


STAR_DICT.py

 This is a program that will create a list of stars and their data into a
 dictionary or dictionaries to make it easier to view data.

Modification history:
 9/19/16  - changed the add_star() function to add_entry(), made it less
             monolithic
"""

star_catalog = {"00key": ["vsini [km/s]", "rotational period [s]",
                "radius [km]"]}


def add_entry(catalog, name, elements):
    """A function that will add an entry a dictionary.

    Parameters
    ----------
    catalog : dictionary
        The dictionary that will be appended
    name : string
        The name of the star/dwarf
    elements : list
        A list of the elements that make up the value of the key

    Returns
    -------
    catalog : dictionary
        The updated dictionary
    """
    catalog[name] = elements
    return catalog
