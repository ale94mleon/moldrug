#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import moldrug, os, sys, inspect


def home(dataDir=None):
    """Return the pathname of the moldrug root directory (or a data subdirectory).
    Parameters
    ----------
    dataDir : str
        If not None, return the path to a specific data directory
    Returns
    -------
    dir : str
        The directory

    Example
    -------
    >>> from moldrug.home import home
    >>> home()                                 # doctest: +ELLIPSIS
    '.../moldrug'
    >>> home(dataDir="test-charge")                   # doctest: +ELLIPSIS
    '.../data/test-charge'
    >>> os.path.join(home(dataDir="test-charge"),"H2O.mol2")  # doctest: +ELLIPSIS
    '.../data/test-charge/H2O.mol2'
    """

    homeDir = os.path.dirname(inspect.getfile(moldrug))
    try:
        if sys._MEIPASS:
            homeDir = sys._MEIPASS
    except Exception:
        pass

    if dataDir:
        return os.path.join(homeDir, "data", dataDir)
    else:
        return homeDir