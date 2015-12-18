# -*- coding: utf-8 -*-
"""This module contains some Sobol functions.
"""

__author__ = 'Fei Xie'
__all__ = []

def getconstant(name='',**kwarg):
    """Get constants from file by name."""

    from . import __path__ as path
    from numpy import load
    from os.path import join
    from os import listdir

    path = path[0]
    if not name.endswith(".npy"):
        name+='.npy'
    if not name in listdir(path):
        raise ValueError("File {0} not exists.".format(name))
    temp = load(join(path, name))

    return temp
