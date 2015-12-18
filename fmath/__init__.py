# -*- coding: utf-8 -*-
__author__ = 'Fei Xie & Wenzhi Mao'
__version__ = '0.0.1'

__release__ = [int(x) for x in __version__.split('.')]
del x
__all__ = []

from . import Random
from .Random import *
__all__.extend(Random.__all__)