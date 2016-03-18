#!/usr/bin/env python

import numpy as np

_FillValueminINF = -9.999e5
def dBinv(x):
    y = 10. ** (x / 10.)
    if np.isscalar(x):
        return 0. if (x == _FillValueminINF) else y
    else:
        return np.where(x == _FillValueminINF, 0., y)
def dB(x):
    y = 10. * np.log10(x)
    if np.isscalar(x):
        return _FillValueminINF if np.isnan(y) else y
    else:
        return np.where(np.isnan(y), _FillValueminINF, y)
