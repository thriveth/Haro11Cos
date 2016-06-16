#!/usr/bin/python
# -*- coding: utf8 -*-

import pandas as pd
import scipy as sp

colours = [  "k", "c", "m", "orange", "k", "c", "m", "orange", "k", "c", "m", "orange", "k", "c", ]
symbol  = [  "s", "s", "s", "s", "o", "o", "o", "o", "*", "*", "*", "*", "^", "^", ]
lstyle  = [  "-", "-", "-", "-", "-", "-", "-", "--", "--", "--", "--", "--", "--", "--", ]
msize   = [  10., 10., 10., 10., 11., 11., 11., 11., 15., 15., 15., 15., 10., 10., ]
lnames  = ['LARS {:02d}'.format(i+1) for i in range(14)]

lbrand = pd.DataFrame(index=range(1, 15))
lbrand['col'] = colours
lbrand['sym'] = symbol
lbrand['sty'] = lstyle
lbrand['siz'] = msize

