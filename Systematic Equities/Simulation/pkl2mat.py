# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pickle
#import scipy.io

name = 'sectors-2018-11-02'  #probe without 1
p = open('%s.pkl'%name,'rb')
data = pickle.load(p)


#dict = {}
#dict['data'] = data
#scipy.io.savemat(name,data)

import numpy
data1 = numpy.load('hist_factor_data 2018-11-03.pkl.npy')
