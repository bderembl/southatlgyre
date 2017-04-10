#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pupynere as netcdf

plt.ion()

dir0 = '/home/bderembl/work/MITgcm/myrun/test_philippines/run/mnc_test_0007/'
file1 = 'state.0000000000.t001.nc'


f = netcdf.netcdf_file(dir0 + file1,'r')

u = f.variables['U'][0,0,:,:]
v = f.variables['V'][0,0,:,:]

u = u[:150,:250]
v = v[:150,:250]

Q = plt.quiver( u[::5,::5], v[::5,::5])
