#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from xmitgcm import open_mdsdataset

import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
plt.ion()

dir0 = '/run/media/bderembl/girtab/southatlgyre/test_southatlgyre04/'
#dir0 = '/home/bderembl/work/MITgcm/myrun/southatlgyre/climserv/'

#ds = open_mdsdataset(dir0,iters=25920)
ds = open_mdsdataset(dir0,prefix=['U','V','T','S','Eta'])
#ds = open_mdsdataset(dir0,prefix=['Eta'])

#print(ds)

#nt = 0
nt = -1

nz = 0

plt.figure()
ax = plt.subplot(projection=ccrs.PlateCarree());
ds.Eta[nt,:,:].plot.pcolormesh('XC', 'YC', ax=ax,vmin=-1);
#ds['U'][nt,nz,:,:].plot.pcolormesh('XG', 'YC', ax=ax);
#ds['T'].where(ds.hFacC>0)[nt,nz,:,:].plot.pcolormesh('XC', 'YC', ax=ax);
ax.coastlines();
gl = ax.gridlines(draw_labels=True, alpha = 0.5, linestyle='--');
gl.xlabels_top = False
gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
