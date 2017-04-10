#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from xmitgcm import open_mdsdataset

import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
plt.ion()

dir0 = '/run/media/bderembl/girtab/southatlgyre/test_southatlgyre03/'

#ds = open_mdsdataset(dir0,iters=25920)
ds = open_mdsdataset(dir0,prefix=['U','V','T','S','Eta'])

#print(ds)

plt.figure()
ax = plt.subplot(projection=ccrs.PlateCarree());
ds.Eta[0,:,:].plot.pcolormesh('XC', 'YC', ax=ax,vmin=-1);
#ds.U[0,0,:,:].plot.pcolormesh('XG', 'YC', ax=ax);
#ds['T'].where(ds.hFacC>0)[-1,0,:,:].plot.pcolormesh('XC', 'YC', ax=ax);
ax.coastlines();
gl = ax.gridlines(draw_labels=True, alpha = 0.5, linestyle='--');
gl.xlabels_top = False
gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
