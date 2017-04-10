#!/usr/bin/env python

import pupynere as netcdf
import numpy as np
import matplotlib.pyplot as plt
import glob, sys
from mpl_toolkits.basemap import Basemap

plt.ion()

dir0 = '/media/bderembl/deneb/HYCOM/philippines/2004/'

file0 = 'state*'


allfiles1 = sorted(glob.glob(dir0 + file0));
nb_files = len(allfiles1);


f1 = netcdf.netcdf_file(allfiles1[0],'r')
lat_i = f1.variables['lat'][:].squeeze()
lon_i = f1.variables['lon'][:].squeeze()
dep_i = f1.variables['dep'][:].squeeze()
f1.close()

si_y = len(lat_i)
si_x = len(lon_i)

lon_t2,lat_t2 = np.meshgrid(lon_i,lat_i) 

lat1 = np.min(lat_i)
lat2 = np.max(lat_i)

lon1 = np.min(lon_i)
lon2 = np.max(lon_i)

# compute mean
ssh_me = np.zeros((si_y,si_x))
n_me = 0

for nfi in range(0,nb_files):
  f1 = netcdf.netcdf_file(allfiles1[nfi],'r')
  ssh = f1.variables['ssh'][:,:,:].squeeze()

  ssh_me = ssh_me + ssh
  n_me = n_me + 1

ssh_me = ssh_me/n_me


plt.figure()
for nfi in range(0,nb_files):

  f1 = netcdf.netcdf_file(allfiles1[nfi],'r')

  ssh = f1.variables['ssh'][:,:,:].squeeze()


  # lat_ts is the latitude of true scale.
  m = Basemap(projection='merc',llcrnrlat=lat1,urcrnrlat=lat2,\
              llcrnrlon=lon1,urcrnrlon=lon2,resolution='i')
  m.drawcoastlines()
  m.fillcontinents(color='gray')
  # draw parallels and meridians.
  m.drawparallels(np.arange(-90.,91.,5.),labels=[True,False,False,False])
  m.drawmeridians(np.arange(-180.,181.,5.),labels=[False,False,False,True])
  m.contour(lon_t2,lat_t2,ssh-ssh_me,colors='k',latlon=True)

  plt.savefig('figures/anim/ssh-anom' + str(nfi) + '.png')
  plt.clf()

