import numpy as np
import struct
import matplotlib.pyplot as plt
import scipy.io.netcdf as netcdf
from scipy import interpolate
import intergrid
import glob,os,re
from scipy import ndimage as nd
from subprocess import call
from scipy.ndimage.filters import gaussian_filter
import f90nml

plt.ion()

binprec = '>f4'

outputdir1 = '/home/bderembl/work/MITgcm/myrun/test_southatlgyre/input/files/'
outputdir2 = '/home/bderembl/work/MITgcm/myrun/test_southatlgyre/input/figures/'

# initial condition
dir0_i = '/home/bderembl/work/MITgcm/myrun/test_southatlgyre/input/files/2000/'
file1 = 'state-20000101.nc'

# bc
dir0_bc = '/home/bderembl/work/MITgcm/myrun/test_southatlgyre/input/files/bc/2000/'

# era data
dir0_era = '/home/bderembl/work/MITgcm/myrun/test_southatlgyre/input/files/era/'

flag_plot = 1
flag_atmos = 1

if not os.path.exists(outputdir1):
  os.makedirs(outputdir1)
if not os.path.exists(outputdir2):
  os.makedirs(outputdir2)


#% ================ GRID =========================================                                  
# dx = 1/10.0
# dy = 1/10.0

#si_x = 960
#si_y = 280
si_x = 1120
si_y = 420

# for mitgcm config
npx = 16
npy = 4

if int(si_x/npx) != si_x/npx:
  si_x = npx*int(si_x/npx)
  print('adjusting si_x to {0}'.format(si_x))
if int(si_y/npy) != si_y/npy:
  si_y = npy*int(si_y/npy)
  print('adjusting si_y to {0}'.format(si_y))

  
#nb of layers
si_z = 60

# max depth
Lz = 5000.0


lat1 = -44
lat2 = -15

lon1 = -66
lon2 = 28

R    = 6400e3
deg1 = 2*np.pi*R/360

# lat = np.arange(lat1,lat2,dy)
# lon = np.arange(lon1,lon2,dx)
# si_x = len(lon)
# si_y = len(lat)

lat = np.linspace(lat1,lat2,si_y)
lon = np.linspace(lon1,lon2,si_x)
lon_g,lat_g = np.meshgrid(lon,lat) 


xx = np.zeros((si_x*si_y,2))
xx[:,0] = lat_g.flat
xx[:,1] = lon_g.flat

dx = lon[1] - lon[0]
dy = lat[1] - lat[0]

dx_m = dx*deg1

print('dx = {0:.1f} m'.format(dx_m))
print('dy = {0:.1f} m'.format(dy*deg1))

# generate vertical grid
si_z1 = si_z + 1


# xf is % of grid points
xf = [0, 0.4, 0.6, 0.8, 0.9, 1]
# yf is % of thickness
yf = [0, 0.05, 0.11, 0.21, 0.4, 1]

hh = np.linspace(0,1,si_z1)
zz = Lz*np.interp(hh,xf,yf)

# smooth
nc = int(si_z/10)
if nc % 2 == 0:
  nc = nc + 1
zz2 = np.convolve(zz, np.ones((nc,))/nc, mode='valid')

zz[int((nc-1)/2):int(-(nc-1)/2)] = zz2

if flag_plot:
  plt.figure()
  plt.plot(hh,zz/Lz,'k')
  plt.plot(hh,hh,'k--')
  plt.plot(xf,yf,'.')
  plt.savefig(outputdir2 + 'vert_res.png')
  plt.close()

dz1 = np.diff(zz)

iz = np.argmin(np.abs(zz-1000.0))

print ('min dz: ', np.min(dz1))
print ('max dz: ', np.max(dz1))
print ('nb layers above 1000m:', iz, '/', si_z)

if np.sum(dz1 < 0) > 0:
  print ('you need you change the polynomial fit!!')

dep2 = zz[0:-1] + 0.5*dz1


dx1 = dx*np.ones((si_x))
dy1 = dy*np.ones((si_y))

dx1.astype(binprec).tofile(outputdir1 + 'dx.box')
dy1.astype(binprec).tofile(outputdir1 + 'dy.box')
dz1.astype(binprec).tofile(outputdir1 + 'dz.box')

print('config: si_x = {0}, si_y = {1}, si_z = {2}'.format(si_x,si_y,si_z))

# ===== TOPOGRAPHY =======

# file_topo = '/home/bderembl/work/data/topo/etopo5.nc'
# lat_name = 'topo_lat'
# lon_name = 'topo_lon'
# topo_name = 'topo'

file_topo = '/home/bderembl/work/data/topo/etopo1.nc'
lat_name = 'y'
lon_name = 'x'
topo_name = 'z'

ft = netcdf.netcdf_file(file_topo,'r')

lat_t = ft.variables[lat_name][:].copy().squeeze()
lon_t = ft.variables[lon_name][:].copy().squeeze()


y1 = np.argmin(np.abs(lat_t-lat1))-10
y2 = np.argmin(np.abs(lat_t-lat2))+10

x1 = np.argmin(np.abs(lon_t-lon1))-10
x2 = np.argmin(np.abs(lon_t-lon2))+10

topo =  ft.variables[topo_name][y1:y2,x1:x2].copy().squeeze()

lat_t = lat_t[y1:y2]
lon_t = lon_t[x1:x2]

lon_t2,lat_t2 = np.meshgrid(lon_t,lat_t) 

func_topo = interpolate.interp2d(lon_t, lat_t, topo, kind='linear')

topo_newgrid = func_topo(lon,lat)

topo_newgrid = np.where(topo_newgrid > 0, 0, topo_newgrid)

topo_smooth = gaussian_filter(topo_newgrid,1)

topo_newgrid.astype(binprec).tofile(outputdir1 + 'htopo.box')
#topo_smooth.astype(binprec).tofile(outputdir1 + 'htopo.box')

mask = np.zeros((si_z,si_y,si_x))

for nz in range(0,si_z):
  mask[nz,:,:] = np.where(-zz[nz]<topo_newgrid,0,1)


if flag_plot:
  plt.figure()
  plt.imshow(topo_smooth[::-1,:],interpolation='nearest')
  plt.colorbar()
  plt.savefig(outputdir2 + 'topo.png')
  plt.clf()

#% ================ initial condition =============================  


f1 = netcdf.netcdf_file(dir0_i + file1,'r')

lat_i = f1.variables['lat'][:].copy().squeeze()
lon_i = f1.variables['lon'][:].copy().squeeze()
dep_i = f1.variables['dep'][:].copy().squeeze()

si_zi = len(dep_i)

lon_i2,lat_i2 = np.meshgrid(lon_i,lat_i) 

u_init = f1.variables['u_vel'][:,:,:].copy().squeeze()
v_init = f1.variables['v_vel'][:,:,:].copy().squeeze()
t_init = f1.variables['theta'][:,:,:].copy().squeeze()
s_init = f1.variables['salt'][:,:,:].copy().squeeze()
e_init = f1.variables['ssh'][:,:].copy().squeeze()

u_init = 1.0*u_init
v_init = 1.0*v_init
t_init = 1.0*t_init
s_init = 1.0*s_init
e_init = 1.0*e_init

# nanfix done after interpolation
# u_init[np.isnan(u_init)] = 0.0
# v_init[np.isnan(v_init)] = 0.0
# #FIXME
# t_init[np.isnan(t_init)] = 0.0
# s_init[np.isnan(s_init)] = 0.0

e_init[np.isnan(e_init)] = 0.0

u_out1 = np.zeros((si_zi,si_y,si_x))
v_out1 = np.zeros((si_zi,si_y,si_x))
t_out1 = np.zeros((si_zi,si_y,si_x))
s_out1 = np.zeros((si_zi,si_y,si_x))
e_out  = np.zeros((si_y,si_x))


lo = np.array([ np.min(lat_i), np.min(lon_i)])  # lowest lat, lowest lon
hi = np.array([ np.max(lat_i), np.max(lon_i) ])   # highest lat, highest lon

interfunc = intergrid.Intergrid( e_init, lo=lo, hi=hi , verbose=0, order=3)
e_out[:,:] = interfunc.at( xx).reshape(si_y,si_x)


for nz in range(0,si_zi):

  interfunc = intergrid.Intergrid( u_init[nz,:,:], lo=lo, hi=hi , verbose=0, order=3)
  u_out1[nz,:,:] = interfunc.at( xx).reshape(si_y,si_x)

  interfunc = intergrid.Intergrid( v_init[nz,:,:], lo=lo, hi=hi , verbose=0, order=3)
  v_out1[nz,:,:] = interfunc.at( xx).reshape(si_y,si_x)

  interfunc = intergrid.Intergrid( t_init[nz,:,:], lo=lo, hi=hi , verbose=0, order=3)
  t_out1[nz,:,:] = interfunc.at( xx).reshape(si_y,si_x)

  interfunc = intergrid.Intergrid( s_init[nz,:,:], lo=lo, hi=hi , verbose=0, order=3)
  s_out1[nz,:,:] = interfunc.at( xx).reshape(si_y,si_x)


# vertical interpolation

u_out = np.zeros((si_z,si_y,si_x))
v_out = np.zeros((si_z,si_y,si_x))
t_out = np.zeros((si_z,si_y,si_x))
s_out = np.zeros((si_z,si_y,si_x))

#store z weights for BC interpolation
wei_z = np.zeros((si_z,3))

for nz in range(0,si_z):
  iz1 = np.argsort(np.abs(dep_i - dep2[nz]))
  nz1 = np.min(iz1[0:2])
  nz2 = np.max(iz1[0:2])
  if dep2[nz] > dep_i[nz1] and dep2[nz] < dep_i[nz2]:
    a1 = (dep2[nz] - dep_i[nz1])/(dep_i[nz2] - dep_i[nz1])
    u_out[nz,:,:] = (1-a1)*u_out1[nz1,:,:] + a1*u_out1[nz2,:,:]
    v_out[nz,:,:] = (1-a1)*v_out1[nz1,:,:] + a1*v_out1[nz2,:,:]
    t_out[nz,:,:] = (1-a1)*t_out1[nz1,:,:] + a1*t_out1[nz2,:,:]
    s_out[nz,:,:] = (1-a1)*s_out1[nz1,:,:] + a1*s_out1[nz2,:,:]
    
    wei_z[nz,0] = a1
    wei_z[nz,1] = nz1
    wei_z[nz,2] = nz2
  else:
    u_out[nz,:,:] = u_out1[iz1[0],:,:]
    v_out[nz,:,:] = v_out1[iz1[0],:,:]
    t_out[nz,:,:] = t_out1[iz1[0],:,:]
    s_out[nz,:,:] = s_out1[iz1[0],:,:]

    wei_z[nz,0] = 0.0
    wei_z[nz,1] = iz1[0]
    wei_z[nz,2] = 0


ind = nd.distance_transform_edt(np.isnan(u_out), return_distances=False, return_indices=True)
u_out = u_out[tuple(ind)]
ind = nd.distance_transform_edt(np.isnan(v_out), return_distances=False, return_indices=True)
v_out = v_out[tuple(ind)]
ind = nd.distance_transform_edt(np.isnan(t_out), return_distances=False, return_indices=True)
t_out = t_out[tuple(ind)]
ind = nd.distance_transform_edt(np.isnan(s_out), return_distances=False, return_indices=True)
s_out = s_out[tuple(ind)]

u_out.astype(binprec).tofile(outputdir1 + 'uinit.box')
v_out.astype(binprec).tofile(outputdir1 + 'vinit.box')
t_out.astype(binprec).tofile(outputdir1 + 'tinit.box')
s_out.astype(binprec).tofile(outputdir1 + 'sinit.box')
e_out.astype(binprec).tofile(outputdir1 + 'einit.box')

if flag_plot:
  plt.subplot(221)
  plt.imshow(u_out[0,::-1,:],interpolation='nearest'); plt.colorbar()
  plt.subplot(222)
  plt.imshow(v_out[0,::-1,:],interpolation='nearest'); plt.colorbar()
  plt.subplot(223)
  plt.imshow(t_out[0,::-1,:],interpolation='nearest'); plt.colorbar()
  plt.subplot(224)
  plt.imshow(s_out[0,::-1,:],interpolation='nearest'); plt.colorbar()
  plt.savefig(outputdir2 + 'init_state_top.png')
  plt.clf()


  iz2 = int(np.floor(si_z/2))
  plt.subplot(221)
  plt.imshow(u_out[iz2,::-1,:],interpolation='nearest'); plt.colorbar()
  plt.subplot(222)
  plt.imshow(v_out[iz2,::-1,:],interpolation='nearest'); plt.colorbar()
  plt.subplot(223)
  plt.imshow(t_out[iz2,::-1,:],interpolation='nearest'); plt.colorbar()
  plt.subplot(224)
  plt.imshow(s_out[iz2,::-1,:],interpolation='nearest'); plt.colorbar()
  plt.savefig(outputdir2 + 'init_state_mid.png')
  plt.clf()

  plt.subplot(221)
  plt.imshow(u_out[-1,::-1,:],interpolation='nearest'); plt.colorbar()
  plt.subplot(222)
  plt.imshow(v_out[-1,::-1,:],interpolation='nearest'); plt.colorbar()
  plt.subplot(223)
  plt.imshow(t_out[-1,::-1,:],interpolation='nearest'); plt.colorbar()
  plt.subplot(224)
  plt.imshow(s_out[-1,::-1,:],interpolation='nearest'); plt.colorbar()
  plt.savefig(outputdir2 + 'init_state_bot.png')
  plt.clf()


#============= boundary conditions ============================

file0 = 'bc_*'

allfiles1 = sorted(glob.glob(dir0_bc + file0));
nb_files = len(allfiles1);


f1 = netcdf.netcdf_file(allfiles1[0],'r')
lat_i = f1.variables['lat'][:].copy().squeeze()
lon_i = f1.variables['lon'][:].copy().squeeze()
dep_i = f1.variables['dep'][:].copy().squeeze()
f1.close()

vars2 = ['u_vel', 'v_vel', 'theta', 'salt']
nvars = len(vars2)


u_s_out = np.zeros((nb_files,si_z,si_x))
v_s_out = np.zeros((nb_files,si_z,si_x))
t_s_out = np.zeros((nb_files,si_z,si_x))
s_s_out = np.zeros((nb_files,si_z,si_x))
 
u_n_out = np.zeros((nb_files,si_z,si_x))
v_n_out = np.zeros((nb_files,si_z,si_x))
t_n_out = np.zeros((nb_files,si_z,si_x))
s_n_out = np.zeros((nb_files,si_z,si_x))
 
u_w_out = np.zeros((nb_files,si_z,si_y))
v_w_out = np.zeros((nb_files,si_z,si_y))
t_w_out = np.zeros((nb_files,si_z,si_y))
s_w_out = np.zeros((nb_files,si_z,si_y))
 
u_e_out = np.zeros((nb_files,si_z,si_y))
v_e_out = np.zeros((nb_files,si_z,si_y))
t_e_out = np.zeros((nb_files,si_z,si_y))
s_e_out = np.zeros((nb_files,si_z,si_y))

psi_s_out1 = np.zeros((si_zi,si_x))
psi_n_out1 = np.zeros((si_zi,si_x))
psi_w_out1 = np.zeros((si_zi,si_y))
psi_e_out1 = np.zeros((si_zi,si_y))

for nfi in range(0,nb_files):
  file1 = allfiles1[nfi]
  f = netcdf.netcdf_file(file1,'r')

  nva = -1
  for myvar in vars2:
    nva = nva + 1

    psi_s = f.variables[myvar + '_s'][:,:,:].copy().squeeze()
    psi_n = f.variables[myvar + '_n'][:,:,:].copy().squeeze()
    psi_w = f.variables[myvar + '_w'][:,:,:].copy().squeeze()
    psi_e = f.variables[myvar + '_e'][:,:,:].copy().squeeze()


    #South
    finterp = interpolate.interp1d(lon_i,psi_s)
    psi_s_out1 = finterp(lon)
    #takes care of nans
    ind = nd.distance_transform_edt(np.isnan(psi_s_out1), return_distances=False, return_indices=True)
    psi_s_out1 = psi_s_out1[tuple(ind)]

    #North
    finterp = interpolate.interp1d(lon_i,psi_n)
    psi_n_out1 = finterp(lon)
    #takes care of nans
    ind = nd.distance_transform_edt(np.isnan(psi_n_out1), return_distances=False, return_indices=True)
    psi_n_out1 = psi_n_out1[tuple(ind)]

    #West
    finterp = interpolate.interp1d(lat_i,psi_w)
    psi_w_out1 = finterp(lat)
    #takes care of nans
    ind = nd.distance_transform_edt(np.isnan(psi_w_out1), return_distances=False, return_indices=True)
    psi_w_out1 = psi_w_out1[tuple(ind)]

    #East
    finterp = interpolate.interp1d(lat_i,psi_e)
    psi_e_out1 = finterp(lat)
    #takes care of nans
    ind = nd.distance_transform_edt(np.isnan(psi_e_out1), return_distances=False, return_indices=True)
    psi_e_out1 = psi_e_out1[tuple(ind)]

    if myvar == 'u_vel':
      u_s_out[nfi,:,:] = (1 - wei_z[:,0].reshape(si_z,1))*psi_s_out1[wei_z[:,1].astype(int),:] + wei_z[:,0].reshape(si_z,1)*psi_s_out1[wei_z[:,2].astype(int),:]
      u_n_out[nfi,:,:] = (1 - wei_z[:,0].reshape(si_z,1))*psi_n_out1[wei_z[:,1].astype(int),:] + wei_z[:,0].reshape(si_z,1)*psi_n_out1[wei_z[:,2].astype(int),:]
      u_w_out[nfi,:,:] = (1 - wei_z[:,0].reshape(si_z,1))*psi_w_out1[wei_z[:,1].astype(int),:] + wei_z[:,0].reshape(si_z,1)*psi_w_out1[wei_z[:,2].astype(int),:]
      u_e_out[nfi,:,:] = (1 - wei_z[:,0].reshape(si_z,1))*psi_e_out1[wei_z[:,1].astype(int),:] + wei_z[:,0].reshape(si_z,1)*psi_e_out1[wei_z[:,2].astype(int),:]
    elif myvar == 'v_vel':              
      v_s_out[nfi,:,:] = (1 - wei_z[:,0].reshape(si_z,1))*psi_s_out1[wei_z[:,1].astype(int),:] + wei_z[:,0].reshape(si_z,1)*psi_s_out1[wei_z[:,2].astype(int),:]
      v_n_out[nfi,:,:] = (1 - wei_z[:,0].reshape(si_z,1))*psi_n_out1[wei_z[:,1].astype(int),:] + wei_z[:,0].reshape(si_z,1)*psi_n_out1[wei_z[:,2].astype(int),:]
      v_w_out[nfi,:,:] = (1 - wei_z[:,0].reshape(si_z,1))*psi_w_out1[wei_z[:,1].astype(int),:] + wei_z[:,0].reshape(si_z,1)*psi_w_out1[wei_z[:,2].astype(int),:]
      v_e_out[nfi,:,:] = (1 - wei_z[:,0].reshape(si_z,1))*psi_e_out1[wei_z[:,1].astype(int),:] + wei_z[:,0].reshape(si_z,1)*psi_e_out1[wei_z[:,2].astype(int),:]
    elif myvar == 'theta':              
      t_s_out[nfi,:,:] = (1 - wei_z[:,0].reshape(si_z,1))*psi_s_out1[wei_z[:,1].astype(int),:] + wei_z[:,0].reshape(si_z,1)*psi_s_out1[wei_z[:,2].astype(int),:]
      t_n_out[nfi,:,:] = (1 - wei_z[:,0].reshape(si_z,1))*psi_n_out1[wei_z[:,1].astype(int),:] + wei_z[:,0].reshape(si_z,1)*psi_n_out1[wei_z[:,2].astype(int),:]
      t_w_out[nfi,:,:] = (1 - wei_z[:,0].reshape(si_z,1))*psi_w_out1[wei_z[:,1].astype(int),:] + wei_z[:,0].reshape(si_z,1)*psi_w_out1[wei_z[:,2].astype(int),:]
      t_e_out[nfi,:,:] = (1 - wei_z[:,0].reshape(si_z,1))*psi_e_out1[wei_z[:,1].astype(int),:] + wei_z[:,0].reshape(si_z,1)*psi_e_out1[wei_z[:,2].astype(int),:]
    elif myvar == 'salt':               
      s_s_out[nfi,:,:] = (1 - wei_z[:,0].reshape(si_z,1))*psi_s_out1[wei_z[:,1].astype(int),:] + wei_z[:,0].reshape(si_z,1)*psi_s_out1[wei_z[:,2].astype(int),:]
      s_n_out[nfi,:,:] = (1 - wei_z[:,0].reshape(si_z,1))*psi_n_out1[wei_z[:,1].astype(int),:] + wei_z[:,0].reshape(si_z,1)*psi_n_out1[wei_z[:,2].astype(int),:]
      s_w_out[nfi,:,:] = (1 - wei_z[:,0].reshape(si_z,1))*psi_w_out1[wei_z[:,1].astype(int),:] + wei_z[:,0].reshape(si_z,1)*psi_w_out1[wei_z[:,2].astype(int),:]
      s_e_out[nfi,:,:] = (1 - wei_z[:,0].reshape(si_z,1))*psi_e_out1[wei_z[:,1].astype(int),:] + wei_z[:,0].reshape(si_z,1)*psi_e_out1[wei_z[:,2].astype(int),:]


u_s_out.astype(binprec).tofile(outputdir1 + 'u_s.box')
u_n_out.astype(binprec).tofile(outputdir1 + 'u_n.box')
u_w_out.astype(binprec).tofile(outputdir1 + 'u_w.box')
u_e_out.astype(binprec).tofile(outputdir1 + 'u_e.box')
                               
v_s_out.astype(binprec).tofile(outputdir1 + 'v_s.box')
v_n_out.astype(binprec).tofile(outputdir1 + 'v_n.box')
v_w_out.astype(binprec).tofile(outputdir1 + 'v_w.box')
v_e_out.astype(binprec).tofile(outputdir1 + 'v_e.box')
                               
t_s_out.astype(binprec).tofile(outputdir1 + 't_s.box')
t_n_out.astype(binprec).tofile(outputdir1 + 't_n.box')
t_w_out.astype(binprec).tofile(outputdir1 + 't_w.box')
t_e_out.astype(binprec).tofile(outputdir1 + 't_e.box')
                               
s_s_out.astype(binprec).tofile(outputdir1 + 's_s.box')
s_n_out.astype(binprec).tofile(outputdir1 + 's_n.box')
s_w_out.astype(binprec).tofile(outputdir1 + 's_w.box')
s_e_out.astype(binprec).tofile(outputdir1 + 's_e.box')

print("{0:d} points in OBC".format(nb_files))

if flag_plot:
  plt.subplot(221)
  plt.imshow(np.mean(u_s_out,0),interpolation='nearest'); plt.colorbar()
  plt.subplot(222)
  plt.imshow(np.mean(v_s_out,0),interpolation='nearest'); plt.colorbar()
  plt.subplot(223)
  plt.imshow(np.mean(t_s_out,0),interpolation='nearest'); plt.colorbar()
  plt.subplot(224)
  plt.imshow(np.mean(s_s_out,0),interpolation='nearest'); plt.colorbar()
  plt.savefig(outputdir2 + 'bc_s.png')
  plt.clf()

  plt.subplot(221)
  plt.imshow(np.mean(u_n_out,0),interpolation='nearest'); plt.colorbar()
  plt.subplot(222)
  plt.imshow(np.mean(v_n_out,0),interpolation='nearest'); plt.colorbar()
  plt.subplot(223)
  plt.imshow(np.mean(t_n_out,0),interpolation='nearest'); plt.colorbar()
  plt.subplot(224)
  plt.imshow(np.mean(s_n_out,0),interpolation='nearest'); plt.colorbar()
  plt.savefig(outputdir2 + 'bc_n.png')
  plt.clf()

  plt.subplot(221)
  plt.imshow(np.mean(u_w_out,0),interpolation='nearest'); plt.colorbar()
  plt.subplot(222)
  plt.imshow(np.mean(v_w_out,0),interpolation='nearest'); plt.colorbar()
  plt.subplot(223)
  plt.imshow(np.mean(t_w_out,0),interpolation='nearest'); plt.colorbar()
  plt.subplot(224)
  plt.imshow(np.mean(s_w_out,0),interpolation='nearest'); plt.colorbar()
  plt.savefig(outputdir2 + 'bc_w.png')
  plt.clf()

  plt.subplot(221)
  plt.imshow(np.mean(u_e_out,0),interpolation='nearest'); plt.colorbar()
  plt.subplot(222)
  plt.imshow(np.mean(v_e_out,0),interpolation='nearest'); plt.colorbar()
  plt.subplot(223)
  plt.imshow(np.mean(t_e_out,0),interpolation='nearest'); plt.colorbar()
  plt.subplot(224)
  plt.imshow(np.mean(s_e_out,0),interpolation='nearest'); plt.colorbar()
  plt.savefig(outputdir2 + 'bc_e.png')
  plt.clf()


#=================== Atmospheric forcing =====================

if flag_atmos:
  
  print('Interpolate atmospheric data')
  
  vars2 = ['u10', 'v10','t2', 'd2', 'ssrd', 'strd']
  nvars = len(vars2)
  
  allfiles1 = sorted(glob.glob(dir0_era + vars2[0] + '*'));
  nb_files = len(allfiles1);
  
  
  f1 = netcdf.netcdf_file(allfiles1[0],'r')
  lat_i = f1.variables['lat'][:].copy().squeeze()
  lon_i = f1.variables['lon'][:].copy().squeeze()
  f1.close()
  
  lat_i = lat_i[::-1]
  lon_i = np.where(lon_i > 180, lon_i - 360,lon_i)

  
  lon_i2,lat_i2 = np.meshgrid(lon_i,lat_i) 
  lo = np.array([ np.min(lat_i), np.min(lon_i)])  # lowest lat, lowest lon
  hi = np.array([ np.max(lat_i), np.max(lon_i) ])   # highest lat, highest lon
  
  si_t_all = 0
  
  for nv in range(0,nvars):
  
    allfiles1 = sorted(glob.glob(dir0_era + vars2[nv] + '*'));
    nb_files = len(allfiles1);
  
  
    for nfi in range(0,nb_files):
      
      file1 = allfiles1[nfi]
      f = netcdf.netcdf_file(file1,'r')
    
      psi = f.variables[vars2[nv]][:,:,:].copy().squeeze()
      f.close()
  
      psi = psi[:,::-1,:]
  
      si_t,naux1,naux2 = psi.shape
      psi_out = np.zeros((si_t,si_y,si_x))
  
      if nv == 0:
        si_t_all = si_t_all + si_t
  
      for nt in range(0,si_t):
  
        interfunc = intergrid.Intergrid( 1.0*psi[nt,:,:], lo=lo, hi=hi , verbose=0, order=3)
        psi_out[nt,:,:] = interfunc.at( xx).reshape(si_y,si_x)
  
    
      #modif unit
      if vars2[nv] == 't2':
        psi_out = psi_out - 273.16
      if vars2[nv] == 'd2':
        psi_out = 6.112*np.exp(17.67*(psi_out - 273.16)/(243.5 + (psi_out - 273.16)))*0.622/1000;
      if vars2[nv] == 'ssrd':
        psi_out = np.delete(np.diff(psi_out,1,0), slice(4, None, 5),0)/6.0/60.0/60.0*0.94 #0.94:albedo (ssrd)
      if vars2[nv] == 'strd':
        psi_out = np.delete(np.diff(psi_out,1,0), slice(4, None, 5),0)/6.0/60.0/60.0
  
      app = str(nfi)
      if nfi < 10:
        app = '0' + app
  
      psi_out.astype(binprec).tofile(outputdir1 + vars2[nv] + '.box' + app)
  
    cmd = 'cat ' + outputdir1 + vars2[nv] + '.box?? > ' + outputdir1 + vars2[nv] + '.box' 
    call(cmd, shell=True)
    call('rm ' + outputdir1 + vars2[nv] + '.box??', shell=True)
  
    if flag_plot:
  
      plt.imshow(np.mean(psi_out[:,::-1,:],0),interpolation='nearest'); plt.colorbar()
      plt.savefig(outputdir2 + vars2[nv] + '.png')
      plt.clf()
  
  
  print("{0:d} points in atm forcing".format(si_t_all))


#=================== Write data files  =====================

with open("data.obcs_0", "r") as sources:
    lines = sources.readlines()
with open("data.obcs", "w") as sources:
    for line in lines:
      line2 = re.sub(r'MY_SI_Y', str(si_y), line)
      line2 = re.sub(r'MY_SI_X', str(si_x), line2)
      sources.write(line2)

with open("../code/SIZE.h_0", "r") as sources:
    lines = sources.readlines()
with open("../code/SIZE.h", "w") as sources:
    for line in lines:
      line2 = re.sub(r'MYSNX', str(int(si_x/npx)), line)
      line2 = re.sub(r'MYSNY', str(int(si_y/npy)), line2)
      line2 = re.sub(r'MYNPX', str(npx), line2)
      line2 = re.sub(r'MYNPY', str(npy), line2)
      line2 = re.sub(r'MYNR', str(si_z), line2)
      sources.write(line2)
      
nml_obcs = f90nml.read('data.obcs_0')
