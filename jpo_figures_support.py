
import numpy as np
from matplotlib import cm
import re, fnmatch
from support_mixing import *
from mpl_toolkits.axes_grid1 import make_axes_locatable

### GENERIC PLOTTING FUNCTIONS
def colorbar(mappable,orientation='vertical',size="5%"):
  ax = mappable.axes
  fig = ax.figure
  divider = make_axes_locatable(ax)
  if orientation=='vertical':
    cax = divider.append_axes(position="top", size=size, pad=0.05)
  else:
    cax = divider.append_axes("top", size=size, pad=0.1)
    #cax.xaxis.set_ticks_position('top')
  return fig.colorbar(mappable, orientation=orientation, cax=cax)

def colorbar(mappable,orientation='vertical',size="5%"):
  ax = mappable.axes
  fig = ax.figure
  divider = make_axes_locatable(ax)
  if orientation=='vertical':
    cax = divider.append_axes("right", size=size, pad=0.05)
  else:
    cax = divider.append_axes("top", size=size, pad=0.1)
    #cax.xaxis.set_ticks_position('top')
  return fig.colorbar(mappable, orientation=orientation, cax=cax)
  
  
def createmap(n,cmap):
  mapvec = []
  for j in range(0,n):
    idd = np.rint(j*256/(n-1))
    idd = idd.astype(int)
    mapvec.append(cmap(idd))
  return mapvec
  
### COLOR LIST
Re_list = np.array([200,400,800])
Diff_list = np.array([[1,1],[1,10],[10,10],[1,30],[1,100],[5,5]])
colormap_list = [createmap(4,cm.Greys_r), createmap(4,cm.Reds_r), createmap(4,cm.Greens_r), createmap(4,cm.Blues_r), createmap(4,cm.Purples_r), createmap(5,cm.PiYG)[1:]] #createmap(4,cm.cool_r)
color_list = [colormap_list[i] for i in range(len(Diff_list))]

### EXPERIMENT COLOR LINE 
def expcolor(Re,Pr,Le):
  iRe = np.argmin(np.abs(Re-Re_list))
  iDiff = np.argmin(np.mean((Diff_list-[Pr,Le])**2,axis=1))
  return color_list[iDiff][iRe]

### STATISTICAL STEADY STATE INDEX AND TIME
def sssindex(t,Re,flag=None):
  if Re*t[-1]>20:
    iss = np.argmin(np.abs(Re*t-15))
  else:
    iss = np.argmin(np.abs(Re*t-(Re*t[-1]-5)))
  tss = t[iss]
  #if flag=='HRSP_Z': 
  #  print('one HRSP detected')
  #  iss = 0; tss = t[iss]
  return iss, tss

### TIME AVERAGE WITH VARIABLE TIME STEP
def mytimeaverage(var,t):
  dt = np.gradient(t)
  return np.sum(var*dt)/np.sum(dt)

### FIND THE DATA FILE
def mydatafile(Re,Pr,Le,label='scalars',flag=None):
  ##
  if label=='scalars': 
    if flag=='wide':
      myfile = 'data-wide_scalars_%i_1_%i_%i_1.npz'%(Re,Pr,Le)  
    elif flag=='HRSP_Z':
      myfile = 'HRSP_Z_data_scalars_%i_1_%i_%i_1.npz'%(Re,Pr,Le) 
    elif flag=='HRSP':
      myfile = 'HRSP_data_scalars_%i_1_%i_%i_1.npz'%(Re,Pr,Le) 
    elif flag=='HXRSP':
      myfile = 'HXRSP_data_scalars_%i_1_%i_%i_1.npz'%(Re,Pr,Le) 
    elif flag=='LRSP':
      myfile = 'LRSP_data_scalars_%i_1_%i_%i_1.npz'%(Re,Pr,Le) 
    elif flag=='LYRSP':
      myfile = 'LYRSP_data_scalars_%i_1_%i_%i_1.npz'%(Re,Pr,Le) 
    elif Le==100 or Le==30:
      myfile = 'SP_data_scalars_%i_1_%i_%i_1.npz'%(Re,Pr,Le)
    elif flag==None:
      myfile = 'data_scalars_%i_1_%i_%i_1.npz'%(Re,Pr,Le)
    else:
      print('CAUTION: No available data with this flag: %s'%flag)
  ##
  elif label=='profiles':
    if flag=='wide':
      myfile = 'data-wide_profiles_%i_1_%i_%i_1.npz'%(Re,Pr,Le)
    elif flag=='HRSP_Z':
      myfile = 'HRSP_Z_data_profiles_%i_1_%i_%i_1.npz'%(Re,Pr,Le) 
    elif flag=='HRSP':
      myfile = 'HRSP_data_profiles_%i_1_%i_%i_1.npz'%(Re,Pr,Le) 
    elif flag=='HXRSP':
      myfile = 'HXRSP_data_profiles_%i_1_%i_%i_1.npz'%(Re,Pr,Le) 
    elif flag=='LRSP':
      myfile = 'LRSP_data_profiles_%i_1_%i_%i_1.npz'%(Re,Pr,Le) 
    elif flag=='LYRSP':
      myfile = 'LYRSP_data_profiles_%i_1_%i_%i_1.npz'%(Re,Pr,Le)
    elif Le==100 or Le==30:
      myfile = 'SP_data_profiles_%i_1_%i_%i_1.npz'%(Re,Pr,Le)
    elif flag==None:
      myfile = 'data_profiles_%i_1_%i_%i_1.npz'%(Re,Pr,Le)
    else:
      print('CAUTION: No available data with this flag: %s'%flag)
  ##
  elif label=='dissipationrate':
    if flag=='HRSP_Z':
      myfile = 'HRSP_Z_data_dissipationrate_%i_1_%i_%i_1.npz'%(Re,Pr,Le) 
    else:
      myfile = 'data_dissipationrate_%i_1_%i_%i_1.npz'%(Re,Pr,Le)
  ##
  elif label=='spectral':
    if flag=='HRSP_Z':
      myfile = 'HRSP_Z_data_spectral_%i_1_%i_%i_1.npz'%(Re,Pr,Le) 
    elif flag=='HRSP':
      myfile = 'HRSP_data_spectral_%i_1_%i_%i_1.npz'%(Re,Pr,Le) 
    elif flag=='HXRSP':
      myfile = 'HXRSP_data_spectral_%i_1_%i_%i_1.npz'%(Re,Pr,Le) 
    elif flag=='LRSP':
      myfile = 'LRSP_data_spectral_%i_1_%i_%i_1.npz'%(Re,Pr,Le) 
    elif flag=='LYRSP':
      myfile = 'LYRSP_data_spectral_%i_1_%i_%i_1.npz'%(Re,Pr,Le) 
    elif flag=='wide':
      myfile = 'data-wide_spectral_%i_1_%i_%i_1.npz'%(Re,Pr,Le)
    elif flag==None:
      myfile = 'data_spectral_%i_1_%i_%i_1.npz'%(Re,Pr,Le) 
    else:
      print('CAUTION: No available data with this flag: %s'%flag)
  return myfile

### FLUXES HISTORY
def fluxes_history(Re,Pr,Le,flag=None):
  # NB: qT and qS in the archives are not the full fluxes because they lack the relaxation terms
  #if flag=='HRSP':
  #  scalars = np.load('HRSP_data_scalars_%i_1_%i_%i_1.npz'%(Re,Pr,Le),'r')
  #else:
  scalars = np.load(mydatafile(Re,Pr,Le,'scalars',flag),'r')
  #print(flag,mydatafile(Re,Pr,Le,'scalars',flag))
  Sc = Pr*Le
  t = scalars['t'][1:] 
  dzTz1 = scalars['dzTz1'][1:]
  dzSz1 = scalars['dzSz1'][1:]
  Tz1 = scalars['Tz1'][1:]	
  Sz1 = scalars['Sz1'][1:]
  u0 = scalars['u0'][1:] # volume averaged u0
  #print(u0.shape)
  #
  iss, tss = sssindex(t,Re,flag)
  #
  qTz1_normalized = -(1/Pr)*dzTz1/((1/Pr)*(1-Tz1[iss:].mean()))
  qSz1_normalized = -(1/Sc)*dzSz1/((1/Sc)*(0-Sz1[iss:].mean()))
  #
  NuT = qTz1_normalized[iss:].mean()
  NuS = qSz1_normalized[iss:].mean()
  CD = 2*(Re/(u0[iss:]).mean())**2
  #
  if (NuT-mytimeaverage(qTz1_normalized[iss:],t[iss:]))/NuT > 1e-2:
    print('time averages must account for variable time steps') 
  #
  return t, iss, qTz1_normalized, qSz1_normalized, NuT, NuS, CD, Tz1, Sz1, u0, dzTz1, dzSz1
  
### FLUXES PROFILES  
def fluxes_profiles(Re,Pr,Le,flag=None):
  profiles = np.load(mydatafile(Re,Pr,Le,'profiles',flag),'r')
  #print(flag,mydatafile(Re,Pr,Le,'profiles',flag))
  Sc = Pr*Le
  z = profiles['z']
  Tu2 = profiles['Tu2']
  Su2 = profiles['Su2']
  dzT = profiles['dzT']
  dzS = profiles['dzS']
  qT = Tu2-(1/Pr)*dzT
  qS = Su2-(1/Sc)*dzS
  Tz1 = profiles['T'][-1]
  Sz1 = profiles['S'][-1]
  qT_normalized = qT/((1/Pr)*(1-Tz1))
  qS_normalized = qS/((1/Sc)*(0-Sz1))
  #
  return z, qT_normalized, qS_normalized
  
### TURBULENT DIFFUSIVITIES PROFILES  
def turbulent_diffusivities_profiles(Re,Pr,Le,flag=None):
  profiles = np.load(mydatafile(Re,Pr,Le,'profiles',flag),'r')
  Sc = Pr*Le
  z = profiles['z']
  Ku0 = profiles['Ku0']
  KT = profiles['KT']
  KS = profiles['KS']
  #
  return z, Ku0, KT, KS

### BL WIDTHS
def deltaps(Re,Pr,Le,flag=None):
  profiles = np.load(mydatafile(Re,Pr,Le,'profiles',flag),'r')
  Sc = Pr*Le
  z = profiles['z']
  u0u2 = profiles['u0u2'][z>0.6]
  Tu2 = profiles['Tu2'][z>0.6]
  Su2 = profiles['Su2'][z>0.6]
  dzu0 = profiles['dzu0'][z>0.6]
  dzT = profiles['dzT'][z>0.6]
  dzS = profiles['dzS'][z>0.6]
  z = z[z>0.6]
  iu0 = np.argmin(np.abs(u0u2+dzu0))
  iT = np.argmin(np.abs(Tu2+dzT/Pr))
  iS = np.argmin(np.abs(Su2+dzS/Sc))
  deltap_u = (1-z[iu0])*Re
  deltap_T = (1-z[iT])*Re
  deltap_S = (1-z[iS])*Re
  #
  return deltap_u, deltap_T, deltap_S
  
### DISSIPATION VARIABLES
def dissipation_history(Re,Pr,Le,flag=None):
  data = np.load(mydatafile(Re,Pr,Le,'dissipationrate',flag),'r')
  Sc = Pr*Le
  #t=t, x=x, y=y, z=z, epsilon05=epsilon05, epsilon07=epsilon07, epsilon09=epsilon09
  x = data['x']; dx = x[1]-x[0]
  y = data['y']; dy = y[1]-y[0]
  z = data['z']
  dz = z[1:]-z[:-1]
  t = data['t'][1:] 
  epsilon = data['epsilon'][1:]
  epsiloncorr = data['epsiloncorr'][1:]
  mke = data['mke'][1:]
  tke = data['tke'][1:]
  #print(z.shape,t.shape,epsilon.shape,data['epsilon05'][1:].shape)
  epsilon05 = np.mean(data['epsilon05'][1:],axis=(1,2))
  epsilon07 = np.mean(data['epsilon07'][1:],axis=(1,2))
  epsilon09 = np.mean(data['epsilon09'][1:],axis=(1,2))
  iss, tss = sssindex(t,Re)
  dt = np.gradient(t[iss:])
  epsilonmean = np.sum(epsilon[iss:]*(dt.reshape((len(t[iss:]),1))),axis=0)/np.sum(dt)
  epsiloncorrmean = np.sum(epsiloncorr[iss:]*(dt.reshape((len(t[iss:]),1))),axis=0)/np.sum(dt)
  mkemean = np.sum(mke[iss:]*(dt.reshape((len(t[iss:]),1))),axis=0)/np.sum(dt)
  tkemean = np.sum(tke[iss:]*(dt.reshape((len(t[iss:]),1))),axis=0)/np.sum(dt)
  #epsilonmean = 0
  epsilon05mean = mytimeaverage(epsilon05[iss:],t[iss:])
  epsilon07mean = mytimeaverage(epsilon07[iss:],t[iss:])
  epsilon09mean = mytimeaverage(epsilon09[iss:],t[iss:])
  diff = (np.mean(epsilon05[iss:])-epsilon05mean)/epsilon05mean
  if diff > 1e-2:
    print(diff,'time averages must account for variable time steps') 
  return t, dx, dy, z, iss, epsilon05, epsilon07, epsilon09, mkemean, tkemean, epsilonmean, epsiloncorrmean, epsilon05mean, epsilon07mean, epsilon09mean

### T & S & Ux PROFILES  
def TSu0_profiles(Re,Pr,Le,flag=None):
  profiles = np.load(mydatafile(Re,Pr,Le,'profiles',flag),'r')
  u0 = profiles['u0']
  T = profiles['T']
  S = profiles['S']
  z = profiles['z']
  return z, u0, T, S
  
### T & S & Ux PROFILES FROM TD SIMULATIONS
def TSu0_profiles_from_TD(Re,Pr):
  profiles = np.load(mydatafile(Re,Pr,1000,'profiles'),'r') #Le=1000 denotes TD results
  z = profiles['z']
  u0_from_TD = profiles['u0']
  TD = profiles['TD']
  TD100000 = 1
  gamma = lambda1*S0/(TD100000/1e5)
  St = Li/(cp*(TD100000/1e5))
  T_from_TD = (St*TD-gamma)/(St-gamma) # from TD definition and Gade's mixing line
  S_from_TD = (TD-1)/(St-gamma)
  return z, u0_from_TD, TD, T_from_TD, S_from_TD

### BULK AVERAGE
def mybulkaverage(z,var):
  dz = np.gradient(z)
  i03 = np.argmin(np.abs(z-0.3))
  i07 = np.argmin(np.abs(z-0.7))
  return np.sum(var[i03:i07]*dz[i03:i07])/np.sum(dz[i03:i07])

### MEAN FLUXES AND STD FROM TD SIMULATIONS
def fluxes_from_TD(Re,Pr,Le,flag=None):
  scalars = np.load(mydatafile(Re,Pr,1000,'scalars',flag),'r') #Le=1000 denotes TD results
  t = scalars['t'][1:]
  dzTDz1 = scalars['dzTDz1'][1:]
  u0_from_TD = scalars['u0'][1:]
  iss, tss = sssindex(t,Re)
  u0m_from_TD = u0_from_TD[iss:].mean()
  dzTDz1m = dzTDz1[iss:].mean()  
  dzTDz1_std_high = np.std(dzTDz1[iss:],where=dzTDz1[iss:]>dzTDz1m) 
  dzTDz1_std_low = np.std(dzTDz1[iss:],where=dzTDz1[iss:]<dzTDz1m) 
  if (dzTDz1m-mytimeaverage(dzTDz1[iss:],t[iss:]))/dzTDz1m > 1e-2:
    print('time averages must account for variable time steps') 
  TD100000 = 1
  NuTD = -1/Pr*dzTDz1[iss:].mean()/(1/Pr)
  gamma = lambda1*S0/(TD100000/1e5)
  St = Li/(cp*(TD100000/1e5))
  dzTz1m_from_dzTDz1m = (St*dzTDz1m)/(St-gamma) # from TD definition and Gade's mixing line
  dzSz1m_from_dzTDz1m = (dzTDz1m)/(St-gamma)
  #
  profiles = np.load(mydatafile(Re,Pr,1000,'profiles'),'r')
  z = profiles['z']
  KTD = profiles['KTD']
  TD = profiles['TD']
  dzTD = profiles['dzTD']
  T = (St*TD-gamma)/(St-gamma) # from TD definition and Gade's mixing line
  S = (TD-1)/(St-gamma)
  qT_from_TD = -(KTD+1/Pr)*St*dzTD/(St-gamma)
  qS_from_TD = -(KTD+1/Pr)*dzTD/(St-gamma)
  qT_from_TD_bulk = mybulkaverage(z,qT_from_TD)   
  qS_from_TD_bulk = mybulkaverage(z,qS_from_TD) 
  Sc = Pr*Le
  dzTz1m_from_TD_bulk = -qT_from_TD_bulk/(1/Pr)
  dzSz1m_from_TD_bulk = -qS_from_TD_bulk/(1/Sc)
  #
  #KTD_bulk = mybulkaverage(z,KTD)
  #dzT_from_TD_bulk = (T[i07]-T[i03])/(z[i07]-z[i03])
  #dzS_from_TD_bulk = (S[i07]-S[i03])/(z[i07]-z[i03])
  #
  return t, iss, dzTz1m_from_dzTDz1m, dzSz1m_from_dzTDz1m, dzTz1m_from_TD_bulk, dzSz1m_from_TD_bulk, dzTDz1_std_low, dzTDz1_std_high, u0m_from_TD, NuTD
  
### MEAN FLUXES AND STD FROM TS SIMULATIONS
def fluxes_from_TS(Re,Pr,Le,flag=None):
  scalars = np.load(mydatafile(Re,Pr,Le,'scalars',flag),'r')
  t = scalars['t'][1:]
  dzTz1 = scalars['dzTz1'][1:]
  dzSz1 = scalars['dzSz1'][1:]
  iss, tss = sssindex(t,Re,flag)
  dzTz1m = dzTz1[iss:].mean() 
  dzSz1m = dzSz1[iss:].mean()   
  dzTz1_std_high = dzTz1[iss:].std(where=dzTz1[iss:]>dzTz1m)
  dzTz1_std_low = dzTz1[iss:].std(where=dzTz1[iss:]<dzTz1m)
  dzSz1_std_high = dzSz1[iss:].std(where=dzSz1[iss:]>dzSz1m)
  dzSz1_std_low = dzSz1[iss:].std(where=dzSz1[iss:]<dzSz1m)
  return dzTz1m, dzSz1m, dzTz1_std_low, dzTz1_std_high, dzSz1_std_low, dzSz1_std_high
  
