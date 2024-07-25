
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
import re, fnmatch
from support_mixing import *
from jpo_figures_support import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable

#############################
########### TABLE/FIGURE 1 SM
#############################

fig, ax = plt.subplots(figsize=(7.2,3),nrows=1,ncols=1)  
plt.rcParams["figure.autolayout"] = True 
ncoldrop = 7

def plot_resolution_table(data,Re,Pr,Le,flag=None):
  t, dx, dy, z, iss, epsilon05, epsilon07, epsilon09, mkemean, tkemean, epsilonmean, epsiloncorrmean, epsilon05mean, epsilon07mean, epsilon09mean = dissipation_history(Re,Pr,Le,flag)
  dz = z[1:]-z[:-1]
  iz05 = np.argmin(np.abs(z-0.5))
  epsilonmean_ave = np.sum(epsilonmean[iz05:]*dz[iz05-1:])/np.sum(dz[iz05-1:])
  epsilonmean_ave2 = np.sum(epsilonmean[1:]*dz)/np.sum(dz)
  #
  #print(dz.max()/epsilonmean_ave)
  ratio = dz*epsilonmean[:-1]
  ratio[z[:-1]<0.5] = 0*ratio[z[:-1]<0.5]
  #ratio = ratio[z[:-1]<0.5]
  iratio = np.argmax(ratio)
  tkemean[z<0.5] = 0*tkemean[z<0.5]
  itke = np.argmax(tkemean)
  #print(ratio.max())
  #
  epsiloncorrmean_ave = np.sum(epsiloncorrmean[1:]*dz)/np.sum(dz)
  eta_k = 1/((epsilonmean_ave)**(1/4))*Re
  eta_k2 = 1/((epsilonmean[iratio])**(1/4))*Re
  eta_k3 = 1/((epsilonmean[itke])**(1/4))*Re
  eta_k4 = 1/((epsilonmean_ave2)**(1/4))*Re
  eta_bT = eta_k/Pr**(1/2)
  Sc = Pr*Le
  eta_bS = eta_k/Sc**(1/2)
  print(r'%1.2f, %1.2f, %1.2f, %1.2f'%(eta_k,eta_k2,eta_k3,eta_k4),r'  %1.2f  '%((eta_k3-eta_k)/eta_k))
  #print(r'%1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f'%(z[itke],z[iratio],dz[itke]*Re/(eta_k3/Sc**(1/2)),dz[iratio]*Re/(eta_k2/Sc**(1/2)),dz.max()*Re/(eta_k/Sc**(1/2)),eta_k3,eta_k2,eta_k))
  #print(np.sum(z*Re<1),z[np.sum(z*Re<1)]*Re,z[np.sum(z*Re<1)-1]*Re)
  #print('premiers points z ',z[:4],z[-1])
  #
  _, _, _, Nx, Ny, Nz, _, _, _ = varnumparams(Re, 1, Pr, Le, 1e-5, 32)
  #
  _, _, T, S = TSu0_profiles(Re,Pr,Le,flag)
  Sc = Pr*Le
  Sc, RaS, Rrho, St, sigma, gamma, Ri  = varphyparams(Re, 1, Pr, Le, 1e-5)
  dTdz = np.gradient(T)/np.gradient(z)
  dSdz = np.gradient(S)/np.gradient(z)
  Reb_var = epsilonmean*Sc/(RaS*(-dSdz+Rrho*dTdz))
  iz02 = np.argmin(np.abs(z-0.2))
  Reb_ave = np.sum(Reb_var[iz02+1:]*dz[iz02:])/np.sum(dz[iz02:])
  #
  deltap_u, deltap_T, deltap_S = deltaps(Re,Pr,Le,flag) 
  #
  local_row = ['%i, %i, %i'%(Re,Pr,Le), '%i, %i, %i'%(Nx,Ny,Nz), '%1.1f, %1.1f'%(eta_k,eta_bS), r'%1.1f, %1.1f, %1.1f'%(dx*Re/eta_bS,dy*Re/eta_bS,dz.max()*Re/eta_bS),r'%1.1f, %1.1f, %1.1f'%(deltap_u,deltap_T,deltap_S),r'%i'%(np.sum(z*Re*Sc**(1/2)<1)),r'%i'%(np.sum(z*Re<deltap_S)),'%1.1f'%(Re*t[-1]-Re*t[iss]), '%1.1e'%epsilonmean_ave, '%1.1e'%epsiloncorrmean_ave, r'%1.1f, %1.1f, %1.2f-%1.1f'%(dx*Re,dy*Re,z[0]*Re,dz.max()*Re), r'%1.1f, %1.1f, %1.2f-%1.1f'%((dx*Re)/eta_k,(dy*Re)/eta_k,(z[0]*Re)/eta_k,(dz.max()*Re)/eta_k),r'%1.1e'%Reb_ave,r'%i'%(np.sum(z*Re<1)),r'%i'%(np.sum(z*Re<10*deltap_S/deltap_u))]
  #
  lastkept = len(local_row)-ncoldrop
  data.append(local_row[:lastkept])
  return data

data = [] # \frac{\Delta_x}{\eta_S},\frac{\Delta_y}{\eta_S},\frac{\Delta_z}{\eta_S}
headers = [r'$(Re_{\tau},Pr,Le)$',r'$(n_x,n_y,n_z)$', r'$(\tilde{\eta}_K^+,\tilde{\eta}_{B}^+)$', r'$(\tilde{\Delta}_x,\tilde{\Delta}_y,\tilde{\Delta}_{z}^{max})\tilde{\eta}_{B}^{-1}$', r'$(\tilde{\delta}_u^+,\tilde{\delta}_T^+,\tilde{\delta}_S^+)$',r'$\mathcal{N}(Sc^{-1/2})$',r'$\mathcal{N}(\tilde{\delta}_S^+)$', r'$Re_{\tau}\tilde{\Delta}_t$', r'$\varepsilon$', r'$\varepsilon_{corr}$', r'$\Delta_x^+,\Delta_y^+,\Delta_z^+$', r'$\frac{\Delta_x}{\eta_k},\frac{\Delta_y}{\eta_k},\frac{\Delta_z}{\eta_k}$',r'$Re_b$',r'$\mathcal{N}(z^+<1)$',r'$\mathcal{N}(z^+<10\delta_S^+/\delta_u^+)$']
lastkept = len(headers)-ncoldrop
data.append(headers[:lastkept])

plot_resolution_table(data,200,1,1)
plot_resolution_table(data,200,1,10)
plot_resolution_table(data,200,1,30)
plot_resolution_table(data,200,1,100)
plot_resolution_table(data,200,5,5)
plot_resolution_table(data,200,10,10)
  
plot_resolution_table(data,400,1,1)
plot_resolution_table(data,400,1,10)
plot_resolution_table(data,400,10,10)
  
plot_resolution_table(data,800,1,1)
plot_resolution_table(data,800,1,10)

column_headers = data.pop(0)
#row_headers = [x.pop(0) for x in data]

the_table = ax.table(cellText=data,rowLoc='center',colLabels=column_headers,loc='center', cellLoc='center') #,rowLabels=row_headers
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)# Hide axes border
plt.box(on=None)
the_table.auto_set_font_size(False)
the_table.set_fontsize(9)
the_table.scale(1, 1.2)
the_table.auto_set_column_width(col=list(range(len(column_headers))))
plt.tight_layout(pad=0.2)
fig.savefig('jpoSMfig1.png', dpi=300, pad_inches=0.1)


####################
########### FIGURE 2
####################
nrows = 1; ncols = 2
fig, ax = plt.subplots(figsize=(4*1.1*ncols,2.5*1.1*nrows),nrows=nrows,ncols=ncols)

#################################
### FIG2b: History of Heat Fluxes

ax[0].axvline(x=40,linestyle=':',color='k')
ax[1].axvline(x=40,linestyle=':',color='k')

def plot_fluxes_history(fig,ax,Re,Pr,Le,flag=None):
  t, iss, qTz1_normalized, qSz1_normalized, NuT, NuS, CD, Tz1, Sz1, u0, dzTz1, dzSz1 = fluxes_history(Re,Pr,Le,flag)
  u0mean = mytimeaverage(u0[iss:],t[iss:])
  if flag=='LYRSP': ls = '--'; cl = 'b'; flag=r'low $n_y$'
  elif flag=='HRSP_Z': ls = ':'; cl = 'y'; flag=r'high $n_z$'
  elif flag=='HXRSP': ls = '--'; cl = 'r'; flag=r'high $n_x,n_z$'
  elif flag=='HRSP': ls = '--'; cl = 'g'; flag=r'high $n_y,n_z$'
  else: ls = '-'; cl = 'm'; flag='ref.'
  ax[0].plot(t*Re,u0,ls=ls,linewidth=2,color=cl,label=flag)
  ax[0].plot([t[iss]*Re,t[-1]*Re],[u0mean,u0mean],color=cl,linewidth=1)
  if flag==None: flag='ref'
  ax[1].plot(t*Re,qTz1_normalized,ls=ls,linewidth=2,color=cl,label=flag)
  ax[1].plot([t[iss]*Re,t[-1]*Re],[NuT,NuT],color=cl,linewidth=1)
  return fig, ax
  
plot_fluxes_history(fig,ax,200,1,100)
plot_fluxes_history(fig,ax,200,1,100,'HRSP')
plot_fluxes_history(fig,ax,200,1,100,'HXRSP')
plot_fluxes_history(fig,ax,200,1,100,'LYRSP')
plot_fluxes_history(fig,ax,200,1,100,'HRSP_Z')

ax[0].set_ylabel(r'$\langle \tilde{u} \rangle$',rotation=0,labelpad=15,y=0.45,fontsize='large')
ax[0].set_xlabel(r'$Re_{\tau}\tilde{t}$',fontsize='large')
ax[1].set_ylabel(r'$\frac{\tilde{Q}_T(\tilde{z}=1)}{\tilde{Q}_T^{diff}}$',rotation=0,fontsize='x-large',labelpad=20,y=0.4)
ax[1].set_xlabel(r'$Re_{\tau}\tilde{t}$',fontsize='large')
ax[0].legend(frameon=False,fontsize='small',ncol=2)

plt.tight_layout(pad=0.02)
plt.savefig('jpoSMfig2.eps')
plt.savefig('jpoSMfig2.png',dpi=300)


#############################
########### TABLE/FIGURE 3 SM
#############################

fig, ax = plt.subplots(figsize=(7,1.6),nrows=1,ncols=1)  
plt.rcParams["figure.autolayout"] = True 

def plot_SMresults_table(data,Re,Pr,Le,flag=None):
  #
  if flag=='LYRSP': nx, ny, nz = (256, 128, 192)
  elif flag=='HRSP_Z':  nx, ny, nz = (256, 256, 384)
  elif flag=='HXRSP':  nx, ny, nz = (384, 256, 256)
  elif flag=='HRSP':  nx, ny, nz = (256, 384, 256)
  else:  nx, ny, nz = (256, 256, 192); flag='reference'
  #
  t, iss, qTz1_normalized, qSz1_normalized, NuT, NuS, CD, Tz1, Sz1, u0, dzTz1, dzSz1 = fluxes_history(Re,Pr,Le,flag)
  u0mean = np.mean(u0[iss:])
  u0std = np.std(u0[iss:])
  NuTmean = np.mean(qTz1_normalized[iss:])
  NuTstd = np.std(qTz1_normalized[iss:])
  NuSmean = np.mean(qSz1_normalized[iss:])
  NuSstd = np.std(qSz1_normalized[iss:])
  Sz1mean = np.mean(Sz1[iss:])
  Sz1std = np.std(Sz1[iss:])
  QTz1mean = -1/Pr*np.mean(dzTz1[iss:])
  QTz1std = 1/Pr*np.std(dzTz1[iss:])
  #
  z, Ku0, KT, KS = turbulent_diffusivities_profiles(Re,Pr,Le,flag)
  Sc = Pr*Le
  Ku0eff = Ku0 + 1
  KTeff = KT + 1/Pr
  KSeff = KS + 1/Sc
  Lee = mybulkaverage(z,KTeff)/mybulkaverage(z,KSeff)
  print(nx,ny,nz,'$Le^{e}$ is ',Lee)
  #
  local_row = ['%i, %i, %i'%(nx,ny,nz), '%i'%u0mean+r'($\pm$%i)'%u0std, '%1.1f'%NuTmean+r'($\pm$%1.1f)'%NuTstd, '%1.1f'%NuSmean+r'($\pm$%1.1f)'%NuSstd, '%1.1e'%Sz1mean+r'($\pm$%1.1e)'%Sz1std, '%1.1f'%QTz1mean+r'($\pm$%1.1f)'%QTz1std, '%1.2f'%Lee,'%1.1f'%(Re*t[-1]-Re*t[iss])]
  #
  data.append(local_row[:lastkept])
  return data

data = [] 
headers = [r'$(n_x,n_y,n_z)$', r'$\langle \overline{\tilde{u}} \rangle$ ($\pm \sigma$)', r'$Nu_T$ ($\pm \sigma$)', r'$Nu_S$ ($\pm \sigma$)', r'$\langle \overline{\tilde{S}}(\tilde{z}=1) \rangle_{\perp}$ ($\pm \sigma$)', r'$\langle \overline{\tilde{Q}_T}(\tilde{z}=1) \rangle$ ($\pm \sigma$)', r'$Le^e$', r'$Re_{\tau}\tilde{\Delta}_t$']
data.append(headers)

plot_SMresults_table(data,200,1,100,'LYRSP')
plot_SMresults_table(data,200,1,100)
plot_SMresults_table(data,200,1,100,'HXRSP')
plot_SMresults_table(data,200,1,100,'HRSP')
plot_SMresults_table(data,200,1,100,'HRSP_Z')

column_headers = data.pop(0)
#row_headers = [x.pop(0) for x in data]

the_table = ax.table(cellText=data,rowLoc='center',colLabels=column_headers,loc='center', cellLoc='center') #,rowLabels=row_headers
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)# Hide axes border
plt.box(on=None)
the_table.auto_set_font_size(False)
the_table.set_fontsize(9)
the_table.scale(1, 1.2)
the_table.auto_set_column_width(col=list(range(len(column_headers))))
plt.tight_layout(pad=0.1)
fig.savefig('jpoSMfig3.png', dpi=300, pad_inches=0.1)




