
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

####################
########### FIGURE 1
####################

####################
########### FIGURE 2
####################

####################
########### FIGURE 3
####################
nrows = 1; ncols = 3
fig, ax = plt.subplots(figsize=(2.5*ncols,2.5*nrows),nrows=nrows,ncols=ncols)

#################################
### FIG3a: History of Heat Fluxes
ax[0].plot(10,10,'-',linewidth=2,color='white',label=r'$Re_{\tau}$, $Pr$, $Le$')

def plot_fluxes_history(fig,ax,Re,Pr,Le,flag=None):
  t, iss, qTz1_normalized, qSz1_normalized, NuT, NuS, CD, _, _, _, _, _ = fluxes_history(Re,Pr,Le,flag)
  ax[0].plot(t*Re,qTz1_normalized,'-',linewidth=2,color=expcolor(Re,Pr,Le),label=r'%i, %i, %i'%(Re,Pr,Le))
  ax[0].plot([t[iss]*Re,t[-1]*Re],[NuT,NuT],'-k',linewidth=1)
  return fig, ax
  
Pr = 1
Le = 10
for Re in [200,400,800]:
  plot_fluxes_history(fig,ax,Re,Pr,Le)
plot_fluxes_history(fig,ax,200,1,100)
plot_fluxes_history(fig,ax,400,10,10)

ax[0].set_yscale('log')
ax[0].set_ylabel(r'$\frac{\tilde{Q}_T(\tilde{z}=1)}{\tilde{Q}_T^{diff}}$',rotation=0,fontsize='large',labelpad=5,y=0.4)
ax[0].set_xlabel(r'$Re_{\tau}\tilde{t}$')
ax[0].legend(frameon=False,fontsize='x-small')

####################################
### FIG3b: z Profiles of Heat Fluxes

def plot_fluxes_profiles(fig,ax,Re,Pr,Le):
  z, qT_normalized, qS_normalized = fluxes_profiles(Re,Pr,Le,flag=None)
  ax[1].plot(qT_normalized,z,'-',linewidth=2,color=expcolor(Re,Pr,Le),label=r'$Re_{\tau}=$%i'%Re)
  ax[1].axvline(qT_normalized[-1],zorder=1,color=expcolor(Re,Pr,Le),linestyle=':',linewidth=1)
  return fig, ax
  
Pr = 1
Le = 10
for Re in [200,400,800]:
  plot_fluxes_profiles(fig,ax,Re,Pr,Le)
plot_fluxes_profiles(fig,ax,200,1,100)
plot_fluxes_profiles(fig,ax,400,10,10)

ax[1].axhline(0.2,zorder=1,color='k',linestyle='--',linewidth=1)
ax[1].set_xlabel(r'$\overline{\tilde{Q}_T}/\tilde{Q}_T^{diff}$')
ax[1].set_ylabel(r'$\tilde{z}$',rotation=0,labelpad=15,y=0.45)
ax[1].set_ylim((0,1))

############################################
### FIG3c: z Profiles of Horizontal Velocity

def plot_mean_flow_profiles(fig,ax,Re,Pr,Le,flag=None):
  z_from_TD, u0_from_TD, TD, T_from_TD, S_from_TD = TSu0_profiles_from_TD(Re,Pr)
  z, u0, T, S = TSu0_profiles(Re,Pr,Le,flag)
  ax[2].plot(u0/Re,z,'-',linewidth=2,color=expcolor(Re,Pr,Le),label=r'%i, %i, %i'%(Re,Pr,Le))
  ax[2].plot(u0_from_TD/Re,z_from_TD,'--',linewidth=2,color=expcolor(Re,Pr,Le))
  return fig, ax

Pr = 1
Le = 10
for Re in [200,400,800]:
  plot_mean_flow_profiles(fig,ax,Re,Pr,Le)
plot_mean_flow_profiles(fig,ax,200,1,100)
plot_mean_flow_profiles(fig,ax,400,10,10)

ax[2].set_xlim((-10,35))
ax[2].set_ylim((0,1))
ax[2].set_ylabel(r'$\tilde{z}$',rotation=0,labelpad=15,y=0.45)
ax[2].set_xlabel(r'$\langle \overline{\tilde{u}}\rangle_{\perp} /Re_{\tau}$')

plt.tight_layout(pad=0.2)
plt.savefig('jpofig3.eps')
plt.savefig('jpofig3.png',dpi=300)

#############################################
########### FIGURE 4: Turbulent Diffusivities
#############################################
nrows = 1; ncols = 3
fig, ax = plt.subplots(figsize=(2.5*ncols,2.5*nrows),nrows=nrows,ncols=ncols)
ax[1].plot(2,0.5,'-',linewidth=2,color='white',label=r'$Re_{\tau}$, $Pr$, $Le$')
axins = inset_axes(ax[2], width="35%", height="35%", loc='right', borderpad=1.5)

def plot_turbulent_diffusivities_profiles(fig,ax,axins,Re,Pr,Le,flag=None):
  z, Ku0, KT, KS = turbulent_diffusivities_profiles(Re,Pr,Le,flag)
  Sc = Pr*Le
  Ku0eff = Ku0 + 1
  KTeff = KT + 1/Pr
  KSeff = KS + 1/Sc
  Lee = mybulkaverage(z,KTeff)/mybulkaverage(z,KSeff)
  print(Re,Pr,Le,'$Le^{e}$ is ',Lee)
  i02 = np.argmin(np.abs(z-0.2))
  ax[0].plot(Ku0,z,'-',linewidth=2,color=expcolor(Re,Pr,Le))
  ax[0].plot(KT[i02:],z[i02:],'-',linewidth=0.5,color=expcolor(Re,Pr,Le))
  ax[1].plot(Ku0eff[i02:]/KTeff[i02:],z[i02:],'-',linewidth=2,color=expcolor(Re,Pr,Le),label=r'%i, %i, %i'%(Re,Pr,Le))
  ax[2].plot(KTeff[i02:]/KSeff[i02:],z[i02:],'-',linewidth=2,color=expcolor(Re,Pr,Le))
  axins.plot(KTeff[i02:]/KSeff[i02:],z[i02:],'-',linewidth=2,color=expcolor(Re,Pr,Le))
  return fig, ax, axins
  
Pr = 1
Le = 10
for Re in [200,400,800]:
  plot_turbulent_diffusivities_profiles(fig,ax, axins,Re,Pr,Le)
  
plot_turbulent_diffusivities_profiles(fig,ax, axins,200,1,100)
plot_turbulent_diffusivities_profiles(fig,ax, axins,400,10,10)

[ax[i].axvline(1,color='k',linestyle='--',linewidth=1) for i in [1,2]]
[ax[i].axhline(0.2,zorder=1,color='k',linestyle='--',linewidth=1) for i in range(3)]
[ax[i].set_ylim((0,1)) for i in range(3)]
ax[0].set_xlim((0,35))
ax[1].set_xlim((0,5))
ax[2].set_xlim((0,5))
#[ax[i].set_yticklabels([]) for i in [1,2]]
ax[1].legend(frameon=False,fontsize='x-small')
ax[0].set_xlabel(r'$\tilde{\nu}^t\; (thick),\; \tilde{\kappa}_T^t\; (thin)$')
ax[1].set_xlabel(r'$Pr^{e}=\tilde{\nu}^e/\tilde{\kappa}_T^e$')
ax[2].set_xlabel(r'$Le^{e}=\tilde{\kappa}_T^e/\tilde{\kappa}_S^e$')
[ax[i].set_ylabel(r'$\tilde{z}$',rotation=0,labelpad=15,y=0.45) for i in range(3)]

axins.set_xscale('log')
axins.set_ylim((0.96,1.0))
#axins.set_xlabel(r'$Le^{e}$',labelpad=1,fontsize='small')
#axins.set_ylabel(r'$\tilde{z}$',rotation=0,labelpad=5,fontsize='small',y=0.45)
axins.set_yticks([0.96,0.98,1.0])
for label in (axins.get_xticklabels() + axins.get_yticklabels()):
    label.set_fontsize('x-small')
pos = axins.get_position()
print(pos)
axins.set_position([pos.x0, pos.y0-50, pos.width, pos.height])

plt.tight_layout(pad=0.2,w_pad=1)
plt.savefig('jpofig4.eps')
plt.savefig('jpofig4.png',dpi=300)

###############################
########### FIGURE SUPP 1: Le^e
###############################

figsupp, axsupp = plt.subplots(figsize=(2.5*1,2.5*1),nrows=1,ncols=1)

def plot_Lee(figsupp,axsupp,Re,Pr,Le,flag=None):
  z, Ku0, KT, KS = turbulent_diffusivities_profiles(Re,Pr,Le,flag)
  Sc = Pr*Le
  Ku0eff = Ku0 + 1
  KTeff = KT + 1/Pr
  KSeff = KS + 1/Sc
  Lee = mybulkaverage(z,KTeff)/mybulkaverage(z,KSeff)
  axsupp.plot(Le, Lee, linestyle=None, marker='o', markersize=7, color=expcolor(Re,Pr,Le), markeredgewidth=0, zorder=1)
  return figsupp,axsupp
  
for Re in [200,400,800]:
  plot_Lee(figsupp,axsupp,Re,1,1)
  plot_Lee(figsupp,axsupp,Re,1,10)
plot_Lee(figsupp,axsupp,200,1,30)
plot_Lee(figsupp,axsupp,200,1,100)
plot_Lee(figsupp,axsupp,200,5,5)
plot_Lee(figsupp,axsupp,200,10,10)
plot_Lee(figsupp,axsupp,400,10,10)

axsupp.legend(frameon=False,fontsize='x-small')
axsupp.set_xscale(r'log')
axsupp.set_xlabel(r'$Le$')
axsupp.set_ylabel(r'$Le^e$',rotation=0,labelpad=15,y=0.45)
figsupp.tight_layout(pad=0.2,w_pad=1)
figsupp.savefig('jposuppfig1.eps')
figsupp.savefig('jposuppfig1.png',dpi=300)

#################################
########### FIGURE 5: T-S Diagram
#################################
nrows = 1; ncols = 1
fig, ax = plt.subplots(figsize=(3.5*ncols,2.5*nrows),nrows=nrows,ncols=ncols)
ax.plot(0,0,'-',linewidth=2,color='white',label=r'$Re_{\tau}$, $Pr$, $Le$')

def plot_TS_profiles(fig,ax,Re,Pr,Le,flag=None):
  z, _, T, S = TSu0_profiles(Re,Pr,Le,flag)
  ax.plot(S,T,'-',linewidth=2,color=expcolor(Re,Pr,Le),label=r'%i, %i, %i'%(Re,Pr,Le))
  return fig, ax

Pr = 1
Le = 10
for Re in [200,400,800]:
  plot_TS_profiles(fig,ax,Re,Pr,Le)
plot_TS_profiles(fig,ax,200,1,100)
plot_TS_profiles(fig,ax,400,10,10)

TD0 = 1e-5
_, _, _, St, sigma, gamma, _ = varphyparams(200, 1, 1, 100, TD0)
Svec = np.linspace(-1.3e-6,0,3)
Tfus = gamma*Svec
TDvec = np.linspace(0,1,100)
Tbulk = ( TDvec-gamma/St ) / ( 1-gamma/St ) # Gade mixing line in the bulk
Sbulk = ( TDvec-1 ) / (St-gamma)
ax.plot(Svec,Tfus,'-k',linewidth=2,zorder=1)
ax.plot(Sbulk[Tbulk>gamma*Sbulk],1+St*Sbulk[Tbulk>gamma*Sbulk],'--k',linewidth=1,zorder=1)

ax.set_xlabel(r'$\langle \overline{\tilde{S}}\rangle_{\perp}$')
ax.set_ylabel(r'$\langle \overline{\tilde{T}}\rangle_{\perp}$',rotation=0,labelpad=15,y=0.45)
ax.legend(frameon=False,fontsize='x-small',loc='upper left')

plt.tight_layout(pad=0.2)
plt.savefig('jpofig5.eps')
plt.savefig('jpofig5.png',dpi=300)

#############################################
########### FIGURE 6: Exchange velocity ratio
#############################################

nrows = 1; ncols = 1
fig, ax = plt.subplots(figsize=(3.5*ncols,2.5*nrows),nrows=nrows,ncols=ncols)

ax.plot([1,10],[1,1*10**(-2/3)],'--k',linewidth=1.5)
ax.plot([10,100],[0.2,0.2*10**(-1/2)],':k',linewidth=1.5)
ax.text(1.5, 0.4, '-2/3', fontsize='x-small')
ax.text(15, 0.1, '-1/2', fontsize='x-small')

def transport_parameters_with_Le(fig,ax,Re,Pr,Le,flag=None):   
  if Re==200: s='o'
  elif Re==400: s='d'
  elif Re==800: s='s'
  s='o'
  t, iss, qTz1_normalized, qSz1_normalized, NuT, NuS, CD, Tz1, Sz1, u0, dzTz1, dzSz1 = fluxes_history(Re,Pr,Le,flag=None)
  ax.loglog(Le,NuS/NuT/Le,linestyle='None', marker=s,markerfacecolor="None", color=expcolor(Re,Pr,Le),markeredgewidth=2,markersize=8, label=r'%i, %i'%(Re,Pr))
  print(Re,Pr,Le,np.std(-1/Pr*dzTz1[iss:]))
  return fig, ax

ax.plot(1e1,1e0,'-',linewidth=2,color='white',label=r'$Re_{\tau}$, $Pr$')
transport_parameters_with_Le(fig,ax,200,1,1)
transport_parameters_with_Le(fig,ax,200,5,5)
transport_parameters_with_Le(fig,ax,200,1,10)
transport_parameters_with_Le(fig,ax,200,1,30)
transport_parameters_with_Le(fig,ax,200,1,100)
transport_parameters_with_Le(fig,ax,400,1,10)
transport_parameters_with_Le(fig,ax,800,1,10)
transport_parameters_with_Le(fig,ax,400,10,10)

ax.set_ylim((0.02,1.2))
ax.legend(frameon=False,fontsize='x-small',loc='lower left',ncol=1,handletextpad=0.1,markerscale=0.7)
ax.set_xlabel(r'$Le$')
ax.set_ylabel(r'$\frac{Nu_S}{Nu_TLe}$', rotation=0, labelpad=5, y=0.5, fontsize='large')

fig.tight_layout(pad=0.02,w_pad=0.1)
fig.savefig('jpofig6.eps',bbox_inches='tight')
fig.savefig('jpofig6.png',dpi=300,bbox_inches='tight')

#######################################
########### FIGURE 7: Model Comparisons
#######################################
nrows = 1; ncols = 4
fig, ax = plt.subplots(figsize=(2*ncols,2.5*nrows),nrows=nrows,ncols=ncols)
ax[0].plot(0,0,'-',linewidth=2,color='white',label=r'$Re_{\tau}$, $Pr$, $Le$')

def plot_and_compare_profiles_from_TD(fig,ax,Re,Pr,Le,flag=None):
  z_from_TD, u0_from_TD, TD, T_from_TD, S_from_TD = TSu0_profiles_from_TD(Re,Pr)
  z, u0, T, S = TSu0_profiles(Re,Pr,Le,flag)
  #
  ax[0].plot(T,z,'-',linewidth=2,color=expcolor(Re,Pr,Le), label=r'%i, %i, %i'%(Re,Pr,Le))
  ax[0].plot(T_from_TD,z_from_TD,'--',linewidth=2,color=expcolor(Re,Pr,Le))
  ax[0].plot(TD,z_from_TD,':',color='k',linewidth=1.5,zorder=3)
  #
  ax[1].plot(S*1e7,z,'-',linewidth=2,color=expcolor(Re,Pr,Le))
  ax[1].plot(S_from_TD*1e7,z_from_TD,'--',linewidth=2,color=expcolor(Re,Pr,Le))
  #
  ax[2].plot(T,(1-z)*Re,'-',linewidth=2,color=expcolor(Re,Pr,Le))
  ax[2].plot(T_from_TD,(1-z_from_TD)*Re,'--',linewidth=2,color=expcolor(Re,Pr,Le))#,alpha=0.85
  ax[2].plot(TD,(1-z_from_TD)*Re,':',color='k',linewidth=1.5,zorder=3)
  #
  ax[3].plot(S*1e7,(1-z)*Re,'-',linewidth=2,color=expcolor(Re,Pr,Le))
  ax[3].plot(S_from_TD*1e7,(1-z_from_TD)*Re,'--',linewidth=2,color=expcolor(Re,Pr,Le))
  #
  return fig, ax

for Re in [200,400,800]:
  plot_and_compare_profiles_from_TD(fig,ax,Re,1,10)
plot_and_compare_profiles_from_TD(fig,ax,200,1,100)
plot_and_compare_profiles_from_TD(fig,ax,400,10,10)

ax[0].legend(frameon=False,fontsize='x-small',loc='lower left')
[ax[i].set_ylim((0,1)) for i in [0,1]]
[ax[i].set_ylim((0,10)) for i in [2,3]]
ax[0].set_ylabel(r'$\tilde{z}$',rotation=0,labelpad=10,y=0.45)
ax[2].set_ylabel(r'$(1-\tilde{z})^+$',rotation=0,labelpad=10,y=0.45)
[ax[i].set_yticklabels([]) for i in [1,3]]
[ax[i].set_xlabel(r'$\langle \overline{\tilde{T}}\rangle_{\perp}$') for i in [0,2]]
[ax[i].set_xlabel(r'$\langle \overline{\tilde{S}}\rangle_{\perp} \times 10^7$') for i in [1,3]]
ax[0].set_xlim([-0.05,1.05])
ax[1].set_xlim([-1.3,0])
ax[2].set_xlim((0.0,0.5))
ax[3].set_xlim((-3,-1e-1))
[ax[i].invert_yaxis() for i in [2,3]]

plt.tight_layout(pad=0.2,w_pad=0)

wcap = 0.075
wocap = 0.02
wspacing = [wcap, wocap, wcap+0.02, wocap]
widths = [0.195-0.005, 0.195-0.005, 0.195-0.005, 0.195-0.005]
x0s = 0.1 + np.cumsum(np.array(wspacing)) + np.cumsum(np.array(widths))-widths[0]
for i in range(4):
  pos = ax[i].get_position()
  newpos = [x0s[i], pos.y0,  widths[i], pos.height]
  ax[i].set_position(newpos)

plt.savefig('jpofig7.eps',bbox_inches='tight')
plt.savefig('jpofig7.png',dpi=300,bbox_inches='tight')

wsum = 0
for i in range(4):
  pos=ax[i].get_position()
  print(pos.width,pos)
  wsum+=pos.width
print(wsum)


#######################################
########### FIGURE 8: Fluxes comparison
#######################################

nrows = 1; ncols = 2
fig, ax = plt.subplots(figsize=(2.5*ncols,2.5*nrows),nrows=nrows,ncols=ncols)
ax[0].plot(10,10,'-',linewidth=2,color='white',label=r'$Re_{\tau}$, $Pr$, $Le$')

def plot_and_compare_fluxes_from_TD(fig,ax,Re,Pr,Le,flag=None):   
  t, iss, dzTz1m_from_dzTDz1m, dzSz1m_from_dzTDz1m, dzTz1m_from_TD_bulk, dzSz1m_from_TD_bulk, dzTDz1_std_low, dzTDz1_std_high, u0m_from_TD, NuTD = fluxes_from_TD(Re,Pr,Le,flag)
  dzTz1m, dzSz1m, dzTz1_std_low, dzTz1_std_high, dzSz1_std_low, dzSz1_std_high = fluxes_from_TS(Re,Pr,Le,flag)
  Sc = Pr*Le
  QT = -(1/Pr)*dzTz1m
  QTdiag = -(1/Pr)*dzTz1m_from_TD_bulk
  QS = -(1/Sc)*dzSz1m
  QSdiag = -(1/Sc)*dzSz1m_from_TD_bulk
  # 
  ax[0].errorbar(-(1/Sc)*dzSz1m_from_dzTDz1m*St, -(1/Pr)*dzTz1m_from_dzTDz1m, xerr=None, yerr=None, linestyle=None, marker='x', markersize=7, color=expcolor(Re,Pr,Le), markeredgewidth=2, zorder=1)
  #
  ax[0].errorbar(QSdiag*St, QTdiag, xerr=None,yerr=None,linestyle=None, marker='s', markersize=7, color=expcolor(Re,Pr,Le),markerfacecolor="None", markeredgewidth=2, zorder=1)
  #
  ax[0].errorbar(QS*St,QT,xerr=None,yerr=None,linestyle="None", marker='o',markersize=8,markerfacecolor="None",markeredgewidth=2,zorder=3,color=expcolor(Re,Pr,Le),alpha=1,label=r'%i, %i, %i'%(Re,Pr,Le))
  #
  ax[1].errorbar(Le,(QT-QTdiag),xerr=None,yerr=[[-1/Pr*dzTz1_std_low],[-1/Pr*dzTz1_std_high]],marker='o',markersize=8,markeredgewidth=2,markerfacecolor="None",color=expcolor(Re,Pr,Le),alpha=1)
  #
  return fig, ax 
  
#for Re in [200,400,800]:
#  plot_and_compare_fluxes_from_TD(fig,ax,Re,1,1)
for Re in [200,400,800]:
  plot_and_compare_fluxes_from_TD(fig,ax,Re,1,10)
plot_and_compare_fluxes_from_TD(fig,ax,200,1,30)
plot_and_compare_fluxes_from_TD(fig,ax,200,1,100)
plot_and_compare_fluxes_from_TD(fig,ax,200,5,5)
plot_and_compare_fluxes_from_TD(fig,ax,200,10,10)
plot_and_compare_fluxes_from_TD(fig,ax,400,10,10)

ax[0].set_ylim((-5,21))
ax[0].set_xlim((-2,24))
ax[0].legend(frameon=False,fontsize='x-small',loc='lower right',ncol=1,handletextpad=0.1,markerscale=0.7)
ax[0].plot([1,20],[1,20],'--k',linewidth=1.5)
ax[0].set_xlabel(r'$\overline{\tilde{Q}_S}(\tilde{z}=1)St$')
ax[0].set_ylabel(r'$\overline{\tilde{Q}_T}(\tilde{z}=1)$')#, rotation=0, labelpad=13, y=0.45)
ax[1].set_xscale('log')
ax[1].set_xlabel(r'$Le$')
ax[1].set_ylabel(r'$\overline{\tilde{Q}_T}(\tilde{z}=1)-\overline{\tilde{Q}_T^{diag}}(\tilde{z}=1)$')

fig.tight_layout(pad=0.02,w_pad=0.1)
fig.savefig('jpofig8.eps',bbox_inches='tight')
fig.savefig('jpofig8.png',dpi=300,bbox_inches='tight')



#########################################################
################## CREATE A TABLE WITH RELEVANT VARIABLES
#########################################################

def plot_results_table(data,Re,Pr,Le,flag=None):
  print(flag)
  if Le==1000:
    Leeff = np.nan
    meltrate = np.nan
    meltratediagfromTD = np.nan
    Reb_ave = np.nan
    Reb_ave2 = np.nan
    Reb_ave3 = np.nan
    #
    t, iss, dzTz1m_from_dzTDz1m, dzSz1m_from_dzTDz1m, dzTz1m_from_TD_bulk, dzSz1m_from_TD_bulk, dzTDz1_std_low, dzTDz1_std_high, u0m_from_TD, NuTD = fluxes_from_TD(Re,Pr,Le,flag)
    #NuT = -dzTz1m_from_dzTDz1m/Pr/(1/Pr)
    NuT = NuTD
    NuS = np.nan
    CD = 2*(Re/(u0m_from_TD))**2 # already volume and time averaged
    Sz1 = np.nan
    QT = NuTD/Pr
  else:
    z, Ku0, KT, KS = turbulent_diffusivities_profiles(Re,Pr,Le,flag)
    Sc = Pr*Le
    Ku0eff = Ku0 + 1
    KTeff = KT + 1/Pr
    KSeff = KS + 1/Sc
    Leeff = mybulkaverage(z,KTeff)/mybulkaverage(z,KSeff)
    #
    t, iss, qTz1_normalized, qSz1_normalized, NuT, NuS, CD, Tz1, Sz1, u0, dzTz1, dzSz1 = fluxes_history(Re,Pr,Le,flag)
    Sz1 = mytimeaverage(Sz1[iss:],t[iss:])
    Sc, RaS, Rrho, St, sigma, gamma, Ri  = varphyparams(Re, 1, Pr, Le, 1e-5)
    meltrate = -mytimeaverage(dzTz1[iss:],t[iss:])*Pr/St
    QT = -1/Pr*mytimeaverage(dzTz1[iss:],t[iss:])
    #
    if flag!='wide' and flag!='HRSP_Z' and flag!='HRSP' and flag!='HXRSP' and flag!='LRSP' and flag!='LYRSP':
      _, _, dzTz1m_from_dzTDz1m, dzSz1m_from_dzTDz1m, dzTz1m_from_TD_bulk, dzSz1m_from_TD_bulk, dzTDz1_std_low, dzTDz1_std_high, u0m_from_TD, NuTD = fluxes_from_TD(Re,Pr,Le,flag)
      #
      _, _, _, St, sigma, gamma, _ = varphyparams(Re, 1, Pr, Le, 1e-5)
      #meltrate = -dzTz1m_from_dzTDz1m*Pr/St
      meltratediagfromTD = -dzTz1m_from_TD_bulk*Pr/St
    else:
      meltrate = np.nan
      meltratediagfromTD = np.nan
    #
  #NuTdiagfromTD = -dzTz1m_from_TD_bulk/Pr # ---> not possible as we don't have access to the normalization T(1) or S(1)
  #NuSdiagfromTD = -dzSz1m_from_TD_bulk/Sc
  #
    t, dx, dy, z, iss, epsilon05, epsilon07, epsilon09, mkemean, tkemean, epsilonmean, epsiloncorrmean, epsilon05mean, epsilon07mean, epsilon09mean = dissipation_history(Re,Pr,Le,flag)
    dz = z[1:]-z[:-1]
    iz05 = np.argmin(np.abs(z-0.5))
    epsilonmean_ave = np.sum(epsilonmean[iz05:]*dz[iz05-1:])/np.sum(dz[iz05-1:])
    epsiloncorrmean_ave = np.sum(epsiloncorrmean[iz05:]*dz[iz05-1:])/np.sum(dz[iz05-1:])
    eta_k = 1/((epsilonmean_ave)**(1/4))*Re
    eta_bT = eta_k/Pr**(1/2)
    Sc = Pr*Le
    eta_bS = eta_k/Sc**(1/2)
  #
    _, _, _, Nx, Ny, Nz, _, _, _ = varnumparams(Re, 1, Pr, Le, 1e-5, 32)
  #
    _, _, T, S = TSu0_profiles(Re,Pr,Le,flag)
    Sc = Pr*Le
    Sc, RaS, Rrho, St, sigma, gamma, Ri  = varphyparams(Re, 1, Pr, Le, 1e-5)
    dTdz = np.gradient(T)/np.gradient(z)
    dSdz = np.gradient(S)/np.gradient(z)
    drhodz = RaS*(-dSdz+Rrho*dTdz)
    drhodz_ave = np.sum(drhodz[iz05:]*dz[iz05-1:])/np.sum(dz[iz05-1:])
    Reb_var = epsilonmean*Sc/drhodz
    Reb_ave = epsilonmean_ave*Sc/drhodz_ave
    Reb_ave2 = mybulkaverage(z,epsilonmean*Sc/drhodz)
    Reb_ave3 = np.sum(Reb_var[iz05:]*dz[iz05-1:])/np.sum(dz[iz05-1:])
  #
  local_row = ['Re %i Pr %i Le %i'%(Re,Pr,Le), '%s'%flag, '%1.3f'%Leeff, '%1.1e'%Reb_ave, '%1.1e'%Reb_ave2, '%1.1e'%Reb_ave3, '%1.1f'%NuT, '%1.1f'%NuS, '%1.1f'%(1000*CD), '%1.1f'%(1e6*meltrate), '%1.1f'%(QT), '%1.1f'%(Re*t[iss]), '%1.1f'%(Re*t[-1]-Re*t[iss])]
  data.append(local_row)
  return data

data = []
data.append([r'flag', r'$Le^{eff}$', r'$Re_b$', r'$Re_b$', r'$Re_b$', r'$Nu_T$', r'$Nu_S$', r'$10^3C_D$', r'$10^6\dot{m}$', r'$Q_T$', r'$Re_{\tau}t_0$', r'$Re_{\tau}\Delta t$']) #,  r'$10^6\langle\overline{S}(z=1)\rangle$'

plot_results_table(data,200,1,1) 
#plot_results_table(data,200,1,10,flag='LRSP') 
plot_results_table(data,200,1,10) 
#plot_results_table(data,200,5,5,flag='LRSP') 
plot_results_table(data,200,5,5) 
#plot_results_table(data,200,1,30,flag='LRSP') 
plot_results_table(data,200,1,30) 
#plot_results_table(data,200,1,30,flag='HRSP_Z') 

#plot_results_table(data,200,1,100,flag='LYRSP') 
plot_results_table(data,200,1,100) 
#plot_results_table(data,200,1,100,flag='HRSP_Z') 
#plot_results_table(data,200,1,100,flag='HRSP') 
#plot_results_table(data,200,1,100,flag='HXRSP') 
plot_results_table(data,200,10,10) 

plot_results_table(data,400,1,1) 
#plot_results_table(data,400,1,10,flag='LRSP')
plot_results_table(data,400,1,10) 
plot_results_table(data,400,1,10,flag='wide')
plot_results_table(data,400,10,10) 

plot_results_table(data,800,1,1) 
plot_results_table(data,800,1,10) 

plot_results_table(data,200,1,1000) 
plot_results_table(data,200,5,1000) 
plot_results_table(data,200,10,1000) 
plot_results_table(data,400,1,1000) 
plot_results_table(data,400,10,1000) 
plot_results_table(data,800,1,1000)

column_headers = data.pop(0)
row_headers = [x.pop(0) for x in data]

plt.rcParams["figure.autolayout"] = True      
fig, ax = plt.subplots(figsize=(12,8),nrows=1,ncols=1)    
      
the_table = ax.table(cellText=data,
                      rowLabels=row_headers,
                      rowLoc='center',
                      colLabels=column_headers,
                      loc='center', cellLoc='center')
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
plt.box(on=None)
the_table.auto_set_font_size(False)
the_table.set_fontsize(9)
the_table.scale(1, 1.2)
the_table.auto_set_column_width(col=list(range(len(column_headers))))
plt.tight_layout()
fig.savefig('jpo_results_table.png', dpi=300 )




##########################################
########### FIGURE 8: Transport parameters
##########################################
nrows = 1; ncols = 2
fig, ax = plt.subplots(figsize=(2.5*ncols,2.5*nrows),nrows=nrows,ncols=ncols)
ax[0].plot(1e3,1e1,'-',linewidth=2,color='white',label=r'$Pr$, $Le$') #,alpha=0.85
ax[1].plot(1e1,1e0,'-',linewidth=2,color='white',label=r'$Re_{\tau}$, $Pr$')

def transport_parameters_with_Re(fig,ax,Re,Pr,Le,flag=None):   
  Sc = Pr*Le
  t, iss, qTz1_normalized, qSz1_normalized, NuT, NuS, CD, Tz1, Sz1, u0, dzTz1, dzSz1 = fluxes_history(Re,Pr,Le,flag=None)
  ax[0].loglog(Re,NuS/Sc, ls='None', marker='o', markerfacecolor="None", color=expcolor(Re,Pr,Le), label=r'%i, %i'%(Pr,Le))
  ax[0].loglog(Re,NuT/Pr,linestyle=None, marker='d',markerfacecolor="None", color=expcolor(Re,Pr,Le))
  #
  t, iss, dzTz1m_from_dzTDz1m, dzSz1m_from_dzTDz1m, dzTz1m_from_TD_bulk, dzSz1m_from_TD_bulk, dzTDz1_std_low, dzTDz1_std_high, u0m_from_TD, NuTD = fluxes_from_TD(Re,Pr,Le)
  #
  Tz1m_in_TD = 0
  Sz1m_in_TD = 0
  NuT_from_dzTDz1m = -(1/Pr)*dzTz1m_from_dzTDz1m/((1/Pr)*(1-Tz1m_in_TD))
  NuS_from_dzTDz1m = -(1/Sc)*dzTz1m_from_dzTDz1m/((1/Sc)*(1-Sz1m_in_TD)) ## !!! CAN'T BE COMPUTED !!!
  #
  ax[0].loglog(Re,NuT/Pr,linestyle=None, marker='x', markersize=8, color=expcolor(Re,Pr,Le), zorder=1) # markerfacecolor
  return fig, ax
  
def transport_parameters_with_Le(fig,ax,Re,Pr,Le,flag=None):   
  if Re==200: s='o'
  elif Re==400: s='d'
  elif Re==800: s='s'
  s='o'
  t, iss, qTz1_normalized, qSz1_normalized, NuT, NuS, CD, Tz1, Sz1, u0, dzTz1, dzSz1 = fluxes_history(Re,Pr,Le,flag=None)
  ax[1].loglog(Le,NuS/NuT/Le,linestyle='None', marker=s,markerfacecolor="None", color=expcolor(Re,Pr,Le), label=r'%i, %i'%(Re,Pr))
  print(Re,Pr,Le,np.std(-1/Pr*dzTz1[iss:]))
  return fig, ax
  
for Re in [200,400,800]:
  transport_parameters_with_Re(fig,ax,Re,1,1)
for Re in [200,400,800]:
  transport_parameters_with_Re(fig,ax,Re,1,10)
transport_parameters_with_Re(fig,ax,200,1,30)
transport_parameters_with_Re(fig,ax,400,10,10)
#transport_parameters_with_Re(fig,ax,200,1,100)

transport_parameters_with_Le(fig,ax,200,1,1)
transport_parameters_with_Le(fig,ax,200,1,10)
transport_parameters_with_Le(fig,ax,200,1,30)
transport_parameters_with_Le(fig,ax,200,1,100)
#transport_parameters_with_Le(fig,ax,400,1,1)
transport_parameters_with_Le(fig,ax,400,1,10)
#transport_parameters_with_Le(fig,ax,800,1,1)
transport_parameters_with_Le(fig,ax,800,1,10)
transport_parameters_with_Le(fig,ax,400,10,10)

ax[0].plot([300,600],[3.5,3.5*2**(1)],':k')
ax[0].plot([300,600],[10,10*2**(2/3)],'--k')
ax[0].set_xlim((100,3000))
ax[0].set_ylim((0.5,30))

ax[1].plot([1,10],[1,1*10**(-2/3)],'--k')
ax[1].plot([10,100],[0.2,0.2*10**(-1/2)],':k')

ax[0].legend(frameon=False,fontsize='x-small',loc='lower right',ncol=1,handletextpad=0.1)
ax[0].set_xlabel(r'$Re_{\tau}$')
ax[0].set_ylabel(r'$Nu_S$')
ax[0].set_ylabel(r'$\frac{Nu_T}{Pr}$ ($ \diamond $)'+'\n'+r'$\frac{Nu_S}{Sc}$ ($\circ$)', rotation=0, labelpad=25, y=0.3, fontsize='large', ha="left")
#ax[0].set_ylabel(r'$\frac{Nu_S}{Sc}$ ($\circ$)\n$\frac{Nu_T}{Pr}$ (\Box)', rotation=0, labelpad=10, y=0.45, fontsize='large', ha="right")
ax[1].set_ylim((0.02,1.2))
ax[1].legend(frameon=False,fontsize='x-small',loc='lower left',ncol=1,handletextpad=0.1)
ax[1].set_xlabel(r'$Le$')
ax[1].set_ylabel(r'$\frac{Nu_S}{Nu_TLe}$', rotation=0, labelpad=1, y=0.45, fontsize='large')

fig.tight_layout(pad=0.1,w_pad=0.25)
fig.savefig('jpofig8bis.eps',bbox_inches='tight')
fig.savefig('jpofig8bis.png',dpi=300,bbox_inches='tight')













