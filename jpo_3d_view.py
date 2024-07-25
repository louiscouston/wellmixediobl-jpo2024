"""
Plot scalars.

Usage:
    jpo_3d_view.py <files>... [--figname=<figname>]

Options:
    --figname=<figname>  figure name [default: jpo_snapshots_april2024]
"""

# python3 jpo_3d_view.py last_200_1_10_1.npz

import numpy as np
import matplotlib
import matplotlib as mpl
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
matplotlib.use('Agg') #backend plotting only--doesn't require visual soft
import matplotlib.pyplot as plt
import h5py, re, pathlib
from docopt import docopt
args = docopt(__doc__)
from support_mixing import *
from jpo_figures_support import *

figname = args['--figname']
files = args['<files>']
files.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))

for thisfile in files:
  fig = plt.figure(figsize=(6,4),tight_layout=True)
  ax = fig.add_subplot(111, projection='3d')
  regex = re.compile(r'\d+')
  Re, At, Pr, Le, TD1000 = [int(x) for x in regex.findall(thisfile)]
  data = np.load(thisfile)
  t = data['t']
  dt = data['dt']
  x = data['x']
  y = data['y']
  z = data['z']
  T = data['T']
  S = data['S']*1e7
  u = data['u']
  w = u[2]/Re
  qT = -(T[:,:,-1]-T[:,:,-2])/(z[-1]-z[-2])/(T[:,:,0].mean()-T[:,:,-1].mean())
  
  X, Y, Z = np.meshgrid(x,y,z)
  
  
  # Plot options
  kwS = {'vmin':-2,'vmax':0,'levels':np.linspace(-2,0,21),'cmap':cm.ocean_r,'extend':'both'}
  kwqT = {'vmin':4,'vmax':40,'levels': np.linspace(0,40,21),'cmap':cm.Reds,'extend':'both'}
  kww = {'vmin':-3,'vmax':3,'levels':np.linspace(-3,3,21),'cmap':cm.gnuplot2,'extend':'both'}

  # Plot contour surfaces
  CM = ax.contourf(X[:,:,-1],Y[:,:,-1],qT,zdir='z',offset=1,**kwqT)
  CS = ax.contourf(X[0,:,:],S[0,:,:],Z[0,:,:],zdir='y',offset=0,**kwS)
  CW = ax.contourf(w[:,0,:],Y[:,0,:],Z[:,0,:],zdir='x',offset=0,**kww) 

  # Set limits of the plot from coord limits
  xmin, xmax = X.min(), X.max()
  ymin, ymax = Y.min(), Y.max()
  zmin, zmax = Z.min(), Z.max()
  ax.set(xlim=[xmin, 2], ylim=[ymin, 1], zlim=[zmin, zmax])

  # Get rid of the panes
  ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
  ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
  ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

  # Plot edges
  edges_kw = dict(color='0.4', linewidth=1, zorder=1e3)
  #ax.plot([xmin, 2], [ymin, 1], [zmin, zmax], **edges_kw)
  #ax.plot([xmin, 2], [ymin, ymin], 1, **edges_kw)
  #ax.plot([xmin, xmin], [ymin, 1], 1, **edges_kw)
  ax.plot([xmin, xmin], [ymin, ymin], [zmin, zmax], **edges_kw)
  ax.plot([xmin, xmax], [ymin, ymin], 1, **edges_kw)
  ax.plot([xmin, xmin], [ymin, ymax], 1, **edges_kw)
  
  # Set labels and zticks
  #if Re==200:
  #  ax.set_yticks([0,0.2,0.4,0.6,0.8])
  #elif Re==400:
  #  ax.set_yticks([0,0.25,0.5])
  #else:
  #  ax.set_yticks([0,0.25,0.5])
  ax.set_xticks([0,0.5,1.0,1.5,2.0])
  ax.set_yticks([0,0.2,0.4,0.6,0.8,1.0])
  #if Re==200:
  ax.set_zticks([0,0.2,0.4,0.6,0.8,1.0])
  #else:
  #  ax.set_zticklabels([])
  #,ylabel='y',zlabel='z',labelpad=10)

  if Re==200:
    cbar=fig.colorbar(CW, ax=ax, fraction=0.03, pad=-0.08, aspect=30, location='top')
    cbar.ax.set_title(r'$\frac{\tilde{w}}{Re_{\tau}}$',fontsize='x-large',y=4.5)
  elif Re==400:
    cbar=fig.colorbar(CM, ax=ax, fraction=0.03, pad=-0.08, aspect=30, location='top')
    cbar.ax.set_title(r'$\tilde{q}_T$')
  elif Re==800:
    cbar=fig.colorbar(CS, ax=ax, fraction=0.03, pad=-0.08, aspect=30, location='top')
    cbar.ax.set_title(r'$\tilde{S}\times 10^7$')

  # Set zoom and angle view
  ax.view_init(30,220,0)
  #ax.set_box_aspect(None, zoom=0.9)
  ax.set_aspect('equal')

  ax.xaxis.set_rotate_label(False)
  ax.yaxis.set_rotate_label(False)
  ax.zaxis.set_rotate_label(False)
  ax.set_xlabel(r'$\tilde{x}$',fontsize='large',labelpad=10)
  ax.set_ylabel(r'$\tilde{y}$',fontsize='large',labelpad=5)
  #if Re==200: 
  ax.set_zlabel(r'$\tilde{z}$',fontsize='large',labelpad=1,y=0.2)
  
  # Show Figure
  fig.tight_layout(pad=2)
  #if Re==200:
  fig.savefig(figname+'_%i_%i_%i.png'%(Re,Pr,Le),dpi=300,bbox_inches='tight',pad_inches=0.2)#
  #else:
  #  fig.savefig(figname+'_%i_%i_%i.png'%(Re,Pr,Le),dpi=300,bbox_inches='tight',pad_inches=0)#
  
  
# Colorbars
figcolorbars = plt.figure(figsize=(7,3.5))
axcolorbars = figcolorbars.add_subplot(111)
#CW = axcolorbars.contourf(w[:,0,:],Y[:,0,:],Z[:,0,:],zdir='x',offset=0,**kww) 
#axcolorbars.set_aspect('equal')
cbar=figcolorbars.colorbar(CW, ax=axcolorbars, fraction=0.04, pad=0, aspect=30, location='top')
cbar.ax.set_title(r'$\frac{\tilde{w}}{Re_{\tau}}$',fontsize='x-large',y=5)
axcolorbars.remove()
figcolorbars.tight_layout(pad=2)
figcolorbars.savefig(figname+'_colorbar1.png',dpi=300,bbox_inches='tight')
#
figcolorbars = plt.figure(figsize=(7,3.5))
axcolorbars = figcolorbars.add_subplot(111)
CM = axcolorbars.contourf(X[:,:,-1],Y[:,:,-1],qT,zdir='z',offset=1,**kwqT)
axcolorbars.set_aspect('equal')
cbar=figcolorbars.colorbar(CM, ax=axcolorbars, fraction=0.04, pad=0, aspect=30, location='top')
cbar.ax.set_title(r'$\tilde{q}_T$')
axcolorbars.remove()
figcolorbars.tight_layout(pad=2)
figcolorbars.savefig(figname+'_colorbar2.png',dpi=300,bbox_inches='tight')
#
figcolorbars = plt.figure(figsize=(7,3.5))
axcolorbars = figcolorbars.add_subplot(111)
axcolorbars.set_aspect('equal')
CS = axcolorbars.contourf(X[0,:,:],S[0,:,:],Z[0,:,:],zdir='y',offset=0,**kwS)
cbar=figcolorbars.colorbar(CS, ax=axcolorbars, fraction=0.04, pad=0, aspect=30, location='top')
cbar.ax.set_title(r'$\tilde{S}\times 10^7$')
axcolorbars.remove()
figcolorbars.tight_layout(pad=2)
figcolorbars.savefig(figname+'_colorbar3.png',dpi=300,bbox_inches='tight')








